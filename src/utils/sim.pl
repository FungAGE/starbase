#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long qw(:config auto_abbrev no_ignore_case);
use Sort::Naturally;
use Storable;
use FileHandle;
$|=1;
# Fishtank lib
use Utils qw(dim_2_hash Parse_known_elements Commandline_check Fasta_hash_many_files Glofish_bed_hash dim_0_hash Open_FH);

sub usage {
	my $message = shift;
	my $usage = qq/
usage: starfish sim [args]

calculate k-mer similarity between captains or elements.
 
Required:
-m, --mode          STR    what features to compare, either 'cap' or 'element'.
-t, --type          STR    sequence type to compare, either 'nucl' or 'prot'.
-b, --bed           FILE   BED file with feature coordinates of mobile elements.
-x, --prefix        STR    prefix for naming all output files.
-o, --outdir        DIR    output directory.

Required, if --type nucl:
-a, --assembly      FILE   2 column tsv: genomeID, path to assembly FASTA.

Required, if --type prot:
-p, --protein       FILE   2 column tsv: genomeID, path to protein FASTA.

Required, with defaults:
--kmer              STR    k-mer size for sourmash.
                           (default: 510 if 'nucl', 17 if 'prot')
--scaled            STR    scaled parameter for sourmash.
                           (default: 100 if 'nucl', 20 if 'prot')

Optional:
-d, --dereplicated  FILE   dereplicated.txt file.
                           (output by starfish dereplicate)
-g, --group         FILE   MCL-formatted group file with elementIDs.
-F, --feat          FILE   elements.feat file.
                           (output by starfish summarize; if -m element)
-r, --restrict      FILE   file where each row contains tab-separated regionIDs.
-h, --help                 print more details and exit.

/;
	if (not defined $message) {
		$message = $usage;
	} else {
		$message = "$message\nuse -h for more details\n\n" ;
	}	
	die($message);
}

main: {

	# parse and check options
	my %opts;
	GetOptions(\%opts, 
		'mode|m=s',
		'type|t=s',
		'assembly|a=s',
		'protein|p=s',
		'bed|b=s',
		'prefix|x=s',
		'outdir|o=s',
		'restrict|r=s',
		'feat|F=s',
		'group|g=s',
		'kmer=s',
		'scaled=s',
		'dereplicated|d=s',
		'h|help');
	Opts_check(\%opts);

	# sourmash info
	# from https://sourmash.readthedocs.io/en/latest/using-sourmash-a-guide.html: 
	# The general rule is that longer k-mer sizes are less prone to false positives.
	# check dependencies
	my @commandlines = ("sourmash");
	Commandline_check(\@commandlines);

	print "\nKey parameters:
--kmer     $opts{'kmer'}
--scaled   $opts{'scaled'}\n";

	#######################################
	#### FORMAT CHECK AND READ IN DATA ####
	#######################################

		my $datestring = localtime();
		print "\n[$datestring] reading in data..\n";
	
		# define the genetic units being compared
		my $comparisonUnit = 'elements';
		$comparisonUnit = 'element captains' if ($opts{'mode'} eq 'cap');

		# Load up all assembly or protein sequences into a single hash
		# Structured: {sequenceID} = sequence
		my $sequencePaths;
		if (defined $opts{'assembly'}) {
			($sequencePaths) = dim_0_hash($opts{'assembly'}, "\t", "1");
		} elsif (defined $opts{'protein'}) {
			($sequencePaths) = dim_0_hash($opts{'protein'}, "\t", "1");
		}
		my ($id2sequence) = Fasta_hash_many_files($sequencePaths);
		
		my ($starshipFeatures, $starship2nodeData, $unitCoords, $unitCount, $restrictedComparisons, $nested);
		if (($opts{'mode'} eq 'cap') || ($opts{'mode'} eq 'element')) {

			# Parse info from Starship regions
			# structured: {contigID}{regionID}{featureID} = [begin, end, strand, tag, annotation]
			($starshipFeatures) = Glofish_bed_hash($opts{'bed'});
			
			# parse family info, if provided
			# structured: {starshipID}{family} = familyID
			($starship2nodeData) = Parse_group_file_by_member($opts{'group'}) if (defined $opts{'group'});

			# If dereplicated file provided, filter out any ship that is not a reference element at a given homologous region
			# restructure starship2nodeData, including number of represented elements and element family and region
			if (defined $opts{'dereplicated'}) {
				my ($ref2region, $ref2count, $ref2fam) = Parse_dereplicated_file_by_member($opts{'dereplicated'});
				($starshipFeatures, $starship2nodeData) = Filter_starships_by_reference($starshipFeatures, $ref2region, $ref2count, $ref2fam);
			}
			
			# Parse the coordinates of any elements with 'cap' DUF3435 genes
			# elementBegin and elementEnd are the upmost and downmost coordinates associated with this captain
			# structured: unitCoords{contigID}{regionID}{capID} = [elementBegin, elementEnd, strand, tag, ann]
			($unitCoords, $unitCount) = Parse_known_elements($starshipFeatures);
	
			# restricted (ie allowed) starship-starship comparisons
			# structured: {ship1}{ship2} = 1
			my ($restrictedComparisons) = dim_2_hash($opts{'restrict'}, "\t", "0:1") if (defined $opts{'restrict'});
		
			# for each element, identify all nested elements it contains
			# structured: {focalElementID}{nestedElementID} = 1
			my ($nested) = Parse_nested_elements($opts{'feat'}) if (defined $opts{'feat'});

		} 
				
	###################################################
	#### PARSE SEQUENCE ACCORDING TO MODE AND TYPE ####
	###################################################
		
		my $simOutfile = "$opts{'outdir'}/$opts{'prefix'}.$opts{'mode'}.$opts{'type'}.sim";
		my $nodeOutfile = "$opts{'outdir'}/$opts{'prefix'}.$opts{'mode'}.$opts{'type'}.nodes";

		if (! -f $simOutfile) {
			$datestring = localtime();
			print "[$datestring] parsing $opts{'type'} sequence from $unitCount $comparisonUnit..\n";

			# structured: {regionID/emptySiteID}{seqID} = sequence
			# structure is a little unnecessarily complicated (don't need two levels to hash) to conform with other possible parsing methods that have since been deprecated, but could be brought back
			my $unitSequence;
			if ($opts{'mode'} eq 'cap') {	
				($unitSequence) = Parse_captain_sequence($starshipFeatures, $id2sequence, $opts{'type'}); # if just comparing captains, no need to figure out nesting
			} elsif ($opts{'mode'} eq 'element') {
				($unitSequence) = Parse_element_sequence($unitCoords, $id2sequence, $nested);
			} 
	
		###################################
		#### CALCULATE MASH SIMILARITY ####
		###################################

			$datestring = localtime();
			print "[$datestring] calculating pairwise mash similarities based on $opts{'mode'} comparisons between $comparisonUnit..\n";
	
			# calculate all possible pairwise similarities
			# filter similarities by comparisons in --restrict, if provided
			# include all possible comparisons, even those that are 0
			my ($similarities) = Calculate_pairwise_mashsim($unitSequence, $restrictedComparisons, $opts{'outdir'}, $opts{'type'}, $opts{'mode'}, $opts{'prefix'}, $opts{'flanking'}, $opts{'kmer'}, $opts{'scaled'});
		
		#######################
		#### PRINT RESULTS ####
		#######################
			
			$datestring = localtime();
			print "[$datestring] printing pairwise mash similarity scores..\n";
			Print_sim_file($similarities, $simOutfile);
			
			if (defined $starship2nodeData) {
				$datestring = localtime();
				print "[$datestring] printing node data..\n";
				Print_node_file($starship2nodeData, $nodeOutfile);
			}
		
		} else {
			$datestring = localtime();
			print "[$datestring] $simOutfile already exists, so no new results to print\n";
		} 

		$datestring = localtime();
		print "[$datestring] done\n";
	
}

sub Print_node_file {
	my ($starship2nodeData, $nodeFile) = @_;
	my ($OUT) = Open_FH($nodeFile);
	
	# retrieve all column headers
	my %columns;
	foreach my $refID (keys %{$starship2nodeData}) {
		foreach my $header (keys %{$starship2nodeData->{$refID}}) {
			$columns{$header} = 1;
		}
	}
	print $OUT "id";
	foreach my $header (nsort keys %columns) {
		print $OUT "\t$header";
	}
	print $OUT "\n";
	
	foreach my $elementID (nsort keys %{$starship2nodeData}) {
		print $OUT "$elementID";
		foreach my $header (nsort keys %columns) {
			if (defined $starship2nodeData->{$elementID}->{$header}) {
				print $OUT "\t$starship2nodeData->{$elementID}->{$header}";
			} else {
				print $OUT "\tNA";
			}
		}
		print $OUT "\n";
	}		
}

sub Print_sim_file {
	my ($similarities, $outfile) = @_;
	my ($OUT) = Open_FH($outfile);
	#print $OUT "from\tto\tweight\n";
	foreach my $refID (nsort keys %{$similarities}) {
		foreach my $queID (nsort keys %{$similarities->{$refID}}) {
			print $OUT "$refID\t$queID\t".sprintf("%.3f", $similarities->{$refID}->{$queID} / 2)."\n";
		}
	}	
}

sub Calculate_pairwise_mashsim {
	my ($unitSequence, $restrictedComparisons, $OUTDIR, $TYPE, $MODE, $PREFIX, $FLANKING, $KMER, $SCALED) = @_;

	# if checkpoint file exists, read it in and skip de novo calculations
	my $datestring = localtime();
	my %similarities;	
	my $matFile = "$OUTDIR/$PREFIX.$MODE.$TYPE.mat";
	
	if (! -f $matFile) {
		# print out a file with all requested regionIDs sequences
		my $seqFile = "$OUTDIR/mashRef.$PREFIX.$MODE.$TYPE.fa";
		my ($seqOUT) = Open_FH($seqFile);

		# since --mode is restricted to 'cap' and 'element', each regionID will at most be associated with 1 sequence feature
		# this means we can build a sequence database, then calculate all sketches for all sequences
		foreach my $regionID (keys %{$unitSequence}) {
			foreach my $featureID (keys %{$unitSequence->{$regionID}}) {
				print $seqOUT ">$regionID\n$unitSequence->{$regionID}->{$featureID}\n";
			}
		}

		# build a kmer signature sketch of all input sequences
		my $sketchFile = "$OUTDIR/mashRef.$PREFIX.$MODE.$TYPE.sig";
		if ($TYPE eq 'nucl') {
			my ($failcheck) = system("sourmash sketch dna --singleton --output $sketchFile -p k=$KMER,scaled=$SCALED,noabund $seqFile 2>> $OUTDIR/mash.$PREFIX.$MODE.$TYPE.err");
			$datestring = localtime();					
			if ($failcheck != 0) { die "\n\n[$datestring] error: could not execute 'sourmash sketch dna' on commandline, exiting..\n$!\n";}
		} elsif ($TYPE eq 'prot') {
			my ($failcheck) = system("sourmash sketch protein --singleton --output $sketchFile -p k=$KMER,scaled=$SCALED,noabund $seqFile 2>> $OUTDIR/mash.$PREFIX.$MODE.$TYPE.err");
			$datestring = localtime();					
			if ($failcheck != 0) { die "\n\n[$datestring] error: could not execute 'sourmash sketch prot' on commandline, exiting..\n$!\n";}
		}
				
		# execute sourmash compute
		$datestring = localtime();
		print "[$datestring] ";
		if ($TYPE eq 'nucl') {
			my ($failcheck) = system("sourmash compare --dna --csv $matFile -k $KMER $sketchFile 2>> $OUTDIR/mash.$PREFIX.$MODE.$TYPE.err");
							
			if ($failcheck != 0) { die "\n\n[$datestring] error: could not execute 'sourmash compare' on commandline, exiting..\n$!\n";}
		} elsif ($TYPE eq 'prot') {
			my ($failcheck) = system("sourmash compare --protein --csv $matFile -k $KMER $sketchFile 2>> $OUTDIR/mash.$PREFIX.$MODE.$TYPE.err");
			if ($failcheck != 0) { die "\n\n[$datestring] error: could not execute 'sourmash compare' on commandline, exiting..\n$!\n";}
		}
		# clean up
		system("rm $sketchFile $seqFile $OUTDIR/mash.$PREFIX.$MODE.$TYPE.err");
	} else {
		print "[$datestring] $matFile exists, skipping de novo mash calculations..\n";
	}

	# read in jaccard similarity matrix
	my %observedComparisons;
	open(my $matIN, '<', $matFile) or usage("\nError: could not open $matFile for reading\n");
	my $header = <$matIN>;
	chomp $header;
	$header =~ s/\R//g;
	my @queRegionIDs = split/,/, $header;
	my @refRegionIDs = @queRegionIDs;
	while (my $line = <$matIN>) {
		chomp $line;
		$line =~ s/\R//g;
		my @similarities = split/,/, $line;
		my $refRegionID = shift @refRegionIDs;
		chomp $refRegionID;
		$refRegionID =~ s/\R//g;
		foreach my $queRegionID (@queRegionIDs) {
			chomp $queRegionID;
			$queRegionID =~ s/\R//g;
			my $sim = shift @similarities;
			next if (($refRegionID eq $queRegionID) || (exists $restrictedComparisons->{$refRegionID}->{$queRegionID}) || (exists $restrictedComparisons->{$queRegionID}->{$refRegionID}));
			next if ($sim == 0.0);
			# sort regionIDs so we can eventually take the average of each pairwise comparison
			my ($ref, $que) = nsort($refRegionID,$queRegionID);
			$similarities{$ref}{$que} += $sim;
			$observedComparisons{$ref}{$que} = 1;
			$observedComparisons{$que}{$ref} = 1;
		}
	} 
	
	# check that all comparisons have been made, even those that are 0
	foreach my $regionID1 (keys %{$unitSequence}) {
		foreach my $regionID2 (keys %{$unitSequence}) {
			next if ($regionID1 eq $regionID2);
			if (not exists $observedComparisons{$regionID1}{$regionID2}) {
				$similarities{$regionID1}{$regionID2} = 0.000;
			}
		}
	}
	
	return(\%similarities);
}

sub Parse_element_sequence {
	my ($unitCoords, $id2sequence, $nested) = @_;

	my %unitSequence;
	foreach my $contigID (keys %{$unitCoords}) {
		foreach my $regionID (keys %{$unitCoords->{$contigID}}) {
		
			# figure out the coordinates of all nested elements so we can ignore them
			my @forbiddenRanges;
			if (defined $nested) {
				if (exists $nested->{$regionID}) {
					foreach my $nestedID (keys %{$nested->{$regionID}}) {
						foreach my $nestedCapID (keys %{$unitCoords->{$contigID}->{$nestedID}}) { # notice we iterate directly through nestedID, since we know it must be on the same contig
							my ($nestedBegin, $nestedEnd, $nestedStrand) = @{$unitCoords->{$contigID}->{$nestedID}->{$nestedCapID}};
							push @forbiddenRanges, "${nestedBegin}-${nestedEnd}";
						}
					}
				}
			}
			
			# consolidate all overlapping ranges into non-overlapping segments
			my @consolidatedRanges;
			if (scalar @forbiddenRanges > 0) {
				@forbiddenRanges = nsort @forbiddenRanges; # should arrange all coordinates nicely
				my $currentRange = shift @forbiddenRanges;
				my ($currentBegin, $currentEnd) = split/-/, $currentRange;
				
				# attempt to extend, or find next non-overlapping range
				while (scalar @forbiddenRanges > 0) {
					my ($nextRange) = shift @forbiddenRanges;
					my ($nextBegin, $nextEnd) = split/-/, $nextRange;
					if (($nextEnd > $currentEnd) && ($nextBegin <= $currentEnd)) { #
						# extend currentEnd if the next range overlaps partially, yet extends further
						$currentEnd = $nextEnd; 
					} elsif ($nextBegin > $currentEnd) {
						# we have a non-overlapping segment, so store the previous segment and continue on from this one
						push @consolidatedRanges, $currentBegin, $currentEnd;
						($currentBegin, $currentEnd) = ($nextBegin, $nextEnd);
					} else {
						# any segment that does not meet these criteria is completely enveloped, so we can ignore it
						next;
					}
				}
				# add the last observed currentBegin and currentEnd
				push @consolidatedRanges, $currentBegin, $currentEnd;
			}
			
			# now iterate through the focal starship of interest
			foreach my $capID (keys %{$unitCoords->{$contigID}->{$regionID}}) {
				my ($shipBegin, $shipEnd, $shipStrand) = @{$unitCoords->{$contigID}->{$regionID}->{$capID}};
				
				# take all starship sequence that is not contained within forbidden coords
				if (exists $id2sequence->{$contigID}) {
					
					if (scalar @consolidatedRanges > 0) {
						# troubleshooting
						#print "$regionID\t";
						#print Dumper(\@consolidatedRanges);
						# add begin and end coordinates so array now contains ordered pairs of coordinates that we should be taking sequence from
						unshift @consolidatedRanges, $shipBegin;
						push @consolidatedRanges, $shipEnd;
						
						# store all segments as individual sequences
						while (scalar @consolidatedRanges > 1) {
							my $segmentBegin = shift @consolidatedRanges;
							my $segmentEnd = shift @consolidatedRanges;
							
							# sanity check
							if (($segmentEnd - $segmentBegin) > 31) { # will also avoid situations where the boundaries were extended past the end of the focal element
								my $segSeq = substr($id2sequence->{$contigID}, $segmentBegin - 1, $segmentEnd - $segmentBegin);
								if (exists $unitSequence{$regionID}{$regionID}) {
									$unitSequence{$regionID}{$regionID} .= "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN$segSeq"; # add 40 Ns to simulate the masking of intervening sequence
								} else {
									$unitSequence{$regionID}{$regionID} = $segSeq; 
								}
							}
						}
						
						
					} else {
						my $seq = substr($id2sequence->{$contigID}, $shipBegin - 1, $shipEnd - $shipBegin);
						$unitSequence{$regionID}{$regionID} = $seq;
					}
					
				} else {
					my $datestring = localtime();
					print "[$datestring] warning: can't find sequence for $contigID in supplied assembly files, skipping..\n";					
				}
			}
		}
	}
	return(\%unitSequence);
}

sub Parse_captain_sequence {
	my ($starshipFeatures, $id2sequence, $TYPE) = @_;

	my %unitSequence;
	foreach my $contigID (keys %{$starshipFeatures}) {
		foreach my $regionID (keys %{$starshipFeatures->{$contigID}}) {
			foreach my $featureID (keys %{$starshipFeatures->{$contigID}->{$regionID}}) {
				my ($begin, $end, $strand, $tag, $ann) = @{$starshipFeatures->{$contigID}->{$regionID}->{$featureID}};
				next if ($tag !~ m/^cap$/); # only look at captain sequences
				if ($TYPE eq 'nucl') { # assume nucleotide id2seq
					if (exists $id2sequence->{$contigID}) {
						my $seq = substr($id2sequence->{$contigID}, $begin - 1, $end - $begin);
						$unitSequence{$regionID}{$featureID} = $seq;
					} else {
						my $datestring = localtime();
						print "[$datestring] warning: can't find $TYPE sequence for $contigID in supplied assembly files, skipping..\n";					
					}
				} elsif ($TYPE eq 'prot') {  # assume amino acid id2seq
					if (exists $id2sequence->{$featureID}) {
						$unitSequence{$regionID}{$featureID} = $id2sequence->{$featureID};
					} else {
						my $datestring = localtime();
						print "[$datestring] warning: can't find $TYPE sequence for $featureID in supplied protein files, skipping..\n";					
					}
				}
			}
		}
	}
	return(\%unitSequence);
}

sub Filter_starships_by_reference {
	my ($starshipFeaturesTemp, $ref2region, $ref2count, $ref2group) = @_;
	
	my (%starshipFeatures, %starship2nodeData);
	# structured: {contigID}{regionID}{featureID} = [begin, end, strand, tag, annotation]
	
	foreach my $contigID (keys %{$starshipFeaturesTemp}) {
		foreach my $starshipID (keys %{$starshipFeaturesTemp->{$contigID}}) {
			if (exists $ref2region->{$starshipID}) { # only keep starships that are references
				foreach my $featureID (keys %{$starshipFeaturesTemp->{$contigID}->{$starshipID}}) {
					@{$starshipFeatures{$contigID}{$starshipID}{$featureID}} = @{$starshipFeaturesTemp->{$contigID}->{$starshipID}->{$featureID}};
					$starship2nodeData{$starshipID}{'regionID'} = $ref2region->{$starshipID};
					$starship2nodeData{$starshipID}{'replicateCount'} = $ref2count->{$starshipID};
					$starship2nodeData{$starshipID}{'groupID'} = $ref2group->{$starshipID};
				}
			}
		}
	}
	return(\%starshipFeatures, \%starship2nodeData);	
}

sub Parse_dereplicated_file_by_member {
	my ($dereplicatedFile) = @_;
	my $datestring = localtime();					
	my (%element2region, %elementRegionCount, %element2groupID);
	open (my $IN, '<', $dereplicatedFile) or usage("\n\n[$datestring] error: cannot read $dereplicatedFile, exiting\n");
	while (my $line = <$IN>) {
		chomp $line;
		my ($regionID, $elementGroupID, $refElementID, $refSiteID, $otherElementIDsString) = split/\t/, $line;
		$element2region{$refElementID} = $regionID;
		$elementRegionCount{$refElementID}++;
		$element2groupID{$refElementID} = $elementGroupID;
		next if ($otherElementIDsString eq '.');
		my @otherElements = split/,/, $otherElementIDsString;
		foreach my $otherElementID (@otherElements) {
			$elementRegionCount{$refElementID}++;
		}
	}
	return(\%element2region, \%elementRegionCount, \%element2groupID);	
}

sub Parse_group_file_by_member {
	my ($clusteringOutfile) = @_;
	my $datestring = localtime();					
	my %element2group;
	open (my $IN, '<', $clusteringOutfile) or usage("\n\n[$datestring] error: cannot read $clusteringOutfile, exiting\n");
	while (my $line = <$IN>) {
		chomp $line;
		my (@elements) = split/\s+/, $line;
		my $group = shift @elements;
		$group =~ s/:$//;
		
		foreach my $element (@elements) {
			$element2group{$element} = $group;
		}
	}
	return(\%element2group);	
}


sub Parse_nested_elements {
	my ($featureFile) = @_;
	my %nested;
	open(my $IN, '<', $featureFile) or usage("\nError: could not open $featureFile for reading\n");
	while (my $line = <$IN>) {
		next if ($line =~ m/^#/);
		chomp $line;
		#acquire nesting info
		my ($contigID, $starshipID, $captainID, $elementBegin, $elementEnd, $elementLength, $strand, $boundaryType, $emptyContig, $emptyBegin, $emptyEnd, $emptySeq, $upTSD, $downTSD, $TSDedit, $upTIR, $downTIR, $TIRedit, $nestedInside, $containNested) = split("\t", $line);
		my @nested = split/,/, $containNested;
		foreach my $nest (@nested) {
			$nested{$starshipID}{$nest} = 1;
		}
	}
	return(\%nested);	
}

#################
### for recycling

# sub Parse_cargo_all_sequence {
# 	my ($starshipFeatures, $unitCoords, $id2sequence, $nested) = @_;
# 
# 	my %unitSequence;
# 	foreach my $contigID (keys %{$unitCoords}) {
# 		foreach my $regionID (keys %{$unitCoords->{$contigID}}) {
# 			
# 
# 			
# 			
# 			
# 			
# 			foreach my $capID (keys %{$unitCoords->{$contigID}->{$regionID}}) {
# 				my ($shipBegin, $shipEnd, $shipStrand) = @{$unitCoords->{$contigID}->{$regionID}->{$capID}};
# 				my ($capBegin, $capEnd, $capStrand) = @{$starshipFeatures->{$contigID}->{$regionID}->{$capID}};
# 				
# 				# take all sequence, excluding captain coordinates (make sure captain coordinates are contained within ship coordinates)
# 				if (($capBegin >= $shipBegin) && ($capEnd <= $shipEnd)) {
# 					if (exists $id2sequence->{$contigID}) {
# 						my @coordinates;
# 						push @coordinates, $shipBegin, $shipEnd, $capBegin, $capEnd;
# 						@coordinates = nsort(@coordinates);
# 						my $seq1 = substr($id2sequence->{$contigID}, $coordinates[0] - 1, $coordinates[1] - $coordinates[0]);
# 						my $seq2 = substr($id2sequence->{$contigID}, $coordinates[2] - 1, $coordinates[3] - $coordinates[2]);
# 						my $seq1header = "$regionID.cargo1";
# 						my $seq2header = "$regionID.cargo2";
# 						$unitSequence{$regionID}{$seq1header} = $seq1 if (length($seq1) > 1);
# 						$unitSequence{$regionID}{$seq2header} = $seq2 if (length($seq2) > 1);
# 					} else {
# 						my $datestring = localtime();
# 						print "[$datestring] warning: can't find sequence for $contigID in supplied assembly files, skipping..\n";					
# 					}
# 				} else {
# 					my $datestring = localtime();
# 					print "[$datestring] warning: captain $capID does not fall within the boundaries of its associated Starship $regionID, skipping cargo-all sequence retrieval\n";
# 				}
# 			}
# 		}
# 	}
# 	return(\%unitSequence);
# }


# sub Parse_cargo_gene_sequence {
# 	my ($starshipFeatures, $id2sequence, $TYPE, $nested) = @_;
# 
# 	my %unitSequence;
# 	foreach my $contigID (keys %{$starshipFeatures}) {
# 		foreach my $regionID (keys %{$starshipFeatures->{$contigID}}) {
# 
# 			# figure out the genes within nested elements so we can ignore them
# 			my %forbiddenGenes;
# 			if (defined $nested) {
# 				if (exists $nested->{$regionID}) {
# 					foreach my $nestedID (keys %{$nested->{$regionID}}) {
# 						foreach my $nestedFeatureID (keys %{$starshipFeatures->{$contigID}->{$nestedID}}) { # notice we iterate directly through nestedID, since we know it must be on the same contig
# 							$forbiddenGenes{$nestedFeatureID} = 1;
# 						}
# 					}
# 				}
# 			}
# 
# 			# now iterate through focal element, and skip any forbidden features
# 			foreach my $featureID (keys %{$starshipFeatures->{$contigID}->{$regionID}}) {
# 				
# 				next if ((scalar keys %forbiddenGenes > 0) && (exists $forbiddenGenes{$featureID}));
# 				
# 				my ($begin, $end, $strand, $tag, $ann) = @{$starshipFeatures->{$contigID}->{$regionID}->{$featureID}};
# 				next if ($tag =~ m/^cap$|^insert$|^flank$|^extend$|^align$/); # only look at cargo gene sequences
# 				if ($TYPE eq 'nucl') { # assume nucleotide id2seq
# 					if (exists $id2sequence->{$contigID}) {
# 						my $seq = substr($id2sequence->{$contigID}, $begin - 1, $end - $begin);
# 						$unitSequence{$regionID}{$featureID} = $seq;
# 					} else {
# 						my $datestring = localtime();
# 						print "[$datestring] warning: can't find $TYPE sequence for $contigID in supplied assembly files, skipping..\n";					
# 					}
# 				} elsif ($TYPE eq 'prot') {  # assume amino acid id2seq
# 					if (exists $id2sequence->{$featureID}) {
# 						$unitSequence{$regionID}{$featureID} = $id2sequence->{$featureID};
# 					} else {
# 						my $datestring = localtime();
# 						print "[$datestring] warning: can't find $TYPE sequence for $featureID in supplied protein files, skipping..\n";					
# 					}
# 				}
# 			}
# 		}
# 	}
# 	return(\%unitSequence);
# }


sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage("\nError: no arguments provided\n") if (scalar keys %{$opts} == 0);
	usage("\nError: please provide a string to --prefix\n") if (not defined $opts->{'prefix'});
	usage("\nError: please provide a directory to --outdir\n") if (not defined $opts->{'outdir'});
	usage("\nError: the directory provided to --outdir does not exist\n") if (! -d $opts->{'outdir'});
	if (defined $opts->{'restrict'}) {
		usage("\nError: the file provided to --restrict does not exist\n") if (! -f $opts->{'restrict'});
	}
	if (defined $opts->{'feat'}) {
		usage("\nError: the file provided to --feat does not exist\n") if (! -f $opts->{'feat'});
	}
	if (defined $opts->{'dereplicated'}) {
		usage("\nError: the file provided to --dereplicated does not exist\n") if (! -f $opts->{'dereplicated'});
	}
	if (not defined $opts->{'mode'}) {
		usage("\nError: you must specify an option for --mode, either 'cap' or 'element'\n");
	} elsif ($opts->{'mode'} !~ m/^cap$|^element$/) {
		usage("\nError: unrecognized option for --mode, must be either 'cap' or 'element'\n");
	} 
	usage("\nError: please provide file to --bed\n") if ((not defined $opts->{'bed'}) && ($opts->{'mode'} =~ m/^cap$|^element$/));
	if (defined $opts->{'bed'}) {
		usage("\nError: the file provided to --bed does not exist\n") if (! -f $opts->{'bed'});
	}
	if (not defined $opts->{'type'}) {
		usage("\nError: you must specify an option for --type, either 'nucl' or 'prot'\n");
	} elsif ($opts->{'type'} !~ m/^nucl$|^prot$/) {
		usage("\nError: unrecognized option for --type, must be either 'nucl' or 'prot'\n");
	} 
	# and celement are restricted to nucleotide-based comparisons (--type nucl)
	if (($opts->{'mode'} =~ m/^element$/) && ($opts->{'type'} =~ m/^prot$/)) {
		usage("\nError: --mode 'element' is restricted to nucleotide-based comparisons of --type nucl\n");
	}
	if (defined $opts->{'group'}) {
		usage("\nError: the file provided to --group does not exist\n") if (! -f $opts->{'group'});
	}
	usage("\nError: please provide a file to --assembly\n") if ((not defined $opts->{'assembly'}) && ($opts->{'type'} eq 'nucl'));
	usage("\nError: the file provided to --assembly does not exist\n") if ((defined $opts->{'assembly'}) && (! -f $opts->{'assembly'}) && ($opts->{'type'} eq 'nucl'));
	usage("\nError: please provide a file to --protein\n") if ((not defined $opts->{'protein'}) && ($opts->{'type'} eq 'prot'));
	usage("\nError: the file provided to --protein does not exist\n") if ((defined $opts->{'protein'}) && (! -f $opts->{'protein'}) && ($opts->{'type'} eq 'prot'));

	if ($opts->{'type'} eq 'nucl' && not defined $opts->{'kmer'}) {
		$opts->{'kmer'} = 510;
	} elsif ($opts->{'type'} eq 'prot' && not defined $opts->{'kmer'}) {
		$opts->{'kmer'} = 17;
	}
	if ($opts->{'type'} eq 'nucl' && not defined $opts->{'scaled'}) {
		$opts->{'scaled'} = 100;
	} elsif ($opts->{'type'} eq 'prot' && not defined $opts->{'scaled'}) {
		$opts->{'scaled'} = 20;
	}
}
