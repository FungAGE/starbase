import tempfile
import subprocess
import logging
from src.utils.blast_utils import hmmsearch, mmseqs_easy_cluster, parse_hmmer, calculate_similarities, write_similarity_file, cluster_sequences, write_cluster_files

logger = logging.getLogger(__name__)

# classifcation pipeline that should be used:
# - classify query sequences for blast page, display results
# - classify submitted sequences from submission page, assign an accession, and input into submissionsdatabase

# accession format:
# - normal ship accession: SBS123456
# - updated ship accession: SBS123456.1

########################################################
# assigning accessions
########################################################
# part 1: check if the sequence is already in the database
# part 2: if it is completely identical, almost identical, or contained within an existing sequence (also flag for review), assign the normal accession
# part 3: if it is novel, assign a new accession

# workflow:
# first check for exact matches
# then check for contained within matches
# then check for almost identical matches
# if no matches, assign a new accession

# if a sequence is a truncated version of a longer sequence, assign the longer sequence accession, flag for review

# fast method for checking for exact matches:
# md5 hash of sequence

# method for checking for contained/highly similar sequences:
# k-mers

########################################################
# classification pipeline
########################################################

def classify_sequence(sequence, db_list, threads=1):
    """Main classification pipeline that coordinates the 3-step process"""   
    # perform the following steps in order:

    # Part 1: Captain gene similarity -> family assignment
    family = classify_family(sequence, db_list, threads)
    
    # Part 2: Cargo jaccard scores -> navis assignment 
    navis = classify_navis(sequence, threads)
    
    # Part 3: Sequence similarity -> haplotype assignment
    haplotype = classify_haplotype(sequence)
    
    return family, navis, haplotype

def classify_family(sequence, db_list, input_eval=0.001, threads=1):
    """Uses hmmsearch to assign family based on captain gene similarity"""
    # first, assign all Starships to a family by searching the candidate captain sequences against the reference captain HMM profile database and identifying the best profile hit to each sequence:
    tmp_hmmer='elementFinder/macpha6_tyr_vs_YRsuperfams.out'
    hmmer_db = db_list["gene"]["captain"]["hmm"]["tyr"]
    # TODO: add this to hmmsearch method?
    hmmer_cmd = f"hmmsearch --noali --notextw --max -o {tmp_hmmer} --cpu {threads} -E {input_eval} {hmmer_db} {sequence}"
    logger.info(f"Running HMMER search: {hmmer_cmd}")
    subprocess.run(hmmer_cmd, shell=True)
    hmmer_results = hmmsearch(db_list=db_list, 
                             query_type="captain",
                             threads=threads,
                             query_fasta=sequence)
    # Parse results and return family assignment
    # perl -p -e 's/ +/\t/g' elementFinder/macpha6_tyr_vs_YRsuperfams.out | cut -f1,3,5 | grep -v '#' | sort -k3,3g | awk '!x[$1]++' > elementFinder/macpha6_tyr_vs_YRsuperfams_besthits.txt
    # TODO: make sure this parses correctly, or just use a simpler method for grabbing the family assignment from the hmmsearch output
    family_assignment = parse_hmmer(hmmer_results)
    return family_assignment

def classify_navis(sequence, threads=1):
    """Uses mmseqs2 to cluster captain sequences"""
    # Need to implement this
    # now group all Starships into naves using mmseqs2 easy-cluster on the captain sequences with a very permissive 50% percent ID/ 25% coverage threshold:
    # mmseqs easy-cluster geneFinder/macpha6_tyr.filt_intersect.fas elementFinder/macpha6_tyr elementFinder/ --threads 2 --min-seq-id 0.5 -c 0.25 --alignment-mode 3 --cov-mode 0 --cluster-reassign

    clusters = mmseqs_easy_cluster(
        sequence,
        min_seq_id=0.5,
        coverage=0.25,
        threads=threads
    )

    # clusters will be a dict like:
    # {
    #    'seq1': 'cluster1_rep',
    #    'seq2': 'cluster1_rep',
    #    'seq3': 'cluster2_rep',
    #    ...
    # }

    # parse results
    # $STARFISHDIR/aux/mmseqs2mclFormat.pl -i elementFinder/macpha6_tyr_cluster.tsv -g navis -o elementFinder/
    return clusters

def classify_haplotype(sequence):
    """Port of starfish sim/group workflow for haplotype assignment"""
    # 1. Calculate similarities (port from sim.pl)
    similarities = calculate_similarities(
        sequence_file,
        mode='element',
        seq_type='nucl',
        threads=threads
    )

    # Write results if needed
    write_similarity_file(similarities, "output.sim")
    
    # 2. Group sequences (port from group.pl)
    # After calculating similarities
    groups, node_data, edge_data = cluster_sequences(
        "similarities.sim",
        group_prefix="HAP",
        inflation=1.5,
        threshold=0.05,
        threads=threads
    )

    # Write results
    write_cluster_files(groups, node_data, edge_data, "output_prefix")    

    return haplotype_assignment

# use methods from starfish pipeline

def kmer_similarity():
    # Sequence similarity calculation using sourmash
    # Handling different sequence types (nucleotide/protein)
    # K-mer based comparison
    # Output formatting
    return None

def haplotype_grouping():
    # MCL clustering interface
    # Similarity value remapping
    # Group naming/formatting
    # Node/edge data generation for visualization
    return None

# use sourmash and mcl to group all elements into haplotypes based on pairwise k-mer similarities of element nucleotide sequences:
# starfish sim -m element -t nucl -b elementFinder/macpha6.elements.bed -x macpha6 -o elementFinder/ -a ome2assembly.txt
# starfish group -m mcl -s elementFinder/macpha6.element.nucl.sim -i hap -o elementFinder/ -t 0.05

