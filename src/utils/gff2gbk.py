#!/usr/bin/env python
"""Convert a GFF and associated FASTA file into GenBank format.

Original script written by Brad Chapman and Hsiao Yi - https://github.com/chapmanb/bcbb/blob/master/gff/Scripts/gff/gff_to_genbank.py
Edited by Kartik Chundru to crop the FASTA sequence between the first and last annotated regions.

Usage:
    gff_to_genbank.py <GFF annotation file> <FASTA sequence file>
"""
from __future__ import print_function

import sys
import tempfile
import logging
from pathlib import Path

from Bio import SeqIO
from Bio import Seq
from Bio import SeqFeature
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

logger = logging.getLogger(__name__)

def main(gff_file, fasta_file, output_file=None):
    """Convert a GFF file to GenBank format using sequence from a FASTA file"""
    try:
        # Convert to strings in case Path objects are passed
        gff_file = str(gff_file)
        fasta_file = str(fasta_file)
        
        if output_file is None:
            output_file = gff_file.rsplit(".", 1)[0] + ".gb"
        else:
            output_file = str(output_file)

        logger.info(f"Converting GFF file: {gff_file}")
        logger.info(f"Using FASTA file: {fasta_file}")
        logger.info(f"Output will be written to: {output_file}")

        # Load the sequence from the FASTA file
        try:
            with open(fasta_file) as handle:
                fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            logger.info(f"Successfully loaded {len(fasta)} sequences from FASTA file")
        except Exception as e:
            logger.error(f"Error loading FASTA file: {str(e)}")
            raise

        # Parse the GFF file and apply transformations
        try:
            with open(gff_file) as handle:
                gff = list(_check_gff(GFF.parse(handle, fasta)))
                gff = list(_fix_ncbi_id(gff))
                gff = list(_extract_regions(gff))
                gff = [_flatten_features(rec) for rec in gff]
            logger.info(f"Successfully processed {len(gff)} records from GFF file")
        except Exception as e:
            logger.error(f"Error processing GFF file: {str(e)}")
            raise

        # Write the GenBank file
        try:
            SeqIO.write(gff, output_file, "genbank")
            logger.info(f"Successfully wrote GenBank file to: {output_file}")
        except Exception as e:
            logger.error(f"Error writing GenBank file: {str(e)}")
            raise
        
        return output_file
        
    except Exception as e:
        logger.error(f"Error converting GFF to GenBank: {str(e)}")
        raise


def _fix_ncbi_id(fasta_iter):
    """GenBank identifiers can only be 16 characters; try to shorten NCBI."""
    for rec in fasta_iter:
        if len(rec.name) > 16 and rec.name.find("|") > 0:
            new_id = [x for x in rec.name.split("|") if x][-1]
            logger.warning(f"Shortening NCBI name {rec.id} to {new_id}")
            rec.id = new_id
            rec.name = new_id
        
        # Shorten feature name length
        for i in range(len(rec.features)):
            if len(rec.features[i].type) > 15:
                old_type = rec.features[i].type
                rec.features[i].type = rec.features[i].type[0:15]
                logger.debug(f"Shortened feature type from {old_type} to {rec.features[i].type}")
        yield rec


def _check_gff(gff_iterator):
    """Generator to check GFF records, replacing unknown sequences with actual sequences"""
    for rec in gff_iterator:
        # Check if sequence is empty or unknown
        if not rec.seq or str(rec.seq).count("N") == len(rec.seq):
            logger.warning(f"Record {rec.id} has unknown sequence")
            continue
            
        # Add required GenBank annotations
        rec.annotations["molecule_type"] = "DNA"
        
        yield rec


def _extract_regions(gff_iterator):
    """Extract regions from the first annotated position to the last annotated position,
    and updates the locations to correspond to the location in the sequence."""
    for rec in gff_iterator:
        try:
            if not rec.features:
                logger.warning(f"No features found in record {rec.id}")
                continue

            loc = min([i.location.start for i in rec.features])
            endloc = max([i.location.end for i in rec.features])
            logger.debug(f"Extracting region {loc}-{endloc} from record {rec.id}")

            # Update feature locations
            for i in range(len(rec.features)):
                rec.features[i].location = SeqFeature.FeatureLocation(
                    SeqFeature.ExactPosition(rec.features[i].location.start - loc),
                    SeqFeature.ExactPosition(rec.features[i].location.end - loc),
                    strand=rec.features[i].strand,
                )
                
                # Update sub-feature locations
                for j in range(len(rec.features[i].sub_features)):
                    rec.features[i].sub_features[j].location = SeqFeature.FeatureLocation(
                        SeqFeature.ExactPosition(
                            rec.features[i].sub_features[j].location.start - loc
                        ),
                        SeqFeature.ExactPosition(
                            rec.features[i].sub_features[j].location.end - loc
                        ),
                        strand=rec.features[i].sub_features[j].strand,
                    )

            # Extract sequence region
            rec.seq = rec.seq[loc:endloc]
            rec.annotations["molecule_type"] = "DNA"
            
            logger.debug(f"Successfully processed record {rec.id}")
            yield rec
            
        except Exception as e:
            logger.error(f"Error processing record {rec.id}: {str(e)}")
            continue


def _flatten_features(rec):
    """Make sub_features in an input rec flat for output.
    GenBank does not handle nested features, so we want to make everything top level.
    Also ensure proper protein translations for CDS features."""
    try:
        out = []
        for f in rec.features:
            cur = [f]
            while len(cur) > 0:
                nextf = []
                for curf in cur:
                    # Fix protein translations for CDS features
                    if curf.type == "CDS":
                        curf = _fix_cds_translation(curf, rec.seq)
                    
                    out.append(curf)
                    if len(curf.sub_features) > 0:
                        nextf.extend(curf.sub_features)
                cur = nextf
        rec.features = out
        return rec
    except Exception as e:
        logger.error(f"Error flattening features for record {rec.id}: {str(e)}")
        raise


def _fix_cds_translation(feature, sequence):
    """Ensure CDS features have correct protein translations."""
    try:
        # Extract the CDS sequence
        cds_seq = feature.extract(sequence)
        
        # Try strict translation first (with start/stop codons)
        try:
            protein_seq = cds_seq.translate(table=1, cds=True)  # Standard genetic code
        except Exception:
            # If strict translation fails, try lenient translation
            # This allows CDS features that don't have proper start/stop codons
            protein_seq = cds_seq.translate(table=1, to_stop=True)
        
        # Update the translation qualifier
        if 'translation' not in feature.qualifiers:
            feature.qualifiers['translation'] = [str(protein_seq)]
        else:
            # Replace existing translation with correct one
            feature.qualifiers['translation'] = [str(protein_seq)]
            
        logger.debug(f"Updated translation for CDS at {feature.location}")
        
    except Exception as e:
        logger.warning(f"Could not fix translation for CDS at {feature.location}: {str(e)}")
        # Keep the original feature if translation fails
        
    return feature


if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    if len(sys.argv) != 3:
        logger.error("Incorrect number of arguments")
        print("Usage: gff2genbank.py <gff_file> <fasta_file>")
        sys.exit(1)
        
    try:
        main(sys.argv[1], sys.argv[2])
    except Exception as e:
        logger.error(f"Program failed: {str(e)}")
        sys.exit(1)