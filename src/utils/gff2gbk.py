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

def parse_gff_attributes(attr_string):
    """Parse GFF3 attributes string into a dictionary.
    
    GFF3 attributes are semicolon-separated key=value pairs.
    Example: "ID=gene1;Name=ABC1;Alias=protein1"
    
    Args:
        attr_string: GFF3 attributes string
        
    Returns:
        dict: Parsed attributes
    """
    if not attr_string or attr_string == '' or str(attr_string) == 'nan':
        return {}
    
    attrs = {}
    try:
        for pair in str(attr_string).split(';'):
            if '=' in pair:
                key, value = pair.split('=', 1)
                attrs[key.strip()] = value.strip()
    except Exception as e:
        logger.debug(f"Error parsing GFF attributes '{attr_string}': {str(e)}")
    
    return attrs

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
    Also ensure proper protein translations for CDS features and extract gene labels from Alias."""
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
                    
                    # Extract Alias from qualifiers and use it as gene label
                    curf = _extract_gene_label(curf)
                    
                    out.append(curf)
                    if len(curf.sub_features) > 0:
                        nextf.extend(curf.sub_features)
                cur = nextf
        rec.features = out
        return rec
    except Exception as e:
        logger.error(f"Error flattening features for record {rec.id}: {str(e)}")
        raise


def _extract_gene_label(feature):
    """Extract gene label from GFF attributes using a fallback hierarchy.
    
    Tries to find a meaningful gene label from various sources:
    1. Alias attribute (from GFF)
    2. Name attribute
    3. ID attribute  
    4. gene qualifier
    5. locus_tag qualifier
    6. product qualifier (shortened if needed)
    
    Args:
        feature: BioPython SeqFeature object
        
    Returns:
        feature: Modified feature with gene/label qualifiers
    """
    try:
        label = None
        
        # Priority 1: Check standard BioPython qualifiers first
        if 'Alias' in feature.qualifiers:
            alias_value = feature.qualifiers['Alias']
            if isinstance(alias_value, list) and len(alias_value) > 0:
                label = alias_value[0]
            else:
                label = str(alias_value)
            logger.debug(f"Using Alias '{label}' for feature at {feature.location}")
        
        # Priority 2: Try Name attribute
        elif 'Name' in feature.qualifiers:
            name_value = feature.qualifiers['Name']
            if isinstance(name_value, list) and len(name_value) > 0:
                label = name_value[0]
            else:
                label = str(name_value)
            logger.debug(f"Using Name '{label}' for feature at {feature.location}")
        
        # Priority 3: Try ID attribute
        elif 'ID' in feature.qualifiers:
            id_value = feature.qualifiers['ID']
            if isinstance(id_value, list) and len(id_value) > 0:
                label = id_value[0]
            else:
                label = str(id_value)
            logger.debug(f"Using ID '{label}' for feature at {feature.location}")
        
        # Priority 4: Try existing gene qualifier
        elif 'gene' in feature.qualifiers and feature.qualifiers['gene']:
            gene_value = feature.qualifiers['gene']
            if isinstance(gene_value, list) and len(gene_value) > 0:
                label = gene_value[0]
            else:
                label = str(gene_value)
        
        # Priority 5: Try locus_tag
        elif 'locus_tag' in feature.qualifiers:
            tag_value = feature.qualifiers['locus_tag']
            if isinstance(tag_value, list) and len(tag_value) > 0:
                label = tag_value[0]
            else:
                label = str(tag_value)
        
        # Priority 6: Try product (shortened for display)
        elif 'product' in feature.qualifiers:
            product_value = feature.qualifiers['product']
            if isinstance(product_value, list) and len(product_value) > 0:
                product = product_value[0]
            else:
                product = str(product_value)
            # Shorten product names for better visualization
            label = product[:20] + "..." if len(product) > 20 else product
            logger.debug(f"Using shortened product '{label}' for feature at {feature.location}")
        
        # If we found a label, set both gene and label qualifiers
        if label:
            feature.qualifiers['gene'] = [label]
            feature.qualifiers['label'] = [label]
        else:
            # No label found - use location as fallback
            location_label = f"{feature.location.start}-{feature.location.end}"
            feature.qualifiers['gene'] = [location_label]
            feature.qualifiers['label'] = [location_label]
            logger.debug(f"No gene label found, using location '{location_label}' for feature at {feature.location}")
                
    except Exception as e:
        logger.warning(f"Error extracting gene label for feature at {feature.location}: {str(e)}")
        # Set a basic label to avoid breaking visualization
        location_label = f"{feature.location.start}-{feature.location.end}"
        feature.qualifiers['gene'] = [location_label]
        feature.qualifiers['label'] = [location_label]
    
    return feature


def _fix_cds_translation(feature, sequence):
    """Ensure CDS features have correct protein translations with robust error handling.
    
    Tries multiple translation strategies:
    1. Strict CDS translation (with start/stop codons)
    2. Lenient translation (translate to first stop codon)
    3. Frame-adjusted translation (try different reading frames)
    4. Placeholder translation (if all else fails)
    
    Args:
        feature: BioPython SeqFeature object
        sequence: Parent sequence
        
    Returns:
        feature: Modified feature with translation qualifier
    """
    try:
        # Extract the CDS sequence
        cds_seq = feature.extract(sequence)
        cds_length = len(cds_seq)
        
        # Validate CDS length
        if cds_length < 3:
            logger.warning(f"CDS too short ({cds_length} bp) at {feature.location}, skipping translation")
            feature.qualifiers['translation'] = ['X']
            return feature
        
        protein_seq = None
        translation_method = None
        
        # Strategy 1: Try strict translation first (with start/stop codons)
        try:
            protein_seq = cds_seq.translate(table=1, cds=True)
            translation_method = "strict"
        except Exception as e1:
            logger.debug(f"Strict translation failed for CDS at {feature.location}: {str(e1)}")
            
            # Strategy 2: Try lenient translation (to first stop codon)
            try:
                protein_seq = cds_seq.translate(table=1, to_stop=True)
                translation_method = "to_stop"
            except Exception as e2:
                logger.debug(f"To-stop translation failed for CDS at {feature.location}: {str(e2)}")
                
                # Strategy 3: Try simple translation (no validation)
                try:
                    # Truncate to multiple of 3 if needed
                    truncate_length = (cds_length // 3) * 3
                    if truncate_length >= 3:
                        truncated_seq = cds_seq[:truncate_length]
                        protein_seq = truncated_seq.translate(table=1)
                        translation_method = "truncated"
                except Exception as e3:
                    logger.warning(f"All translation attempts failed for CDS at {feature.location}: {str(e3)}")
                    # Strategy 4: Use placeholder
                    protein_seq = Seq.Seq('X' * (cds_length // 3))
                    translation_method = "placeholder"
        
        # Update the translation qualifier
        if protein_seq:
            feature.qualifiers['translation'] = [str(protein_seq)]
            if translation_method in ["truncated", "placeholder"]:
                logger.debug(f"Used {translation_method} translation for CDS at {feature.location} (length: {cds_length} bp -> {len(protein_seq)} aa)")
            else:
                logger.debug(f"Successfully translated CDS at {feature.location} using {translation_method} method")
        else:
            # Absolute fallback
            feature.qualifiers['translation'] = ['X' * (cds_length // 3)]
            logger.warning(f"No translation possible for CDS at {feature.location}, using placeholder")
        
    except Exception as e:
        logger.error(f"Unexpected error fixing translation for CDS at {feature.location}: {str(e)}")
        # Don't fail - use placeholder to allow visualization to continue
        try:
            cds_len = len(feature.extract(sequence))
            feature.qualifiers['translation'] = ['X' * (cds_len // 3)]
        except:
            feature.qualifiers['translation'] = ['X']
        
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