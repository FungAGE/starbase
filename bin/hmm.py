import argparse
import pandas as pd
from Bio import SearchIO
from Bio import AlignIO
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Parse HMMER output and extract gene sequences from full starship sequence"
)
parser.add_argument(
    "-i",
    "--hmmer_output_file",
    dest="hmmer_output_file",
    help="Path to the 'hmmer3-text' format file",
    required=True,
)
parser.add_argument(
    "-o", "--parsed", dest="parsed_file", help="Output file", required=True
)
args = parser.parse_args()


# Parse the HMMER results
def parse_hmmer(hmmer_output_file, parsed_file):
    with open(parsed_file, "w") as tsv_file:
        tsv_file.write(
            "query_id\thit_IDs\taln_length\tquery_start\tquery_end\tgaps\tquery_seq\tsubject_seq\tevalue\tbitscore\n"
        )
        for record in SearchIO.parse(hmmer_output_file, "hmmer3-text"):
            for hit in record.hits:
                for hsp in hit.hsps:
                    query_seq = str(hsp.query.seq)
                    subject_seq = str(hsp.hit.seq)
                    aln_length = hsp.aln_span
                    query_start = hsp.query_start
                    query_end = hsp.query_end
                    gaps = str("N/A")
                    bitscore = hsp.bitscore
                    evalue = hsp.evalue
                    tsv_file.write(
                        f"{record.id}\t{hit.id}\t{aln_length}\t{query_start}\t{query_end}\t{gaps}\t{query_seq}\t{subject_seq}\t{evalue}\t{bitscore}\n"
                    )


parse_hmmer(args.hmmer_output_file, args.parsed_file)

# Create a dictionary to store sequences
sequences = {}


# TODO: just return top hit
def extract_hmmer(parsed_file):
    # Read the TSV file into a DataFrame
    data = pd.read_csv(parsed_file, sep="\t")

    # Get rows with the lowest e-value for each unique entry in Query
    min_evalue_rows = data.loc[data.groupby("query_id")["evalue"].idxmin()]

    # Create individual FASTA files
    for index, row in min_evalue_rows.iterrows():
        query = row["query_id"]
        query_sequence = row["query_seq"]

        # Create a SeqRecord
        sequence = SeqIO.SeqRecord(SeqIO.Seq(query_sequence), id=query, description="")

        # Write the SeqRecord to a FASTA file
        # TODO: create directory for output
        output_filename = os.path.dirname(parsed_file) + f"/{query}_best_hsp.fa"

        SeqIO.write(sequence, output_filename, "fasta")


extract_hmmer(args.parsed_file)
