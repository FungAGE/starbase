import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Parse HMMER output and extract gene sequences from full starship sequence"
)
parser.add_argument(
    "-i",
    "--hmmer_output_file",
    dest="hmmer_output_file",
    help="Path to output file from HMMER search, in `domtblout` format",
    required=True,
)
parser.add_argument(
    "-q",
    "--ship_seq",
    dest="ship_seq",
    help="Path to query sequence FASTA to extract gene sequences from",
    required=False,
)
args = parser.parse_args()


# Parse the HMMER results
def parse_hmmer(hmmer_output_file):
    with open(hmmer_output_file, "r") as hmmer_output:
        for line in hmmer_output:
            if line.startswith("#"):
                continue
            fields = line.strip().split()
            hit_name = fields[0]
            length = int(fields[2])
            query_name = fields[3]
            eval = fields[11]
            score = fields[13]

            df = pd.DataFrame(
                {
                    "query_ID": query_name,
                    "hit_IDs": hit_name,
                    "hit_length": length,
                    "bitscore": score,
                    "eval": eval,
                },
                index=[0],
            )
            df.to_csv("tmp/hmmer-parsed.csv")


parse_hmmer(args.hmmer_output_file)

# Create a dictionary to store sequences
sequences = {}


def extract_hmmer(hmmer_output_file, ship_seq):
    with open(hmmer_output_file, "r") as hmmer_output:
        for line in hmmer_output:
            if line.startswith("#"):
                continue
            fields = line.strip().split()
            hit_name = fields[0]
            query_name = fields[3]
            start = int(fields[19])
            end = int(fields[20])
            # Load the query sequence
            query_sequences = SeqIO.to_dict(SeqIO.parse(ship_seq, "fasta"))
            if query_name in query_sequences:
                hit_sequence = query_sequences[query_name]
                gene_sequence = hit_sequence[
                    start - 1 : end
                ]  # Adjust for 0-based indexing
                sequences[query_name] = gene_sequence
                # TODO: SeqIO write fasta
                # print(sequences[query_name])


if args.ship_seq:
    extract_hmmer(args.hmmer_output_file, args.ship_seq)
