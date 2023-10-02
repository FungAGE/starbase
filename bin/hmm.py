import argparse
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO

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
parser.add_argument(
    "-a",
    "--alignment",
    dest="hmmer_aln_file",
    help="Path to alignment from output of HMMER search",
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


def hmmer_aln(hmmer_aln_file):
    # Parse the Stockholm file and open the CSV file for writing
    with open(hmmer_aln_file, "r") as stockholm_handle, open(
        csv_file, "w", newline=""
    ) as csv_handle:
        # Create a Stockholm parser
        stockholm_alignments = AlignIO.parse(stockholm_handle, "stockholm")

        # Initialize a CSV writer
        csv_writer = csv.writer(csv_handle)

        # Write the header row to the CSV file
        csv_writer.writerow(["Sequence Name", "Aligned Sequence"])

        # Iterate through each alignment in the Stockholm file
        for alignment in stockholm_alignments:
            for record in alignment:
                # Write the sequence name and aligned sequence to the CSV file
                csv_writer.writerow([record.id, str(record.seq)])

    # Load the Stockholm format alignment file
    alignment = AlignIO.parse("alignment.sto", "stockholm")

    # Create an empty DataFrame to store the alignment data
    alignment_data = []

    # Iterate through the alignment and extract sequence data
    for record in alignment:
        sequence_data = {
            "ID": record.id,
            "Description": record.description,
            "Sequence": str(record.seq),
        }
        alignment_data.append(sequence_data)

    # Convert the alignment data to a DataFrame
    alignment_df = pd.DataFrame(alignment_data)

    # Save the alignment as a CSV file
    alignment_df.to_csv("tmp/alignment.csv", index=False)


hmmer_aln(args.hmmer_aln_file)

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
