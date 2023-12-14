import argparse
import pandas as pd
import os
from Bio import SearchIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Parse HMMER output and extract useful information. Options for extracting gene sequences from full starship sequence as well as placement into most likely gene family"
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


def extract_hmmer(parsed_file):
    # Read the TSV file into a DataFrame
    data = pd.read_csv(parsed_file, sep="\t")

    # Get rows with the lowest e-value for each unique entry in Query
    min_evalue_rows = data.loc[data.groupby("query_id")["evalue"].idxmin()]

    # Use os.path.join to construct file paths correctly
    top_hit_out_path = os.path.join(
        os.path.dirname(parsed_file), f"{os.path.basename(parsed_file)}.besthit.txt"
    )

    with open(top_hit_out_path, "w") as top_hit_out:
        # Write the header line
        top_hit_out.write(
            "query_id\thit_IDs\taln_length\tquery_start\tquery_end\tgaps\tquery_seq\tsubject_seq\tevalue\tbitscore\n"
        )

        # Write the rows to the file using to_csv
        min_evalue_rows.to_csv(top_hit_out, sep="\t", header=False, index=False)

        for index, row in min_evalue_rows.iterrows():
            # Create a SeqRecord
            query = row["query_id"]
            qseq = row["query_seq"]
            sequence = SeqRecord(Seq(qseq), id=query, description="")

            # Write the SeqRecord to a FASTA file
            # Use os.path.join for constructing output file paths
            output_filename = os.path.join(
                os.path.dirname(parsed_file), f"{query}_best_hsp.fa"
            )

            SeqIO.write(sequence, output_filename, "fasta")


extract_hmmer(args.parsed_file)
