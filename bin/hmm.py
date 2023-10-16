import argparse
from Bio import SearchIO
from Bio import AlignIO

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
    "-o",
    "--parsed",
    dest="parsed_file",
    help="Output file",
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
def parse_hmmer(hmmer_output_file, parsed_file):
    with open(parsed_file, "w") as tsv_file:
        tsv_file.write("Query\tSubject\tQuery Sequence\tSubject Sequence\tevalue\n")
        for record in SearchIO.parse(hmmer_output_file, "hmmer3-text"):
            for hit in record.hits:
                for hsp in hit.hsps:
                    query_seq = str(hsp.query.seq)
                    subject_seq = str(hsp.hit.seq)
                    evalue = hsp.evalue
                    tsv_file.write(
                        f"{record.id}\t{hit.id}\t{query_seq}\t{subject_seq}\t{evalue}\n"
                    )


parse_hmmer(args.hmmer_output_file, args.parsed_file)

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
                # Write the SeqRecord objects to a FASTA file
                # TODO: where to save these files
                output_filename = (
                    os.path.dirname(hmmer_output_file) + query_name + ".fa"
                )
                with open(output_filename, "w") as output_handle:
                    SeqIO.write(sequences[query_name], output_handle, "fasta")


if args.ship_seq:
    extract_hmmer(args.hmmer_output_file, args.ship_seq)
