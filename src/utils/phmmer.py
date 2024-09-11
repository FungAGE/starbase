import os
import pandas as pd
from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from Bio.SeqUtils import nt_search, six_frame_translations
import subprocess
import argparse


def load_fasta_to_dict(fasta_file):
    """
    Loads a multi-FASTA file into a dictionary with headers as keys and sequences as values.
    """
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}


def find_start_codons(seq):
    """
    Finds the start codon positions (ATG) in the nucleotide sequence.
    Returns the first start codon position or None if not found.
    """
    start_codon = "ATG"
    start_positions = nt_search(str(seq), start_codon)[1:]  # Skip the sequence itself

    if start_positions:
        return start_positions[0]
    return None


def translate_seq(seq):
    """
    Translate the nucleotide sequence in three forward and three reverse frames.
    """
    s = Seq(seq)
    frames = []
    for phase in range(3):
        fwd = s[phase:].translate(to_stop=True)
        rev = s.reverse_complement()[phase:].translate(to_stop=True)
        frames.append(fwd)
        frames.append(rev)
    return frames


def find_longest_orf(aa_seq, min_length=50):
    """
    Finds the longest ORF in the translated amino acid sequence, ignoring stretches of 'X'.
    """
    longest_orf = ""
    current_orf = []

    for aa in aa_seq:
        if aa == "*":
            orf = "".join(current_orf)
            if (
                len(orf) >= min_length
                and "X" not in orf
                and len(orf) > len(longest_orf)
            ):
                longest_orf = orf
            current_orf = []
        else:
            current_orf.append(aa)

    # In case there's an ORF at the end of the sequence
    orf = "".join(current_orf)
    if len(orf) >= min_length and "X" not in orf and len(orf) > len(longest_orf):
        longest_orf = orf

    return longest_orf


def get_protein_sequence(header, nuc_sequence, protein_output_filename):
    """
    Translates the nucleotide sequence to protein sequence and writes to a FASTA file.
    """
    # start_codon = find_start_codons(nuc_sequence)
    # if start_codon is not None:
    #     protein_seqs = translate(nuc_sequence[start_codon:])
    # else:
    protein_seqs = translate_seq(nuc_sequence)

    longest_protein = max((find_longest_orf(frame) for frame in protein_seqs), key=len)

    # Write the longest ORF to the output file
    SeqIO.write(
        SeqRecord(Seq(longest_protein), id=header, description=""),
        protein_output_filename,
        "fasta",
    )

    return len(Seq(longest_protein))


def run_hmmersearch(query_fasta, hmm_db, output_file, threads=2, eval_threshold=1):
    """
    Runs the HMMER search using the query FASTA file and HMM database.
    """
    hmmer_cmd = [
        "phmmer",
        "-o",
        output_file,
        "--cpu",
        str(threads),
        "--domE",
        str(eval_threshold),
        hmm_db,
        query_fasta,
    ]
    subprocess.run(hmmer_cmd, check=True)


def parse_hmmer(hmmer_output_file):
    """
    Parses HMMER output and extracts relevant information.
    """
    parsed_hits = []

    for record in SearchIO.parse(hmmer_output_file, "hmmer3-text"):
        for hit in record.hits:
            for hsp in hit.hsps:
                parsed_hits.append(
                    {
                        "query_id": record.id,
                        "hit_IDs": hit.id,
                        "aln_length": hsp.aln_span,
                        "query_start": hsp.query_start,
                        "query_end": hsp.query_end,
                        "gaps": "N/A",  # You can change this as needed
                        "query_seq": str(hsp.query.seq),
                        "subject_seq": str(hsp.hit.seq),
                        "evalue": hsp.evalue,
                        "bitscore": hsp.bitscore,
                    }
                )

    return parsed_hits


def main():
    parser = argparse.ArgumentParser(
        description="Given a nucleotide fasta file, use all ORFs for searching HMMs."
    )
    parser.add_argument(
        "-i", "--input_fasta", required=True, help="Path to the input fasta file."
    )
    parser.add_argument("-d", "--hmmer_db", required=True, help="Path to HMM database.")
    parser.add_argument(
        "-o", "--output_dir", required=True, help="Path to output directory."
    )

    args = parser.parse_args()

    threads = 2
    input_eval = 0.01

    # Load the sequences
    seq_dict = load_fasta_to_dict(args.input_fasta)

    all_results = []
    for header, seq in seq_dict.items():
        print(f"Processing {header}...")
        protein_output = os.path.join(args.output_dir, f"{header}.fa")
        hmmer_output = os.path.join(args.output_dir, f"{header}.txt")

        # Get the protein sequence
        len = get_protein_sequence(header, seq, protein_output)

        if len > 1:

            # Run HMMER search
            run_hmmersearch(
                protein_output, args.hmmer_db, hmmer_output, threads, input_eval
            )

            # Parse the HMMER output
            parsed_hits = parse_hmmer(hmmer_output)
            all_results.extend(parsed_hits)

    # Create a DataFrame from all parsed HMMER results
    df = pd.DataFrame(all_results)
    parsed_file = os.path.join(
        args.output_dir, f"{os.path.basename(args.input_fasta)}.tsv"
    )
    df.to_csv(parsed_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
