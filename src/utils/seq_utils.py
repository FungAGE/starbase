import re
import warnings
import logging
import tempfile
from io import StringIO
import base64
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from Bio import SeqIO
from dash import html

warnings.filterwarnings("ignore")

logger = logging.getLogger(__name__)

def clean_sequence(seq):
    valid_nucleotides = {"A", "T", "C", "G"}
    seq = re.sub(r"\(.*?\)", "", seq)
    seq = seq.upper()
    if all(nuc in valid_nucleotides for nuc in seq):
        return seq
    else:
        return None

def clean_lines(queries):
    try:
        cleaned_seq = ""
        for query in queries:
            cleaned_seq += re.sub("[^A-Za-z]", "", str(query))

        return cleaned_seq
    except Exception as e:
        logger.error(f"Error cleaning lines: {e}")
        return None

def guess_seq_type(query_seq):
    nucl_char = set("ATGC")
    prot_char = set("ARNDBCEQZGHILKMFPSTWYV")

    try:
        if query_seq is None:
            logger.error("query_seq is None.")
            return None

        cleaned_seq = clean_lines(query_seq)

        nucl_count = sum(1 for nt in query_seq.upper() if nt in nucl_char)
        prot_count = sum(1 for aa in query_seq.upper() if aa in prot_char)

        query_type = "nucl" if nucl_count >= (0.8 * len(cleaned_seq)) else "prot"
        logger.debug(
            f"Cleaned sequence length: {len(cleaned_seq)}, Sequence type guessed: {query_type}, Nucleotide count: {nucl_count}, Protein count: {prot_count}"
        )

        return query_type
    except Exception as e:
        logger.error(f"Error in guessing sequence type: {e}")
        return None

def write_temp_fasta(header, sequence):
    try:
        cleaned_query_seq = SeqRecord(Seq(sequence), id=header, description="")
        tmp_query_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
        SeqIO.write(cleaned_query_seq, tmp_query_fasta, "fasta")
        logger.debug(f"Temporary FASTA file written: {tmp_query_fasta}")
        return tmp_query_fasta
    except Exception as e:
        logger.error(f"Error writing temporary FASTA: {e}")
        return None


def check_input(query_text_input, query_file_contents):
    try:
        if query_text_input in ("", None) and query_file_contents is None:
            raise ValueError(
                "No input provided. Both text and file contents are empty."
            )
        elif query_text_input and query_file_contents:
            logger.warning(
                "Both text input and file contents are provided. Only one will be processed."
            )
            return "both", None, None
        elif query_text_input:
            input_type = "text"
            header, query = parse_fasta_from_text(query_text_input)
        elif query_file_contents:
            input_type = "file"
            header, query, error_message = parse_fasta_from_file(query_file_contents)
        logger.debug(
            f"Input type: {input_type}, Header: {header}, Query Length: {len(query) if query else 'None'}"
        )
        return input_type, header, query
    except Exception as e:
        logger.error(f"Error in check_input: {e}")
        return None, None, None

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
        fwd = s[phase:].translate(to_stop=False)
        rev = s.reverse_complement()[phase:].translate(to_stop=False)
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


def get_protein_sequence(header, nuc_sequence):
    """
    Translates the nucleotide sequence to protein sequence and writes to a multi-sequence FASTA file.

    Parameters:
        headers (list): List of headers for each protein sequence.
        nuc_sequence (str): The nucleotide sequence to be translated.

    Returns:
        str: The filename of the generated FASTA file.
    """
    protein_output_filename = tempfile.NamedTemporaryFile(
        suffix=".faa", delete=False
    ).name
    protein_seqs = translate_seq(nuc_sequence)

    if protein_seqs is not None:
        if isinstance(protein_seqs, list):
            records = []
            for i, protein_seq in enumerate(protein_seqs):
                new_header = f"{header}_{i}"
                records.append(
                    SeqRecord(Seq(protein_seq), id=new_header, description="")
                )

            SeqIO.write(records, protein_output_filename, "fasta")
        return protein_output_filename
    else:
        return None

def ensure_fasta_header(text, default_header=">query"):
    """
    Ensures that the input contains a valid FASTA header.
    If no header is found, adds a default one.
    """
    lines = text.strip().splitlines()

    if not lines or not lines[0].startswith(">"):
        logger.warning("No FASTA header found. Adding default header.")
        text = default_header + "\n" + text

    return text


def parse_fasta_from_text(text, format="fasta"):
    """
    Parses a FASTA sequence from a text string.
    Ensures a valid FASTA header is present.
    Returns the header and sequence if successful, otherwise (None, None).
    """
    try:
        if not text:
            logger.error("Input text is empty.")
            raise
        text = ensure_fasta_header(text)
        queries = SeqIO.parse(StringIO(text), format)
        query = next(queries, None)  # Use `None` to avoid StopIteration

        if query is None:
            logger.error("No sequence found in the provided text.")
            raise

        header, seq = str(query.id), str(query.seq)
        logger.info(f"Parsed sequence: {header}")
        return header, seq

    except ValueError as ve:
        logger.error(f"Value error: {ve}")
        return None, None
    except Exception as e:
        logger.error(f"Error parsing text: {e}")
        return None, None


def parse_fasta_from_file(file_contents):
    """
    Parses a FASTA sequence from a file's contents (base64 encoded).
    Ensures a valid FASTA header is present and only one sequence is present.
    Returns the header and sequence if successful, otherwise (None, None).
    """
    try:
        if not file_contents:
            raise ValueError("No file contents provided.")

        split_contents = file_contents.split(",")
        file_type = split_contents[0].strip()
        sequence = "".join(split_contents[1:])

        decoded_sequence = base64.b64decode(sequence).decode("utf-8")

        fasta_io = StringIO(decoded_sequence)
        sequences = list(SeqIO.parse(fasta_io, "fasta"))

        if len(sequences) == 0:
            raise ValueError("No valid FASTA sequence found.")
        elif len(sequences) > 1:
            raise ValueError("More than one sequence found in the FASTA file.")

        query = sequences[0]
        header, seq = str(query.id), str(query.seq)
        logger.info(f"Parsed sequence: {header}")
        return header, seq, None

    except ValueError as ve:
        logger.error(f"Value error: {ve}")
        return header, seq, str(ve)
    except Exception as e:
        logger.error(f"Error parsing file: {e}")
        return header, seq, str(e)


def parse_gff(contents, filename):
    """
    Parses a GFF file from base64-encoded contents and returns a brief summary.
    """
    try:
        content_type, content_string = contents.split(",")
        decoded = base64.b64decode(content_string)

        # Assuming GFF is tab-separated
        gff = pd.read_csv(StringIO(decoded.decode("utf-8")), sep="\t")

        nanno = len(gff)
        logger.info(f"Parsed {nanno} annotations from {filename}")
        return [
            html.Div(
                [
                    html.H6(f"File name: {filename}"),
                    html.H6(f"Number of annotations: {nanno}"),
                ]
            )
        ]

    except Exception as e:
        logger.error(f"Error parsing GFF file {filename}: {e}")
        return [
            html.Div(
                [
                    html.H6(f"File name: {filename}"),
                    html.H6(f"Error parsing GFF file."),
                ]
            )
        ]


def parse_fasta(contents, filename):
    """
    Parses a multi-sequence FASTA file and returns a brief summary of sequences.
    """
    try:
        sequences = SeqIO.parse(StringIO(contents), "fasta")
        records = []
        nseq = 0

        for sequence in sequences:
            records.append({"ID": sequence.id, "Sequence": str(sequence.seq)})
            nseq += 1

        logger.info(f"Parsed {nseq} sequences from {filename}")
        return [
            html.Div(
                [
                    html.H6(f"File name: {filename}"),
                    html.H6(f"Number of sequences: {nseq}"),
                ]
            )
        ]

    except Exception as e:
        logger.error(f"Error parsing FASTA file {filename}: {e}")
        return [
            html.Div(
                [
                    html.H6(f"File name: {filename}"),
                    html.H6(f"Error parsing FASTA file."),
                ]
            )
        ]


def clean_shipID(s):
    return re.sub(
        r"_(\d)",
        "",
        re.sub(
            "-Captain_.*",
            "",
            s.replace("|-", "").replace("|+", "").replace("|", "").replace("fam_", "_"),
        ),
    )


def clean_contigIDs(string):
    """Removing omes from contigIDs."""
    if string is None:
        return None
    else:
        parts = string.split("_", 1)
        if len(parts) > 1:
            prefix = parts[0]
            suffix = parts[1]
            if 7 <= len(prefix) <= 9:
                return suffix
            else:
                return string
