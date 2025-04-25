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
import dash_mantine_components as dmc
from typing import Dict

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

def ensure_fasta_header(text, default_header=">query_sequence"):
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
    Returns the header and sequence if successful, otherwise (None, None, Alert).
    """
    try:
        if not text or not isinstance(text, str):
            logger.error("Input text is empty or invalid type")
            return None, None, dmc.Alert(
                title="Invalid Input",
                color="yellow",
                children="Please enter a valid FASTA sequence.",
            )

        # Ensure proper FASTA format with header
        text = ensure_fasta_header(text)
        
        try:
            fasta_io = StringIO(text)
            sequences = list(SeqIO.parse(fasta_io, format))
        except Exception as e:
            logger.error(f"FASTA parse error: {e}")
            return None, None, dmc.Alert(
                title="FASTA Parse Error",
                color="red",
                children="Could not parse the input as FASTA format. Please check the sequence format.",
            )

        if len(sequences) == 0:
            logger.warning("No sequences found in input")
            return None, None, dmc.Alert(
                title="Empty Sequence",
                color="yellow",
                children="No valid sequence was found in the input.",
            )
        elif len(sequences) > 1:
            logger.warning("Multiple sequences found")
            return None, None, dmc.Alert(
                title="Multiple Sequences",
                color="yellow",
                children=[
                    "Multiple sequences detected in the input. ",
                    "Please enter only one sequence.",
                ],
            )

        query = sequences[0]
        header, seq = str(query.id), str(query.seq)
        
        if not seq:
            logger.error("Empty sequence after parsing")
            return None, None, dmc.Alert(
                title="Empty Sequence",
                color="yellow",
                children="The sequence appears to be empty.",
            )

        logger.info(f"Successfully parsed sequence: {header} ({len(seq)} bp)")
        return header, seq, None

    except ValueError as ve:
        logger.error(f"Value error in parse_fasta_from_text: {ve}")
        return None, None, dmc.Alert(
            title="Invalid Format",
            color="red",
            children="The input sequence appears to be incorrectly formatted.",
        )
    except Exception as e:
        logger.error(f"Unexpected error in parse_fasta_from_text: {str(e)}")
        return None, None, dmc.Alert(
            title="Error Processing Sequence",
            color="red",
            children=[
                "An unexpected error occurred while processing the sequence. ",
                "Please check the format and try again.",
            ],
        )


def parse_fasta_from_file(contents):
    header = None
    seq = None
    fasta_error = None
    try:
        if not contents:
            logger.warning("No file contents provided")
            fasta_error = dmc.Alert(
                title="No File Contents",
                color="yellow",
                children="Please select a valid FASTA file.",
            )

        # Validate content format
        if "," not in contents:
            logger.error("Invalid content format")
            fasta_error = dmc.Alert(
                title="Invalid File Format",
                color="red",
                children="The file content appears to be corrupted. Please try uploading again.",
            )

        split_contents = contents.split(",")
        if len(split_contents) < 2:
            logger.error("Invalid content split")
            fasta_error = dmc.Alert(
                title="Invalid File Format",
                color="red",
                children="The file content is not in the expected format.",
            )

        file_type = split_contents[0].strip()
        sequence = "".join(split_contents[1:])

        try:
            decoded_sequence = base64.b64decode(sequence).decode("utf-8")
        except Exception as e:
            logger.error(f"Base64 decode error: {e}")
            fasta_error = dmc.Alert(
                title="Decoding Error",
                color="red",
                children="Could not decode the file contents. Please ensure you're uploading a valid FASTA file.",
            )

        fasta_io = StringIO(decoded_sequence)
        try:
            sequences = list(SeqIO.parse(fasta_io, "fasta"))
        except Exception as e:
            logger.error(f"FASTA parse error: {e}")
            fasta_error = dmc.Alert(
                title="FASTA Parse Error",
                color="red",
                children="Could not parse the file as FASTA format. Please check the file format.",
            )

        if len(sequences) == 0:
            logger.warning("No sequences found in file")
            fasta_error = dmc.Alert(
                title="Empty FASTA File",
                color="yellow",
                children="No sequences were found in the uploaded file.",
            )
        elif len(sequences) > 1:
            logger.warning("Multiple sequences found")
            fasta_error = dmc.Alert(
                title="Multiple Sequences",
                color="yellow",
                children=[
                    "Multiple sequences detected in the file. ",
                    "Please upload a file with only one sequence.",
                ],
            )

        query = sequences[0]
        header, seq = str(query.id), str(query.seq)
        logger.info(f"Successfully parsed sequence: {header} ({len(seq)} bp)")
        return header, seq, fasta_error

    except Exception as e:
        logger.error(f"Unexpected error in parse_fasta_from_file: {str(e)}")
        return None, None, dmc.Alert(
            title="Error Processing File",
            color="red",
            children=[
                "An unexpected error occurred while processing the file. ",
                "Please try refreshing the page and uploading again.",
            ],
        )


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


def write_fasta(sequences: Dict[str, str], fasta_path: str):
    """Write a FASTA file from a dictionary of sequences.
    
    Args:
        sequences: Dictionary mapping sequence names to sequences
        fasta_path: Path to output FASTA file
    """
    with open(fasta_path, "w") as fasta_file:
        for name, sequence in sequences.items():
            fasta_file.write(f">{name}\n{sequence}\n")


def write_multi_fasta(sequences, fasta_path, sequence_col, id_col):
    records = []
    for idx, row in sequences.iterrows():
        name = row[id_col] if id_col else str(idx)
        sequence = row[sequence_col]
        records.append(SeqRecord(
            Seq(sequence),
            id=str(name),  # ensure ID is string
            description=""
        ))

    # Write to temporary file
    SeqIO.write(records, fasta_path, "fasta")
    logger.debug(f"Combined FASTA file written: {fasta_path}")

def write_combined_fasta(new_sequence: str, 
                        existing_sequences: pd.DataFrame,
                        fasta_path: str,
                        sequence_col: str = "sequence",
                        id_col: str = None) -> str:
    """Write a temporary FASTA file combining new sequence with existing sequences.
    
    Args:
        new_sequence: New sequence to classify
        existing_sequences: DataFrame containing existing sequences
        sequence_col: Name of column containing sequences
        id_col: Name of column to use as sequence IDs. If None, uses DataFrame index
    
    Returns:
        str: Path to temporary FASTA file
    """
    try:
        records = []
        
        records.append(SeqRecord(
            Seq(new_sequence),
            id="query_sequence",
            description=""
        ))
        
        write_multi_fasta(existing_sequences, fasta_path,sequence_col, id_col)
        
        
    except Exception as e:
        logger.error(f"Error writing combined FASTA: {e}")
        return None

def create_tmp_fasta_dir(fasta: str, existing_ships: pd.DataFrame) -> str:
    """Create a temporary directory for FASTA files."""
    import os
    
    # load sequence from fasta file
    fasta_sequences = load_fasta_to_dict(fasta)
    
    # append existing sequences to dict
    if existing_ships is not None and not existing_ships.empty:
        sequences = {
            **fasta_sequences,
            **dict(zip(existing_ships['accession_tag'], existing_ships['sequence']))
        }
    else:
        sequences = fasta_sequences

    # save each as a fasta file in a temporary directory
    tmp_fasta_dir = tempfile.mkdtemp()
    for seq_id, seq in sequences.items():
        tmp_fasta = os.path.join(tmp_fasta_dir, f"{seq_id}.fa")
        write_fasta(
            {seq_id: seq},
            tmp_fasta
        )
    logger.debug(f"Created temporary dir for FASTA files: {tmp_fasta_dir}")
    return tmp_fasta_dir
