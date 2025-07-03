import re
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
import os
from src.config.logging import get_logger

logger = get_logger(__name__)


def clean_sequence(seq):
    if seq is None:
        return None

    # Remove any content within parentheses
    seq = re.sub(r"\(.*?\)", "", seq)
    seq = seq.upper()

    # Remove any non-nucleotide characters (keep IUPAC codes)
    valid_nucleotides = {
        "A",
        "T",
        "C",
        "G",
        "N",
        "R",
        "Y",
        "S",
        "W",
        "K",
        "M",
        "B",
        "D",
        "H",
        "V",
    }
    cleaned_seq = "".join(nuc for nuc in seq if nuc in valid_nucleotides)

    # Return the cleaned sequence if it's not empty
    if cleaned_seq:
        return cleaned_seq
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


def check_input(query_text_input, query_file_contents, max_sequences=10):
    """
    Checks the input type and parses the sequences.

    Args:
        query_text_input: Text input
        query_file_contents: Base64-encoded file contents
        max_sequences: Maximum number of sequences to return (default: 10)

    Returns:
        input_type: Type of input (text, file, or both [throws error])
        seq_list: List of dictionaries with sequence metadata
        n_seqs: Number of sequences
        error: Error message (if any)
    """
    error = None
    try:
        if query_text_input in ("", None) and query_file_contents is None:
            raise ValueError(
                "No input provided. Both text and file contents are empty."
            )
        elif query_text_input and query_file_contents:
            logger.warning(
                "Both text input and file contents are provided. Only one will be processed."
            )
            return "both", None, None, error
        elif query_text_input:
            input_type = "text"
            seq_list, n_seqs, error = parse_fasta_from_text(
                query_text_input, max_sequences=max_sequences
            )
            if error and isinstance(error, dmc.Alert) and error.color == "red":
                # Only treat red alerts as errors that should prevent processing
                logger.error(f"Error parsing text input: {error}")
                return None, None, None, error
        elif query_file_contents:
            input_type = "file"
            seq_list, n_seqs, warning_or_error = parse_fasta_from_file(
                query_file_contents, max_sequences=max_sequences
            )

            # Only treat real errors (red alerts) as blocking errors
            if (
                warning_or_error
                and isinstance(warning_or_error, dmc.Alert)
                and warning_or_error.color == "red"
            ):
                logger.error(f"Error parsing file contents: {warning_or_error}")
                return None, None, None, warning_or_error

        # Return early if seq_list is None
        if seq_list is None:
            logger.error("No sequences were parsed")
            return None, None, None, error if error else "No sequences could be parsed"

        # Check sequence type for each sequence in the list
        for seq_data in seq_list:
            seq = seq_data.get("sequence", "")
            seq_type = guess_seq_type(seq)
            if seq_type is None:
                logger.error(
                    f"Error guessing sequence type for {seq_data.get('header', 'unknown')}"
                )
                return None, None, None, error
            seq_data["type"] = seq_type  # Set the type in the seq_data dictionary

        # Check if all sequences are the same type
        seq_types = set(seq_data.get("type") for seq_data in seq_list)
        if len(seq_types) > 1:
            error_msg = "Sequences are of different types"
            logger.error(error_msg)
            return None, None, None, error_msg
        else:
            seq_type = list(seq_types)[0]
            logger.debug(f"All sequences are {seq_type}")

        return input_type, seq_list, n_seqs, error
    except Exception as e:
        logger.error(f"Error in check_input: {e}")
        return None, None, None, str(e)


def load_fasta_to_dict(fasta_input):
    """
    Loads sequence data into a dictionary with headers as keys and sequences as values.

    Parameters:
        fasta_input: One of:
            - A file path to a FASTA file
            - A string containing FASTA-formatted data
            - An iterable of SeqRecord objects

    Returns:
        dict: Dictionary mapping sequence IDs to sequence strings
    """
    if isinstance(fasta_input, str):
        if os.path.isfile(fasta_input):
            return {
                record.id: str(record.seq)
                for record in SeqIO.parse(fasta_input, "fasta")
            }
        else:
            fasta_io = StringIO(fasta_input)
            return {
                record.id: str(record.seq) for record in SeqIO.parse(fasta_io, "fasta")
            }
    else:
        return {record.id: str(record.seq) for record in fasta_input}


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


def seq_processing_error_alert(error):
    return dmc.Alert(
        title="Error Processing Sequence",
        children=error,
        color="red",
        variant="filled",
    )


def parse_fasta_from_text(text, format="fasta", max_sequences=10):
    """
    Parses a (multi) FASTA sequence from a text string.
    Ensures a valid FASTA header is present.
    Returns a list of dictionaries with sequence metadata, number of sequences parsed, and error message if any.
    """
    seq_list = []
    try:
        if not text or not isinstance(text, str):
            logger.error("Input text is empty or invalid type")
            return (
                None,
                None,
                dmc.Alert(
                    title="Invalid Input",
                    color="yellow",
                    children="Please enter a valid FASTA sequence.",
                ),
            )

        # Ensure proper FASTA format with header
        # only needed for a single sequence text input
        text = ensure_fasta_header(text)

        try:
            fasta_io = StringIO(text)
            sequences = list(SeqIO.parse(fasta_io, format))
        except Exception as e:
            logger.error(f"FASTA parse error: {e}")
            return (
                None,
                None,
                dmc.Alert(
                    title="FASTA Parse Error",
                    color="red",
                    children="Could not parse the input as FASTA format. Please check the sequence format.",
                ),
            )

        if len(sequences) == 0:
            logger.warning("No sequences found in input")
            return (
                None,
                None,
                dmc.Alert(
                    title="Empty Sequence",
                    color="yellow",
                    children="No valid sequence was found in the input.",
                ),
            )

        if len(sequences) > max_sequences:
            logger.warning(
                f"File contains {len(sequences)} sequences, limiting to {max_sequences}"
            )
            sequences = sequences[:max_sequences]

        for seq in sequences:
            seq_type = guess_seq_type(str(seq.seq))
            seq_list.append(
                {
                    "header": seq.id,
                    "sequence": str(seq.seq),
                    "type": seq_type,
                    "processed": False,
                }
            )

        if not seq_list:
            logger.error("Empty sequence(s) after parsing")

            return (
                None,
                None,
                dmc.Alert(
                    title="Empty Sequence",
                    color="yellow",
                    children="The sequence appears to be empty.",
                ),
            )

        return seq_list, len(sequences), None

    except ValueError as ve:
        logger.error(f"Value error in parse_fasta_from_text: {ve}")
        return (
            None,
            None,
            dmc.Alert(
                title="Invalid Format",
                color="red",
                children="The input sequence appears to be incorrectly formatted.",
            ),
        )
    except Exception as e:
        logger.error(f"Unexpected error in parse_fasta_from_text: {str(e)}")
        return (
            None,
            None,
            seq_processing_error_alert(
                "An unexpected error occurred while processing the sequence. "
                "Please check the format and try again."
            ),
        )


def parse_fasta_from_file(contents, max_sequences=10):
    """
    Parses a FASTA file that may contain multiple sequences.
    Returns a list of dictionaries with sequence metadata up to max_sequences.

    Args:
        contents: Base64-encoded file contents
        max_sequences: Maximum number of sequences to return (default: 10)

    Returns:
        List of dictionaries with sequence metadata
        Number of sequences parsed
        Warning/Error message (if any)
    """
    seq_list = []  # Initialize as empty list
    warning = None  # Use warning instead of error for non-critical issues
    n_seqs = 0

    try:
        if not contents:
            logger.warning("No file contents provided")
            return (
                None,
                None,
                dmc.Alert(
                    title="No File Contents",
                    color="yellow",
                    children="Please select a valid FASTA file.",
                ),
            )

        # Validate content format
        if "," not in contents:
            logger.error("Invalid content format")
            return (
                None,
                None,
                dmc.Alert(
                    title="Invalid File Format",
                    color="red",
                    children="The file content appears to be corrupted. Please try uploading again.",
                ),
            )

        split_contents = contents.split(",")
        if len(split_contents) < 2:
            logger.error("Invalid content split")
            return (
                None,
                None,
                dmc.Alert(
                    title="Invalid File Format",
                    color="red",
                    children="The file content is not in the expected format.",
                ),
            )

        # file_type = split_contents[0].strip()
        sequence = "".join(split_contents[1:])

        try:
            decoded_sequence = base64.b64decode(sequence).decode("utf-8")
        except Exception as e:
            logger.error(f"Base64 decode error: {e}")
            return (
                None,
                None,
                dmc.Alert(
                    title="Decoding Error",
                    color="red",
                    children="Could not decode the file contents. Please ensure you're uploading a valid FASTA file.",
                ),
            )

        fasta_io = StringIO(decoded_sequence)
        try:
            sequences = list(SeqIO.parse(fasta_io, "fasta"))
        except Exception as e:
            logger.error(f"FASTA parse error: {e}")
            return (
                None,
                None,
                dmc.Alert(
                    title="FASTA Parse Error",
                    color="red",
                    children="Could not parse the file as FASTA format. Please check the file format.",
                ),
            )

        n_seqs = len(sequences)
        if n_seqs == 0:
            logger.warning("No sequences found in file")
            return (
                None,
                None,
                dmc.Alert(
                    title="Empty FASTA File",
                    color="yellow",
                    children="No sequences were found in the uploaded file.",
                ),
            )

        # Limit number of sequences - but set a WARNING not an error
        if n_seqs > max_sequences:
            logger.warning(
                f"File contains {n_seqs} sequences, limiting to {max_sequences}"
            )
            sequences = sequences[:max_sequences]
            warning = f"Only the first {max_sequences} sequences will be processed."

        # Process sequences
        for seq in sequences:
            seq_type = guess_seq_type(str(seq.seq))
            seq_list.append(
                {
                    "header": seq.id,
                    "sequence": str(seq.seq),
                    "type": seq_type,
                    "processed": False,
                }
            )

        # If we reach this point with sequences processed successfully, return them
        # Make sure to use a yellow alert for warnings so they aren't treated as errors
        if warning:
            warning_alert = dmc.Alert(
                title="Sequence Limit Applied",
                color="yellow",
                variant="light",
                children=warning,
            )
            return seq_list, n_seqs, warning_alert

        return seq_list, n_seqs, None

    except Exception as e:
        logger.error(f"Unexpected error in parse_fasta_from_file: {str(e)}")
        return (
            None,
            None,
            dmc.Alert(
                title="Error Processing File",
                color="red",
                children=[
                    "An unexpected error occurred while processing the file. ",
                    "Please try refreshing the page and uploading again.",
                ],
            ),
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
        logger.debug(f"Parsed {nanno} annotations from {filename}")
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
                    html.H6("Error parsing GFF file."),
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

        logger.debug(f"Parsed {nseq} sequences from {filename}")
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
                    html.H6("Error parsing FASTA file."),
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


def create_ncbi_style_header(row):
    try:
        clean_contig = clean_contigIDs(row["contigID"])

        # Use accession_display if available, otherwise combine accession_tag and version_tag
        if "accession_display" in row and pd.notnull(row["accession_display"]):
            accession_with_version = row["accession_display"]
        elif (
            "version_tag" in row
            and pd.notnull(row["version_tag"])
            and row["version_tag"] != ""
        ):
            accession_with_version = f"{row['accession_tag']}.{row['version_tag']}"
        else:
            accession_with_version = row["accession_tag"]

        return (
            f">{accession_with_version} "
            f"[organism={row['name']}] "
            f"[lineage=Fungi; {row['order']}; {row['family']}] "
            f"[location={clean_contig}:{row['elementBegin']}-{row['elementEnd']}] "
            + (
                f"[assembly={row['assembly_accession']}] "
                if row["assembly_accession"]
                else ""
            )
            + f"[family={row['familyName']}]"
        )
    except Exception as e:
        logger.warning(
            f"Failed to create NCBI-style header for {row.get('accession_tag', 'unknown')}: {str(e)}"
        )
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
        records.append(
            SeqRecord(
                Seq(sequence),
                id=str(name),  # ensure ID is string
                description="",
            )
        )

    # Write to temporary file
    SeqIO.write(records, fasta_path, "fasta")
    logger.debug(f"Combined FASTA file written: {fasta_path}")


def write_combined_fasta(
    new_sequence: str,
    existing_sequences: pd.DataFrame,
    fasta_path: str,
    sequence_col: str = "sequence",
    id_col: str = None,
) -> str:
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

        records.append(
            SeqRecord(Seq(new_sequence), id="query_sequence", description="")
        )

        write_multi_fasta(existing_sequences, fasta_path, sequence_col, id_col)

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
        # Use accession_display if available, otherwise fall back to accession_tag
        if "accession_display" in existing_ships.columns:
            id_col = "accession_display"
        else:
            id_col = "accession_tag"

        sequences = {
            **fasta_sequences,
            **dict(zip(existing_ships[id_col], existing_ships["sequence"])),
        }
    else:
        sequences = fasta_sequences

    # save each as a fasta file in a temporary directory
    tmp_fasta_dir = tempfile.mkdtemp()
    for seq_id, seq in sequences.items():
        tmp_fasta = os.path.join(tmp_fasta_dir, f"{seq_id}.fa")
        write_fasta({seq_id: seq}, tmp_fasta)
    logger.debug(f"Created temporary dir for FASTA files: {tmp_fasta_dir}")
    return tmp_fasta_dir
