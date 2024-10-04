import warnings

warnings.filterwarnings("ignore")

import logging

logging.basicConfig(level=logging.DEBUG)

from dash import html

from io import StringIO
import base64
import pandas as pd
import tempfile
from Bio import SeqIO


def ensure_fasta_header(text, default_header=">query"):
    """
    Ensures that the input contains a valid FASTA header.
    If no header is found, adds a default one.
    """
    lines = text.strip().splitlines()

    if not lines or not lines[0].startswith(">"):
        logging.warning("No FASTA header found. Adding default header.")
        text = default_header + "\n" + text

    return text


def extract_fasta_headers(fasta_file, out_type="list"):
    """Extract headers from the provided FASTA file and saves them in a temporary file."""

    headers = []
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            headers.append(record.id)
        if out_type == "file":
            tmp_headers = tempfile.NamedTemporaryFile(suffix=".headers", delete=False)
            with open(tmp_headers.name, "w") as f:
                for header in headers:
                    f.write(header + "\n")
            # logging.debug(f"Headers extracted and saved to: {tmp_headers.name}")
            return tmp_headers.name
        if out_type == "list":
            return headers
    except Exception as e:
        logging.error(f"Error extracting headers: {e}")
        return None


def parse_fasta_from_text(text, format="fasta"):
    """
    Parses a FASTA sequence from a text string.
    Ensures a valid FASTA header is present.
    Returns the header and sequence if successful, otherwise (None, None).
    """
    try:
        if not text:
            logging.error("Input text is empty.")
            raise
        text = ensure_fasta_header(text)
        queries = SeqIO.parse(StringIO(text), format)
        query = next(queries, None)  # Use `None` to avoid StopIteration

        if query is None:
            logging.error("No sequence found in the provided text.")
            raise

        header, seq = str(query.id), str(query.seq)
        logging.info(f"Parsed sequence: {header}")
        return header, seq

    except ValueError as ve:
        logging.error(f"Value error: {ve}")
        return None, None
    except Exception as e:
        logging.error(f"Error parsing text: {e}")
        return None, None


def parse_fasta_from_file(file_contents):
    """
    Parses FASTA sequences from a file's contents (base64 encoded).
    Ensures valid FASTA headers are present.
    Returns a list of tuples containing headers and sequences if successful,
    otherwise an empty list.
    """
    try:
        if not file_contents:
            raise ValueError("No file contents provided.")

        split_contents = file_contents.split(",")
        file_type = split_contents[0].strip()
        sequence = "".join(split_contents[1:])

        decoded_sequence = ensure_fasta_header(
            base64.b64decode(sequence).decode("utf-8")
        )

        fasta_sequences = SeqIO.parse(StringIO(decoded_sequence), "fasta")

        parsed_sequences = []
        for query in fasta_sequences:
            header, seq = str(query.id), str(query.seq)
            parsed_sequences.append((header, seq))
            logging.info(f"Parsed sequence: {header}")

        return parsed_sequences
    except ValueError as ve:
        logging.error(f"Value error: {ve}")
        return []
    except Exception as e:
        logging.error(f"Error parsing file: {e}")
        return []


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
        logging.info(f"Parsed {nanno} annotations from {filename}")
        return [
            html.Div(
                [
                    html.H6(f"File name: {filename}"),
                    html.H6(f"Number of annotations: {nanno}"),
                ]
            )
        ]

    except Exception as e:
        logging.error(f"Error parsing GFF file {filename}: {e}")
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

        logging.info(f"Parsed {nseq} sequences from {filename}")
        return [
            html.Div(
                [
                    html.H6(f"File name: {filename}"),
                    html.H6(f"Number of sequences: {nseq}"),
                ]
            )
        ]

    except Exception as e:
        logging.error(f"Error parsing FASTA file {filename}: {e}")
        return [
            html.Div(
                [
                    html.H6(f"File name: {filename}"),
                    html.H6(f"Error parsing FASTA file."),
                ]
            )
        ]
