import warnings

warnings.filterwarnings("ignore")

import logging

from dash import html

from io import StringIO
import base64
import re
import pandas as pd

from Bio import SeqIO

logger = logging.getLogger(__name__)


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
