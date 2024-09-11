import warnings

warnings.filterwarnings("ignore")

import dash_bootstrap_components as dbc
from dash import html

from io import StringIO
import base64
import pandas as pd

from Bio import SeqIO


def parse_fasta_from_text(text):
    try:
        queries = SeqIO.parse(StringIO(text), "fasta")
        query = next(queries)
        header, seq = str(query.id), str(query.seq)
        return header, seq
    except Exception as e:
        print(f"Error parsing text: {e}")
        return None, None


def parse_fasta_from_file(file_contents):
    try:
        if not file_contents:
            raise ValueError("No file contents provided.")
        split_contents = file_contents.split(",")
        header = split_contents[0].strip()
        sequence = "".join(split_contents[1:])
        decoded_sequence = base64.b64decode(sequence).decode("utf-8")
        query = SeqIO.read(StringIO(decoded_sequence), "fasta")
        header, seq = str(query.id), str(query.seq)
        return header, seq

    except Exception as e:
        print(f"Error parsing file: {e}")
        return None, None

        # nseq = 0
        # for query in queries:
        #     nseq += 1
        #     # Stop if there are more than one sequence records
        #     if nseq > 1:
        #         raise ValueError(
        #             "More than one sequence record found in the FASTA file."
        #         )

        #     # Extract header and sequence from the single record
        #     header, seq = query.id, str(query.seq)
        #     return header, seq

        # # Handle the case where no sequence record was found
        # raise ValueError("No sequence record found in the FASTA file.")


def parse_gff(contents, filename):
    content_type, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)
    gff = pd.read_csv(StringIO(decoded.decode("utf-8")), sep="\t")
    # cols = [
    #     "seqid",
    #     "source",
    #     "type",
    #     "start",
    #     "end",
    #     "score",
    #     "strand",
    #     "phase",
    #     "attributes",
    # ]

    # annotations = pd.DataFrame(gff)
    nanno = len(gff)

    return [
        html.Div(
            [
                html.H6(f"File name: {filename}"),
                html.H6(f"Number of annotations: {nanno}"),
            ]
        )
    ]


def parse_fasta(contents, filename):
    sequences = SeqIO.parse(StringIO(contents), "fasta")

    records = []
    nseq = 0
    for sequence in sequences:
        records.append({"ID": sequence.id, "Sequence": str(sequence.seq)})
        nseq += 1

    return [
        html.Div(
            [
                html.H6(f"File name: {filename}"),
                html.H6(f"Number of sequences: {nseq}"),
            ],
        )
    ]
