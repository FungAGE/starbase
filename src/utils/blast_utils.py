import warnings

warnings.filterwarnings("ignore")

import logging

logging.basicConfig(level=logging.DEBUG)

import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html
import dash_bio as dashbio

import os
import re
import tempfile
import subprocess
import json
import pandas as pd

from Bio.Blast.Applications import (
    NcbiblastnCommandline,
    NcbitblastnCommandline,
    NcbiblastpCommandline,
    NcbiblastxCommandline,
)
from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


from src.utils.parsing import parse_fasta_from_file, parse_fasta_from_text


def write_temp_fasta(header, sequence):
    try:
        cleaned_query_seq = SeqRecord(Seq(sequence), id=header, description="")
        tmp_query_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
        SeqIO.write(cleaned_query_seq, tmp_query_fasta, "fasta")
        logging.debug(f"Temporary FASTA file written: {tmp_query_fasta}")
        return tmp_query_fasta
    except Exception as e:
        logging.error(f"Error writing temporary FASTA: {e}")
        return None


def check_input(query_text_input, query_file_contents):
    try:
        if query_text_input in ("", None) and query_file_contents is None:
            raise ValueError(
                "No input provided. Both text and file contents are empty."
            )
        elif query_text_input and query_file_contents:
            logging.warning(
                "Both text input and file contents are provided. Only one will be processed."
            )
            return "both", None, None
        elif query_text_input:
            input_type = "text"
            header, query = parse_fasta_from_text(query_text_input)
        elif query_file_contents:
            input_type = "file"
            header, query = parse_fasta_from_file(query_file_contents)
        logging.debug(
            f"Input type: {input_type}, Header: {header}, Query Length: {len(query) if query else 'None'}"
        )
        return input_type, header, query
    except Exception as e:
        logging.error(f"Error in check_input: {e}")
        return None, None, None


nucl_char = set("ATGC")
prot_char = set("ARNDBCEQZGHILKMFPSTWYV")


def clean_lines(queries):
    try:
        cleaned_seq = ""
        for query in queries:
            cleaned_seq += re.sub("[^A-Za-z]", "", str(query))

        return cleaned_seq
    except Exception as e:
        logging.error(f"Error cleaning lines: {e}")
        return None


def guess_seq_type(query_seq):
    try:
        if query_seq is None:
            logging.error("query_seq is None.")
            return None

        cleaned_seq = clean_lines(query_seq)

        nucl_count = sum(1 for nt in query_seq.upper() if nt in nucl_char)
        prot_count = sum(1 for aa in query_seq.upper() if aa in prot_char)

        query_type = "nucl" if nucl_count >= (0.8 * len(cleaned_seq)) else "prot"
        logging.debug(
            f"Cleaned sequence length: {len(cleaned_seq)}, Sequence type guessed: {query_type}, Nucleotide count: {nucl_count}, Protein count: {prot_count}"
        )

        return query_type
    except Exception as e:
        logging.error(f"Error in guessing sequence type: {e}")
        return None


def run_blast(
    db_list=None,
    query_type=None,
    query_fasta=None,
    tmp_blast=None,
    input_eval=None,
    threads=None,
):
    try:
        db_type = "nucl"

        blastdb = db_list["ship"][db_type]
        if db_type == "nucl":
            if query_type == "nucl":
                blast_program = NcbiblastnCommandline
            else:
                blast_program = NcbitblastnCommandline
        else:
            if query_type == "prot":
                blast_program = NcbiblastpCommandline
            else:
                blast_program = NcbiblastxCommandline

        if not os.path.exists(blastdb) or os.path.getsize(blastdb) == 0:
            raise ValueError(f"BLAST database {blastdb} not found or is empty.")

        logging.info(f"Running BLAST with query: {query_fasta}, Database: {blastdb}")
        blast_cline = blast_program(
            query=query_fasta,
            db=blastdb,
            evalue=input_eval,
            out=tmp_blast,
            outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq",
            num_threads=threads,
        )
        stdout, stderr = blast_cline()
        logging.debug(f"BLAST stdout: {stdout}, stderr: {stderr}")

        # Optionally stitch BLAST results
        # if stitch:
        #     stitched_blast_tmp = tempfile.NamedTemporaryFile(suffix=".stitch").name
        #     stitch_blast_cmd = f"python bin/BLASTstitcher.py -i {tmp_blast} -o {stitched_blast_tmp}"
        #     subprocess.run(stitch_blast_cmd, shell=True)
        #     ship_blast_out = stitched_blast_tmp

        # with open(tmp_query_fasta, "r") as file:
        #     file_contents = file.read()
        # logging.info(file_contents)

        # with open(tmp_blast, "r") as file:
        #     file_contents = file.read()
        # logging.info(file_contents)

        df = pd.read_csv(
            tmp_blast,
            sep="\t",
            names=[
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
                "qseq",
                "sseq",
            ],
        )
        
        df["qseqid"] = df["qseqid"].replace("|-", "").replace("|+", "").replace("|", "")
        
        logging.info(f"BLAST results parsed with {len(df)} hits.")

        return df
    except Exception as e:
        logging.error(f"Error during BLAST search: {e}")
        return None


def run_hmmer(
    db_list=None,
    query_type=None,
    input_genes="tyr",
    input_eval=None,
    query_fasta=None,
    tmp_hmmer=None,
    tmp_hmmer_parsed=None,
    threads=None,
):
    try:
        hmmer_db = db_list["gene"][input_genes]["hmm"][query_type]
        if not os.path.exists(hmmer_db) or os.path.getsize(hmmer_db) == 0:
            raise ValueError(f"HMMER database {hmmer_db} not found or is empty.")

        hmmer_cmd = f"hmmsearch -o {tmp_hmmer} --cpu {threads} --domE {input_eval} {hmmer_db} {query_fasta}"
        logging.info(f"Running HMMER search: {hmmer_cmd}")
        subprocess.run(hmmer_cmd, shell=True)

        parse_hmmer(tmp_hmmer, tmp_hmmer_parsed)
        subject_seq = extract_gene_from_hmmer(tmp_hmmer_parsed)

        hmmer_results = pd.read_csv(tmp_hmmer_parsed, sep="\t")
        logging.debug(f"HMMER results parsed: {hmmer_results.shape[0]} rows.")

        return hmmer_results, subject_seq
    except Exception as e:
        logging.error(f"Error in HMMER search: {e}")
        return None, None


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
                    subject_seq = re.sub("-Captain_.*", "", str(hsp.hit.seq).replace("|-", "").replace("|+", "").replace("|", ""))
                    aln_length = hsp.aln_span
                    query_start = hsp.query_start
                    query_end = hsp.query_end
                    gaps = str("N/A")
                    bitscore = hsp.bitscore
                    evalue = hsp.evalue
                    tsv_file.write(
                        f"{hit.id}\t{record.id}\t{aln_length}\t{query_start}\t{query_end}\t{gaps}\t{query_seq}\t{subject_seq}\t{evalue}\t{bitscore}\n"
                    )


def extract_gene_from_hmmer(parsed_file):
    # Read the TSV file into a DataFrame
    data = pd.read_csv(parsed_file, sep="\t")

    # Get rows with the lowest e-value for each unique entry in Query
    min_evalue_rows = data.loc[data.groupby("query_id")["evalue"].idxmin()]

    # Use os.path.join to construct file paths correctly
    top_hit_out_path = os.path.join(
        os.path.dirname(parsed_file), f"{os.path.splitext(parsed_file)[0]}.besthit.txt"
    )

    output_filename = None
    with open(top_hit_out_path, "w") as top_hit_out:
        # Write the header line
        top_hit_out.write(
            "hit_IDs\tquery_id\taln_length\tquery_start\tquery_end\tgaps\tquery_seq\tsubject_seq\tevalue\tbitscore\n"
        )

        # Write the rows to the file using to_csv
        min_evalue_rows.to_csv(top_hit_out, sep="\t", header=False, index=False)

        for index, row in min_evalue_rows.iterrows():
            # Create a SeqRecord
            query = row["query_id"]
            qseq = re.sub(r"\.", "", str(row["query_seq"]))
            sequence = SeqRecord(Seq(qseq), id=query, description="")

            # Write the SeqRecord to a FASTA file
            # Use os.path.join for constructing output file paths
            output_filename = os.path.join(
                os.path.dirname(parsed_file), f"{query}_best_hsp.fa"
            )

            SeqIO.write(sequence, output_filename, "fasta")

    output = html.Div(
        [
            dbc.Button(
                "Download best captain hit",
                id="subject-seq-button",
                n_clicks=0,
                className="d-grid gap-2 col-6 mx-auto",
                style={"fontSize": "1rem"},
            ),
            dcc.Download(id="subject-seq-dl-package"),
        ]
    )
    return output


def circos_prep(blast_output, links_output, layout_output):
    best_blast = blast_output.nsmallest(10, "evalue")
    links = [
        {
            "source": {
                "id": row["qseqid"],
                "start": int(row["qstart"]),
                "end": int(row["qend"]),
            },
            "target": {
                "id": row["sseqid"],
                "start": int(row["sstart"]),
                "end": int(row["send"]),
            },
            # Use pident as the single numeric value for simplicity
            "value": row["pident"],
        }
        for index, row in best_blast.iterrows()
    ]
    # Save to JSON file
    with open(links_output, "w") as json_file:
        json.dump(links, json_file, indent=4)

    # Determine unique sequence IDs and their lengths
    sequence_lengths = {}
    for index, row in best_blast.iterrows():
        if row["qseqid"] not in sequence_lengths:
            sequence_lengths[row["qseqid"]] = int(row["qend"])
        else:
            sequence_lengths[row["qseqid"]] = max(
                sequence_lengths[row["qseqid"]], int(row["qend"])
            )

        if row["sseqid"] not in sequence_lengths:
            sequence_lengths[row["sseqid"]] = int(row["send"])
        else:
            sequence_lengths[row["sseqid"]] = max(
                sequence_lengths[row["sseqid"]], int(row["send"])
            )

    # Define colors for each chromosome block
    colors = [
        "#FF5733",
        "#33FF57",
        "#3357FF",
        "#FF33A6",
        "#33FFF1",
        "#F1FF33",
        "#B633FF",
        "#33FFB6",
        "#FFA633",
        "#33A6FF",
    ]

    # Create Circos layout with colors
    layout = [
        {
            "id": seq_id,
            "len": length,
            "color": colors[i % len(colors)],
            "label": {"display": seq_id, "color": "#000000", "size": 14},
        }
        for i, (seq_id, length) in enumerate(sequence_lengths.items())
    ]

    # Save layout to JSON file
    with open(layout_output, "w") as json_file:
        json.dump(layout, json_file, indent=4)


def blast_chords(blast_output):
    if blast_output is not None and not blast_output.empty:
        tmp_links_json = tempfile.NamedTemporaryFile(suffix=".json").name
        tmp_layout_json = tempfile.NamedTemporaryFile(suffix=".json").name

        # Prepare data for Circos plot
        circos_prep(blast_output, tmp_links_json, tmp_layout_json)

        # Load the prepared JSON data
        try:
            with open(tmp_links_json) as f:
                circos_graph_data = json.load(f)
            with open(tmp_layout_json) as f:
                circos_layout = json.load(f)
        except Exception as e:
            logging.error(f"Error loading JSON data: {e}")
            return html.Div(["Error loading plot data."])

        # Check if the loaded data is not empty
        if not circos_graph_data or not circos_layout:
            return html.Div(["No valid data found for the BLAST search."])

        logging.info("Circos graph data:", circos_graph_data)
        logging.info("Circos layout data:", circos_layout)

        layout_config = {
            "innerRadius": 100,
            "outerRadius": 200,
            "cornerRadius": 4,
            "labels": {
                "size": 10,
                "color": "#4d4d4d",
            },
        }

        # Minimal Circos plot configuration
        try:
            circos_plot = dashbio.Circos(
                layout=circos_layout,
                config=layout_config,
                tracks=[
                    {
                        "type": "CHORDS",
                        "data": circos_graph_data,
                        "config": {
                            "opacity": 0.7,
                            "color": {
                                "name": "source",
                                "field": "id",
                            },
                            "tooltipContent": {
                                "source": "source",
                                "sourceID": "id",
                                "target": "target",
                                "targetID": "id",
                                "targetEnd": "end",
                                "fields": [
                                    {
                                        "field": "value.pident",
                                        "name": "Identity (%)",
                                    },
                                    {"field": "value.evalue", "name": "E-value"},
                                    {
                                        "field": "value.bitscore",
                                        "name": "Bit Score",
                                    },
                                ],
                            },
                        },
                    },
                ],
                id="blast-chord",
            )

            return html.Div(circos_plot)

        except Exception as e:
            logging.error(f"Error creating Circos plot: {e}")
            return html.Div(["Error creating plot."])
    else:
        return html.Div(["No results found for the BLAST search."])


# TODO: link in ship classification information for subjects here
def blast_table(ship_blast_results):
    tbl = html.Div(
        [
            dash_table.DataTable(
                columns=[
                    {
                        "name": i,
                        "id": i,
                        "deletable": False,
                        "selectable": True,
                    }
                    for i in ship_blast_results.columns
                ],
                data=ship_blast_results.to_dict("records"),
                hidden_columns=["qseqid", "qseq", "sseq"],
                id="ship-blast-table",
                editable=False,
                sort_action="native",
                sort_by=[{"column_id": "evalue", "direction": "asc"}],
                sort_mode="single",
                row_selectable="single",
                selected_rows=[0],
                row_deletable=False,
                selected_columns=[],
                page_action="native",
                page_current=0,
                page_size=10,
                export_format="tsv",
                css=[{"selector": ".show-hide", "rule": "display: none"}],
                style_table={
                    "overflow": "hidden",
                    "overflowX": "auto",
                    "maxWidth": "100%",
                    "padding": "10px",
                },
                style_cell={
                    "whiteSpace": "normal",
                },
            ),
            dbc.Button(
                "Download BLAST results",
                id="blast-dl-button",
                n_clicks=0,
                style={"textAlign": "center", "fontSize": "1rem"},
                className="d-grid gap-2 col-3 mx-auto",
            ),
            dcc.Download(id="blast-dl"),
        ]
    )
    return tbl


def select_ship_family(hmmer_results):
    hmmer_results["evalue"] = pd.to_numeric(hmmer_results["evalue"], errors="coerce")
    hmmer_results.dropna(subset=["evalue"], inplace=True)

    if hmmer_results.empty:
        logging.warning("HMMER results DataFrame is empty after dropping NaNs.")
        return None, None, None

    try:
        hmmer_results_sorted = hmmer_results.sort_values(by=["query_id", "evalue"])
        best_matches = hmmer_results_sorted.drop_duplicates(
            subset="query_id", keep="first"
        )

        if not best_matches.empty:
            best_match = best_matches.iloc[0]
            superfamily = best_match["hit_IDs"]
            aln_length = best_match["aln_length"]
            evalue = best_match["evalue"]

            return superfamily, aln_length, evalue

        logging.warning("No valid rows found in hmmer_results DataFrame.")
        return None, None, None

    except KeyError as e:
        logging.error(f"KeyError encountered: {e}")
        return None, None, None
    except IndexError as e:
        logging.error(f"IndexError encountered: {e}")
        return None, None, None
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        return None, None, None


def run_lastz(query_type, seqs, output_file):
    """
    Runs LASTZ to align two sequences and writes the output to a specified file.
    """
    lastdb_output = tempfile.NamedTemporaryFile(suffix=".db", delete=True)
    if query_type == "nucl":
        db_command = f"lastdb {lastdb_output.name} {seqs}"
    elif query_type == "prot":
        db_command = f"lastdb -p -c {lastdb_output.name} {seqs}"

    else:
        db_command = None

    if db_command is not None:
        subprocess.run(db_command, shell=True, check=True)
        command = f"lastal {lastdb_output.name} {seqs} -f BLASTTAB > {output_file}"
        subprocess.run(command, shell=True, check=True)


def parse_lastz_output(output_file):
    """
    Parses the LASTZ output file and returns a DataFrame with relevant data.
    """
    columns = [
        "query_id",
        "subject_id",
        "pident",
        "aln_len",
        "mismatch",
        "gap_opens",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]
    df = pd.read_csv(output_file, sep="\t", header=0, comment="#", names=columns)
    df["qstart"] = pd.to_numeric(df["qstart"], errors="coerce")
    df["sstart"] = pd.to_numeric(df["sstart"], errors="coerce")
    df["qend"] = pd.to_numeric(df["qend"], errors="coerce")
    df["send"] = pd.to_numeric(df["send"], errors="coerce")

    return df


def diamond(query, db, threads=2):
    threads = 2
    diamond_cmd = f"/usr/bin/diamond blastp --db {db} -q {query} -f 6 qseqid pident evalue qseq -e 0.001 --strand both -p {threads} -k 1 --skip-missing-seqids"
    subprocess.run(diamond_cmd, shell=True, check=True)
