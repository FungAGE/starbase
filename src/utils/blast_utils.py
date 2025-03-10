import dash_mantine_components as dmc
from dash import dcc, html
import dash_bio as dashbio
import dash_table

import os
import tempfile
import subprocess
import json
import pandas as pd

from Bio import SearchIO

from src.utils.seq_utils import get_protein_sequence, parse_fasta_from_text, clean_shipID
from src.database.blastdb import blast_db_exists, create_dbs
from src.components.tables import create_ag_grid
from src.components.error_boundary import create_error_alert
import logging

logger = logging.getLogger(__name__)


def print_table(list_thing, output_name, header=""):
    with open(output_name, "w") as ofile:
        if header != "":
            ofile.write(header + "\n")
        for line in list_thing:
            hitstring = "\t".join(str(tab) for tab in line)
            ofile.write(hitstring + "\n")


def stitch_blast(tabfile, output_name):
    # ------------------------------------------------------
    # Read file into list
    # ------------------------------------------------------
    with open(tabfile, "r") as tabopen:
        unsorted_tabs = [line.rstrip("\n").split("\t") for line in tabopen]

    # ------------------------------------------------------
    # Sort file in case some subjects have multiple hits in different places
    # ------------------------------------------------------
    # Make a list of all unique queries
    queries = list(set(tab[0] for tab in unsorted_tabs))

    # tabs = [] # the future sorted list
    THRESHOLD = (
        2500  # A threshold for the distance between two hits to be distinct loci
    )
    stitchedtab = []

    # For every query, sort locally
    for query in queries:
        subject_dic = {}
        right_subject_order = []  # To keep the right order
        for tab in unsorted_tabs:
            if query == tab[0]:  # Just for this query
                subject = tab[1]
                if subject not in subject_dic.keys():
                    subject_dic[subject] = [tab]
                    right_subject_order.extend([subject])
                else:
                    subject_dic[subject].append(tab)

        # Sort the resulting list for that particular query, check every subject and sort it locally
        for sub in right_subject_order:
            currentsubject_sorted = sorted(
                subject_dic[sub], key=lambda x: int(x[9])
            )  # Sort by the subject_start
            # tabs.extend(currentsubject_sorted) # Save it into the final tabs list

            # ------------------------------------------------------
            # Stitch pieces together
            # ------------------------------------------------------
            # print("Subject:", sub)

            maxalign = 0
            queryStarts = []
            queryEnds = []
            subjectStarts = []
            subjectEnds = []

            for i in range(0, len(currentsubject_sorted)):  # Find the main piece
                # print(i, len(currentsubject_sorted), currentsubject_sorted[i] # For debugging)
                # query_id, subject_id, percent_identity, alignment_length, N_mismatches, N_gaps, query_start, query_end, subject_start, subject_end, evalue, bit_score = currentsubject_sorted[i]
                query_start = int(currentsubject_sorted[i][6])
                query_end = int(currentsubject_sorted[i][7])
                query_seq = currentsubject_sorted[i][12]

                alignment_length = int(currentsubject_sorted[i][3])
                subject_start = int(currentsubject_sorted[i][8])
                subject_end = int(currentsubject_sorted[i][9])
                subject_seq = currentsubject_sorted[i][13]

                queryStarts.append(query_start)
                queryEnds.append(query_end)
                subjectStarts.append(subject_start)
                subjectEnds.append(subject_end)

                # Is this last piece one better than the previous ones?
                if alignment_length > maxalign:
                    upper_hit = currentsubject_sorted[i]
                    maxalign = alignment_length

                # There is only one clean hit
                if len(currentsubject_sorted) == 1:
                    stitchedtab.append(currentsubject_sorted[i])

                # The hit is broken
                elif i < len(currentsubject_sorted) - 1:
                    next_subject_start = int(currentsubject_sorted[i + 1][8])

                    if (
                        abs(subject_end - next_subject_start) > THRESHOLD
                    ):  # It's probably not part of the same hit # It was subject_start - next_subject_start before
                        # So write down the previous one
                        # -----------------
                        # The new values for the subject
                        # -----------------
                        # query_id, subject_id, percent_identity, alignment_length, N_mismatches, N_gaps, query_start, query_end, subject_start, subject_end, evalue, bit_score
                        new_tab = upper_hit[0:4]

                        # These are our new values of the "unbroken" hit, but I'm not sure how to retrieve the no. of gaps and mistmatches
                        new_tab.extend([".", ".", min(queryStarts), max(queryEnds)])

                        # subject_start > subject_end for the upper_hit
                        if int(upper_hit[8]) > int(upper_hit[9]):  # hit is reversed
                            new_tab.extend([max(subjectStarts), min(subjectEnds)])
                        else:
                            new_tab.extend([min(subjectStarts), max(subjectEnds)])

                        # Let's leave the e-val and Bit score the same as the upper hit
                        new_tab.extend(upper_hit[10:11])

                        gap_length = abs(next_subject_start - subject_end)

                        # Concatenate subject sequences with '-' characters filling the gap
                        stitched_subject_seq = (
                            currentsubject_sorted[i][13]
                            + "-" * gap_length
                            + currentsubject_sorted[i + 1][13]
                        )

                        # Update the subject_seq in the upper_hit
                        upper_hit[13] = stitched_subject_seq

                        stitchedtab.append(new_tab)  # Write it in the final output

                        # -----------------
                        # Reset for the next hit
                        # -----------------
                        queryStarts = []
                        queryEnds = []
                        subjectStarts = []
                        subjectEnds = []

                        maxalign = int(currentsubject_sorted[i + 1][3])
                        upper_hit = currentsubject_sorted[i + 1]

                else:
                    # The last hit in that subject
                    # -----------------
                    # The new values for the subject
                    # -----------------
                    new_tab = upper_hit[0:4]

                    # These are our new values of the "unbroken" hit, but I'm not sure how to retrieve the no. of gaps and mistmatches
                    new_tab.extend([".", ".", min(queryStarts), max(queryEnds)])

                    # subject_start > subject_end for the upper_hit
                    if int(upper_hit[8]) > int(upper_hit[9]):  # hit is reversed
                        new_tab.extend([max(subjectStarts), min(subjectEnds)])
                    else:
                        new_tab.extend([min(subjectStarts), max(subjectEnds)])
                    # Let's leave the e-val and Bit score the same as the upper hit
                    new_tab.extend(upper_hit[10:])

                    stitchedtab.append(new_tab)  # Write it in the final output

            # print # For separating the subjects

    # ------------------------------------------------------
    # Print filtered tab file
    # ------------------------------------------------------
    print_table(stitchedtab, output_name)
    df = pd.read_csv(
        output_name,
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
    df = df.dropna()

    logger.info(f"BLAST results parsed with {len(df)} hits.")
    return df


def run_blast(db_list, query_type, query_fasta, tmp_blast, input_eval=0.01, threads=2):
    try:
        # Add input size check
        max_input_size = 10 * 1024 * 1024  # 10MB
        if os.path.getsize(query_fasta) > max_input_size:
            logger.error(f"Input FASTA file too large: {os.path.getsize(query_fasta)} bytes")
            return None
        logger.debug(f"db_list contents: {db_list}")
        
        if not isinstance(db_list, dict):
            logger.error(f"db_list must be a dictionary, got {type(db_list)}")
            raise ValueError("Invalid database configuration")
            
        ship_config = db_list.get('ship')
        if ship_config is None:
            logger.error("'ship' key not found in db_list")
            raise ValueError("Database path for 'ship' not configured")
            
        db_path = ship_config.get(query_type)
        if db_path is None:
            logger.error(f"No database path found for query type: {query_type}")
            raise ValueError(f"Database path for {query_type} not configured")
            
        logger.debug(f"Using database path: {db_path}")
        
        if not blast_db_exists(db_path):
            logger.info("BLAST database not found. Creating new database...")
            create_dbs()
            
            if not blast_db_exists(db_path):
                logger.error("Failed to create BLAST database")
                raise ValueError("Failed to create BLAST database")
            
        if not isinstance(db_path, (str, bytes, os.PathLike)):
            logger.error(f"db_path must be a path-like object, got {type(db_path)}")
            raise ValueError("Invalid database path type")
            
        blast_cmd = [
            "blastn" if query_type == "nucl" else "blastp",
            "-query", str(query_fasta),
            "-db", str(db_path),
            "-out", str(tmp_blast),
            "-evalue", str(input_eval),
            "-num_threads", str(threads),
            "-max_target_seqs", "100",
            "-max_hsps", "1",
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
        ]
        
        subprocess.run(blast_cmd, check=True, timeout=300)  # Add 5-minute timeout
        
        blast_results = pd.read_csv(
            tmp_blast,
            sep="\t",
            names=[
                "qseqid", "sseqid", "pident", "length", "mismatch",
                "gapopen", "qstart", "qend", "sstart", "send",
                "evalue", "bitscore", "qseq", "sseq"
            ]
        )
        
        return blast_results

    except subprocess.TimeoutExpired:
        logger.error("BLAST search timed out after 5 minutes")
        return None
    except Exception as e:
        logger.error(f"Error during BLAST search: {e}")
        return None

def hmmsearch(
    db_list=None,
    query_type=None,
    input_genes="tyr",
    input_eval=None,
    query_fasta=None,
    threads=None,
):
    tmp_hmmer = tempfile.NamedTemporaryFile(suffix=".hmmer.txt").name
    hmmer_db = db_list["gene"][input_genes]["hmm"][query_type]
    if not os.path.exists(hmmer_db) or os.path.getsize(hmmer_db) == 0:
        raise ValueError(f"HMMER database {hmmer_db} not found or is empty.")

    # Calculate optimal thread count based on system resources
    if threads is None:
        # Use 1 thread by default to avoid resource contention
        # since multiple workers/threads might run HMMER simultaneously
        threads = 1
    else:
        # Cap threads to avoid oversubscription
        # Maximum 4 threads per HMMER process to leave resources for other concurrent requests
        threads = min(int(threads), 4)

    hmmer_cmd = f"hmmsearch -o {tmp_hmmer} --cpu {threads} --domE {input_eval} {hmmer_db} {query_fasta}"
    logger.info(f"Running HMMER search: {hmmer_cmd}")
    subprocess.run(hmmer_cmd, shell=True)
    return tmp_hmmer


def run_hmmer(
    db_list=None,
    query_type=None,
    input_genes="tyr",
    input_eval=None,
    query_fasta=None,
    threads=None,
):
    try:
        # family classification should only use protein hmm
        # first run hmmer using nucl hmm
        tmp_hmmer = hmmsearch(
            db_list,
            query_type,
            input_genes,
            input_eval,
            query_fasta,
            threads,
        )

        tmp_hmmer_parsed = parse_hmmer(tmp_hmmer)
        # extract the gene sequence
        gene_header, gene_seq = extract_gene_from_hmmer(tmp_hmmer_parsed)

        # translate nucl queries
        if query_type == "nucl" and gene_header is not None and gene_seq is not None:
            tmp_protein = get_protein_sequence(gene_header, gene_seq)
            if tmp_protein is not None:
                # run hmmer using protein sequence of extracted gene
                tmp_hmmer_protein = hmmsearch(
                    db_list,
                    "prot",
                    input_genes,
                    input_eval,
                    tmp_protein,
                    threads,
                )
                tmp_hmmer_protein_parsed = parse_hmmer(tmp_hmmer_protein)
                logger.info(
                    f"Second hmmersearch run stored at: {tmp_hmmer_protein_parsed}"
                )
                hmmer_results = pd.read_csv(tmp_hmmer_protein_parsed, sep="\t")
        else:
            hmmer_results = pd.read_csv(tmp_hmmer_parsed, sep="\t")

        logger.debug(f"HMMER results parsed: {hmmer_results.shape[0]} rows.")

        return hmmer_results.to_dict("records")
    except Exception as e:
        logger.error(f"Error in HMMER search: {e}")
        return None


# Parse the HMMER results
def parse_hmmer(hmmer_output_file):
    parsed_file = tempfile.NamedTemporaryFile(suffix=".hmmer.parsed.txt").name
    with open(parsed_file, "w") as tsv_file:
        tsv_file.write(
            "query_id\thit_IDs\taln_length\tquery_start\tquery_end\tgaps\tquery_seq\tsubject_seq\tevalue\tbitscore\n"
        )
        for record in SearchIO.parse(hmmer_output_file, "hmmer3-text"):
            for hit in record.hits:
                for hsp in hit.hsps:
                    query_seq = str(hsp.query.seq)
                    subject_seq = clean_shipID(str(hsp.hit.seq))
                    aln_length = hsp.aln_span
                    query_start = hsp.query_start
                    query_end = hsp.query_end
                    gaps = str("N/A")
                    bitscore = hsp.bitscore
                    evalue = hsp.evalue
                    tsv_file.write(
                        f"{hit.id}\t{record.id}\t{aln_length}\t{query_start}\t{query_end}\t{gaps}\t{query_seq}\t{subject_seq}\t{evalue}\t{bitscore}\n"
                    )
    return parsed_file


def extract_gene_from_hmmer(parsed_file):
    data = pd.read_csv(parsed_file, sep="\t")

    # Get rows with the lowest e-value for each unique entry in Query
    min_evalue_rows = data.loc[data.groupby("query_id")["evalue"].idxmin()]
    # Reset the index
    min_evalue_rows = min_evalue_rows.reset_index(drop=True)

    # top_hit_out_path = tempfile.NamedTemporaryFile(suffix=".besthit.txt").name
    # top_hit_out_path = tempfile.NamedTemporaryFile(suffix=".best_hsp.fa").name
    # logger.info(f"Best hit for gene sequence: {top_hit_out_path}")

    query = min_evalue_rows.loc[0, "query_id"]
    qseq = min_evalue_rows.loc[0, "query_seq"].replace(".", "")

    logger.info(f"Sequence has length {len(qseq)}")

    return query, qseq

    # output = html.Div(
    #     [
    #         dbc.Button(
    #             "Download best captain hit",
    #             id="subject-seq-button",
    #             n_clicks=0,
    #             className="d-grid gap-2 col-6 mx-auto",
    #             style={"fontSize": "1rem"},
    #         ),
    #         dcc.Download(id="subject-seq-dl-package"),
    #     ]
    # )
    # return output


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
            logger.error(f"Error loading JSON data: {e}")
            return html.Div(["Error loading plot data."])

        # Check if the loaded data is not empty
        if not circos_graph_data or not circos_layout:
            return html.Div(["No valid data found for the BLAST search."])

        logger.info("Circos graph data:", circos_graph_data)
        logger.info("Circos layout data:", circos_layout)

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
            logger.error(f"Error creating Circos plot: {e}")
            return html.Div(["Error creating plot."])
    else:
        return html.Div(["No results found for the BLAST search."])


def blast_table(ship_blast_results):
    """Creates an AG Grid table for displaying BLAST results."""
    try:
        # Validation checks...
        if not isinstance(ship_blast_results, pd.DataFrame):
            logger.error("Invalid input type for blast_table")
            return html.Div("Error: Invalid data format")
            
        if ship_blast_results.empty:
            logger.warning("Empty DataFrame passed to blast_table")
            return html.Div("No results to display")
            
        # Column checks...
        required_cols = ["sseqid", "pident", "length", "evalue", "bitscore"]
        missing_cols = [col for col in required_cols if col not in ship_blast_results.columns]
        if missing_cols:
            logger.error(f"Missing required columns: {missing_cols}")
            return html.Div(f"Error: Missing columns: {', '.join(missing_cols)}")

        columns = [
            {
                "field": "accession_tag",
                "headerName": "Accession",
                "flex": 1,
                "cellStyle": {"cursor": "pointer", "color": "#1976d2"},
                "tooltipField": "accession_tag"
            },
            {
                "field": "familyName",
                "headerName": "Starship Family",
                "flex": 1,
                "tooltipField": "familyName"
            },
            {
                "field": "pident",
                "headerName": "Percent Identity",
                "flex": 1,
                "valueFormatter": {"function": "value.toFixed(2)"},
                "type": "numericColumn",
                "filter": "agNumberColumnFilter"
            },
            {
                "field": "length",
                "headerName": "Hit Length",
                "flex": 1,
                "type": "numericColumn",
                "filter": "agNumberColumnFilter"
            },
            {
                "field": "evalue",
                "headerName": "E-value",
                "flex": 1,
                "valueFormatter": {"function": "value.toExponential(2)"},
                "type": "numericColumn",
                "filter": "agNumberColumnFilter"
            },
            {
                "field": "bitscore",
                "headerName": "Bitscore",
                "flex": 1,
                "valueFormatter": {"function": "value.toFixed(2)"},
                "type": "numericColumn",
                "filter": "agNumberColumnFilter"
            }
        ]

        return html.Div(
            html.Div(
                create_ag_grid(
                    df=ship_blast_results,
                    id="blast-table",
                    columns=columns,
                    select_rows=False,
                    pg_sz=15,
                ),
                style={
                    "height": "400px",
                    "width": "100%",
                    "overflow": "auto",  # Changed from 'hidden' to 'auto'
                    "position": "relative"
                }
            ),
            style={
                "marginBottom": "20px",
                "width": "100%"
            }
        )
        
    except Exception as e:
        logger.error(f"Error in blast_table: {e}")
        return html.Div(f"Error creating table: {str(e)}")

def blast_download_button():
    """Creates the download button for BLAST results."""
    return html.Div([
        dmc.Space(h="xl"),
        dmc.Center(
            dmc.Button(
                "Download BLAST Results",
                id="blast-dl-button",
                variant="gradient",
                gradient={"from": "indigo", "to": "cyan"},
                size="lg",
                leftSection=[html.I(className="bi bi-download")]
            )
        ),
        dcc.Download(id="blast-dl")
    ])

def select_ship_family(hmmer_results):
    hmmer_results["evalue"] = pd.to_numeric(hmmer_results["evalue"], errors="coerce")
    hmmer_results.dropna(subset=["evalue"], inplace=True)

    if hmmer_results.empty:
        logger.warning("HMMER results DataFrame is empty after dropping NaNs.")
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
            try:
                evalue_num = float(best_match["evalue"])
                evalue = "{:.2e}".format(evalue_num)
            except (ValueError, TypeError, KeyError) as e:
                logger.error(f"Error formatting e-value: {e}")
                evalue = str(best_match.get("evalue", "N/A"))

            return superfamily, aln_length, evalue

        logger.warning("No valid rows found in hmmer_results DataFrame.")
        return None, None, None

    except KeyError as e:
        logger.error(f"KeyError encountered: {e}")
        return None, None, None
    except IndexError as e:
        logger.error(f"IndexError encountered: {e}")
        return None, None, None
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
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


def run_diamond(
    db_list=None,
    query_type=None,
    input_genes="tyr",
    input_eval=None,
    query_fasta=None,
    threads=2,
):

    diamond_out = tempfile.NamedTemporaryFile(suffix=".fa").name

    header, seq = parse_fasta_from_text(query_fasta)

    diamond_db = db_list["gene"][input_genes]["prot"]
    if not os.path.exists(diamond_db) or os.path.getsize(diamond_db) == 0:
        raise ValueError(f"HMMER database {diamond_db} not found or is empty.")

    if query_type == "nucl":
        blast_type = "blastx"
        out_fmt = "6 qseqid sseqid length qstart qend gaps qseq_translated sseq evalue bitscore"
    else:
        blast_type = "blastp"
        out_fmt = "6 qseqid sseqid length qstart qend gaps qseq sseq evalue bitscore"

    subprocess.run("/home/adrian/anaconda3/bin/diamond help", shell=True, check=True)
    diamond_cmd = f"diamond {blast_type} --db {diamond_db} -q {query_fasta} -f {out_fmt} -e 0.001 --strand both -p {threads} -k 1 --skip-missing-seqids | sed '1i >{header}' > {diamond_out}"

    subprocess.run(diamond_cmd, shell=True, check=True)

    column_names = out_fmt.split()[1:]
    diamond_results = pd.read_csv(diamond_out, sep="\t", names=column_names)

    return diamond_results.to_dict("records")

def make_captain_alert(family, aln_length, evalue, search_type="blast"):
    try:
        # Validate inputs
        if not family or not isinstance(family, str):
            logger.error(f"Invalid family parameter: {family}")
            return create_error_alert("Invalid family name")
            
        if not aln_length:
            logger.error(f"Invalid alignment length: {aln_length}")
            return create_error_alert("Invalid alignment length")
            
        try:
            evalue_num = float(evalue)
            if evalue_num == 0.0:
                formatted_evalue = "0.00e+00"
            else:
                formatted_evalue = "{:.2e}".format(evalue_num)
        except (ValueError, TypeError) as e:
            logger.error(f"Error formatting e-value {evalue}: {e}")
            return create_error_alert("Invalid e-value format")
            
        if search_type not in ["blast", "hmmsearch"]:
            logger.error(f"Invalid search type: {search_type}")
            return create_error_alert("Invalid search type")

        if search_type == "blast":
            return dmc.Alert(
                title="Family Classification via BLAST Search",
                children=[
                    f"Your sequence is likely in Starship family: {family}",
                    dmc.Space(h=5),
                    dmc.Text(
                        f"BLAST Search: Alignment length = {aln_length}, E-value = {formatted_evalue}",
                        size="sm",
                        c="dimmed"
                    ),
                ],
                color="blue",
                variant="light",
                withCloseButton=False,
            )
        else:  # hmmsearch
            return dmc.Alert(
                title="Family Classification via HMMER Search",
                children=[
                    f"Your sequence is likely in Starship family: {family}",
                    dmc.Space(h=5),
                    dmc.Text(
                        f"HMMER Search: Alignment length = {aln_length}, E-value = {formatted_evalue}",
                        size="sm",
                        c="dimmed"
                    ),
                ],
                color="blue",
                variant="light",
                withCloseButton=False,
            )
            
    except Exception as e:
        logger.error(f"Error in make_captain_alert: {e}")
        return create_error_alert(str(e))

def process_captain_results(captain_results_dict):
    no_captain_alert = dmc.Alert(
        "No captain sequence found (e-value threshold 0.01).",
        color="warning",
    )

    if not captain_results_dict:
        return no_captain_alert
        
    try:
        captain_results_df = pd.DataFrame(captain_results_dict)
        if len(captain_results_df) > 0:
            superfamily, family_aln_length, family_evalue = select_ship_family(captain_results_df)
            if superfamily:
                return make_captain_alert(superfamily, family_aln_length, family_evalue, search_type="hmmsearch")
                
        return no_captain_alert
    except Exception as e:
        logger.error(f"Error processing captain results: {str(e)}")
        return dmc.Alert(
            title="Error",
            children="Failed to process captain results",
            color="red",
            variant="light",
        )
    
def create_no_matches_alert():
    return dmc.Alert(
        title="No Matches Found",
        children=[
            dmc.Text("Your sequence did not match any Starships in our database."),
            dmc.Space(h=10),
            dmc.Text("Suggestions:", size="sm", fw=500),
            dmc.Space(h=5),
            dmc.List(
                [
                    dmc.ListItem(
                        "Check if your sequence is in the correct format",
                        size="sm"
                    ),
                    dmc.ListItem(
                        "Try searching with a different region of your sequence",
                        size="sm"
                    ),
                    dmc.ListItem(
                        "Consider using a less stringent E-value threshold",
                        size="sm"
                    ),
                ],
                withPadding=True,
                spacing="xs",
                size="sm",
                type="ordered"
            ),
        ],
        color="yellow",
        variant="light",
        withCloseButton=False,
    )