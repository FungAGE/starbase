from dash import callback
from dash.dependencies import Output, Input

import re
import tempfile
import subprocess
import pandas as pd

import plotly.graph_objects as go


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


@callback(
    Output("lastz-plot", "figure"),
    [
        Input("ship-blast-table", "derived_virtual_data"),
        Input("ship-blast-table", "derived_virtual_selected_rows"),
    ],
)
def create_alignment_plot(ship_blast_results, selected_row):
    """
    Creates a Plotly scatter plot from the alignment DataFrame and saves it as a PNG.
    """
    tmp_fasta_clean = tempfile.NamedTemporaryFile(suffix=".fa", delete=True)
    lastz_output = tempfile.NamedTemporaryFile(suffix=".tsv", delete=True)

    # Convert ship_blast_results to a Pandas DataFrame
    ship_blast_results_df = pd.DataFrame(ship_blast_results)

    if selected_row is None:
        return None
    else:
        row = ship_blast_results_df.iloc[selected_row]
        qseq = re.sub("-", "", row["qseq"].iloc[0])
        qseqid = row["qseqid"].iloc[0]
        sseq = re.sub("-", "", row["sseq"].iloc[0])
        sseqid = row["sseqid"].iloc[0]

        # Write the specific row to the file
        with open(tmp_fasta_clean.name, "w") as f:
            f.write(f">{qseqid}\n")
            f.write(f"{qseq}\n")
            f.write(f">{sseqid}\n")
            f.write(f"{sseq}\n")

        # Run LASTZ alignment
        run_lastz(tmp_fasta_clean.name, lastz_output.name)

        # Parse LASTZ output
        lastz_df = parse_lastz_output(lastz_output.name)

        # Extract individual columns and prepare the data for plotting
        x_values = []
        y_values = []

        for _, row in lastz_df.iterrows():
            x_values.append(row["qstart"])
            x_values.append(row["qend"])
            y_values.append(row["sstart"])
            y_values.append(row["send"])

        # Create the scatter plot using Plotly Graph Objects
        fig = go.Figure(data=go.Scatter(x=x_values, y=y_values, mode="lines"))

        # Add layout details
        fig.update_layout(
            title="LASTZ Alignment",
            xaxis_title=qseqid,
            yaxis_title=sseqid,
        )

        return fig
