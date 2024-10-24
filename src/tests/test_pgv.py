import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import dcc, html, callback
from dash.dependencies import Output, Input, State

import tempfile
import pandas as pd

from sqlalchemy import create_engine

import warnings

warnings.filterwarnings("ignore")

import os

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Gff
from matplotlib.lines import Line2D
from pygenomeviz.utils import ColorCycler
from pygenomeviz.align import Blast, AlignCoord
from jinja2 import Template

# Initialize the SQLite engine
engine = create_engine("sqlite:///src/data/db/starbase.sqlite")

query = """
SELECT a.accession_tag, f.familyName, t.species
FROM joined_ships j
JOIN taxonomy t ON j.taxid = t.id
JOIN family_names f ON j.ship_family_id = f.id
JOIN accessions a ON j.ship_id = a.id
JOIN ships s on s.accession = a.id
JOIN gff g ON a.id = g.ship_id
WHERE s.sequence is NOT NULL AND g.ship_id is NOT NULL AND j.orphan IS NULL
"""
table_df = pd.read_sql_query(query, engine)
table_df = table_df.drop_duplicates(subset=["accession_tag"])

ColorCycler.set_cmap("tab10")

layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        dmc.Grid(
            justify="center",
            align="top",
            children=[
                dmc.GridCol(
                    style={"padding": "40px"},
                    span={
                        "sm": 12,
                        "lg": 8,
                    },
                    children=[
                        dcc.Loading(
                            id="loading",
                            type="circle",
                            children=[
                                html.Div(id="pgv-table", className="center-content"),
                                dbc.Button(
                                    "Show Starship(s) in Viewer",
                                    id="update-button",
                                    n_clicks=0,
                                    className="d-grid gap-2 mx-auto",
                                    style={"fontSize": "1rem"},
                                ),
                            ],
                        ),
                    ],
                ),
                dmc.GridCol(
                    style={"padding": "40px"},
                    span={
                        "sm": 12,
                        "lg": 8,
                    },
                    children=[
                        dcc.Loading(
                            id="loading-1",
                            type="circle",
                            children=[
                                html.Div(
                                    id="pgv-figure",
                                    style={
                                        "height": "800px",
                                        "width": "100%",
                                        "overflow": "auto",
                                        "border": "1px solid #ccc",
                                    },
                                ),
                            ],
                        )
                    ],
                ),
            ],
        ),
    ],
)


def plot_legend(gv):
    fig = gv.plotfig(fast_render=False)

    # Plot legend for groups
    _ = fig.legend(
        handles=[
            Line2D(
                [], [], marker=">", color="#dc267fff", label="captain", ms=12, ls="none"
            ),
            Line2D(
                [],
                [],
                marker=">",
                color="#785ef0ff",
                label="auxillary",
                ms=12,
                ls="none",
            ),
            Line2D(
                [], [], marker=">", color="#ffb000ff", label="Others", ms=12, ls="none"
            ),
        ],
        fontsize=12,
        title="Starship Gene",
        title_fontsize=12,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        handlelength=1.0,
    )
    return fig


def add_gene_feature(gene, track, idx=None):

    start = int(gene.location.start)
    end = int(gene.location.end)
    strand = gene.location.strand
    gene_name = str(gene.qualifiers.get("Alias", [""])[0])
    attributes = gene.qualifiers

    # BUG: there is still an issue with SeqFeature end coordinates being less than start coordinates
    if (end < start) and idx is not None:
        gene.reverse_complement(id=True, name=True, description=True)
        # start, end = end, start

    if "tyr" in gene_name:
        color = "#dc267fff"
    elif any(substring in gene_name for substring in ["fre", "DUF3723", "nlr", "plp"]):
        color = "#785ef0ff"
    else:
        color = "#ffb000ff"

    track.add_feature(
        start=start,
        end=end,
        strand=strand,
        lw=1,
        color=color,
        label_type="gene",
        extra_tooltip=attributes,
    )

    return track


def inject_svg_to_html(svg_file, html_template_file, output_html_file):
    """Inject SVG content into an HTML template."""

    # Read the SVG content from file
    with open(svg_file, "r") as svg:
        svg_content = svg.read()

    # Read the HTML template file
    with open(html_template_file, "r") as template_file:
        template_content = template_file.read()

    # Use jinja2 or string replacement to inject the SVG into the HTML template
    template = Template(template_content)
    rendered_html = template.render(svg_content=svg_content)

    # Save the rendered HTML to an output file
    with open(output_html_file, "w") as output_file:
        output_file.write(rendered_html)


def write_tmp(df, seqid, type=None):
    tmp = tempfile.NamedTemporaryFile(suffix=f".{type}", delete=False).name
    if type == "gff":
        df.to_csv(tmp, sep="\t", header=False, index=False)
    elif type == "fa":
        seq = df["sequence"][0]
        with open(tmp, "w") as f:
            f.write(f">{seqid}\n{seq}\n")
    return tmp


def fetch_gff(accession):
    query = """
    SELECT g.*
    FROM gff g
    JOIN accessions a ON g.ship_id = a.id
    JOIN ships s on s.accession = a.id
    JOIN joined_ships j ON j.ship_id = a.id
    JOIN taxonomy t ON j.taxid = t.id
    JOIN family_names f ON j.ship_family_id = f.id
    WHERE a.accession_tag = :accession_tag AND j.orphan IS NULL
    """

    # print("Executing GFF query with accession:", accession)
    df = pd.read_sql_query(query, engine, params={"accession_tag": accession})

    if df.empty:
        print(f"No GFF records found for accession_tag: {accession}")

    return df


def fetch_fa(accession):
    query = """
    SELECT s.*
    FROM ships s
    LEFT JOIN accessions a ON s.accession = a.id
    WHERE a.accession_tag = :accession_tag
    """

    # print("Executing FA query with accession:", accession)
    df = pd.read_sql_query(query, engine, params={"accession_tag": accession})

    if df.empty:
        print(f"No FA records found for accession_tag: {accession}")

    return df


def single_pgv(gff_file, tmp_file):
    gff = Gff(gff_file)

    gv = GenomeViz()
    gv.set_scale_xticks()
    # gv.set_scale_bar(ymargin=0.5)

    for seqid, size in gff.get_seqid2size().items():
        track = gv.add_feature_track(seqid, size, labelsize=15)
        track.add_sublabel(size=10, color="grey")
        gene_features = gff.get_seqid2features(feature_type="gene")[seqid]
        for gene in gene_features:
            add_gene_feature(gene, track)
    fig = plot_legend(gv)
    gv.savefig_html(tmp_file, fig)


def multi_pgv(gff_files, seqs, tmp_file):
    gff_list = list(map(Gff, gff_files))
    gv = GenomeViz(track_align_type="center", fig_track_height=0.7)
    gv.set_scale_bar()

    for gff in gff_list:
        # Check GFF SeqIDs
        gff_seqids = gff.get_seqid2size().keys()
        print(f"GFF SeqIDs: {gff_seqids}")

        for seqid, features in gff.get_seqid2features("gene").items():
            print(f"Processing seqid: {seqid}")

            if seqid not in gff.get_seqid2size():
                print(f"Error: SeqID {seqid} not found in GFF sizes")
                continue

            track = gv.add_feature_track(seqid, gff.get_seqid2size(), align_label=False)
            segment = track.get_segment(seqid)

            for idx, gene in enumerate(features):
                add_gene_feature(gene, segment, idx)

    for seq_file in seqs:
        if not os.path.isfile(seq_file):
            print(f"Error: {seq_file} is not a valid file path.")
            return  # Exit the function if an invalid file is found
        if os.path.getsize(seq_file) == 0:
            print(f"Error: {seq_file} is an empty file.")
            return  # Exit the function if an empty file is found

    blast_cmd = Blast(seqs, seqtype="nucleotide")
    align_coords = blast_cmd.run()
    align_coords = AlignCoord.filter(align_coords, length_thr=50, identity_thr=20)

    print(f"Found {len(align_coords)} alignments")

    if len(align_coords) > 0:
        min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
        print(f"Minimum identity: {min_ident}")

        color, inverted_color = "blue", "orange"
        for ac in align_coords:
            print(f"Adding link between {ac.query_link} and {ac.ref_link}")
            gv.add_link(
                ac.query_link,
                ac.ref_link,
                color=color,
                inverted_color=inverted_color,
                v=ac.identity,
                vmin=min_ident,
            )
        gv.set_colorbar([color, inverted_color], vmin=min_ident)

    fig = plot_legend(gv)
    gv.savefig_html(tmp_file, fig)


@callback(
    Output("pgv-figure", "children"),
    Input("update-button", "n_clicks"),
    [
        State("pgv-table", "derived_virtual_selected_rows"),
        State("pgv-table", "derived_virtual_data"),
    ],
)
def update_pgv(n_clicks, selected_rows, table_data):
    if n_clicks > 0:
        tmp_pgv = tempfile.NamedTemporaryFile(suffix=".html", delete=False).name

        if table_data and selected_rows is not None:
            try:
                table_df = pd.DataFrame(table_data)
                # Debug: Check columns and selected_rows
                # print("Columns in table_df:", table_df.columns.tolist())
                # print("Selected rows:", selected_rows)

                if isinstance(selected_rows, list) and all(
                    isinstance(idx, int) for idx in selected_rows
                ):
                    try:
                        rows = table_df.iloc[selected_rows]
                    except IndexError as e:
                        return html.Div("Index out of bounds")

                    tmp_gffs = []
                    tmp_fas = []
                    for index, row in rows.iterrows():
                        accession = row["accession_tag"]
                        # print("Fetching GFF for accession:", accession)
                        gff_df = fetch_gff(accession)

                        tmp_gff = write_tmp(gff_df, accession, "gff")
                        tmp_gffs.append(tmp_gff)
                        # print("Fetching FA for accession:", accession)
                        fa_df = fetch_fa(accession)
                        tmp_fa = write_tmp(fa_df, accession, "fa")
                        tmp_fas.append(str(tmp_fa))
                        output = html.P("Select up to four Starships to compare.")
                    if len(selected_rows) > 1 and len(selected_rows) <= 4:
                        multi_pgv(tmp_gffs, tmp_fas, tmp_pgv)
                    elif len(selected_rows) == 1:
                        single_pgv(tmp_gffs[0], tmp_pgv)
                    else:
                        output = html.P("No valid selection.")
                    try:
                        with open(tmp_pgv, "r") as file:
                            pgv_content = file.read()
                    except IOError:
                        output = html.P("Failed to read the temporary file.")

                    output = html.Iframe(
                        srcDoc=pgv_content,
                        style={
                            "width": "100%",
                            "height": "100%",
                            "border": "none",
                        },
                    )
                else:
                    output = html.H4("Invalid row selection.")
            except Exception as e:
                print("Exception:", e)
                output = html.H4("Error processing data.")
        else:
            output = html.H4("Select Starship(s) to visualize.")
    else:
        output = html.H4(
            "Select the Starship(s) in the table above and click the button to visualize."
        )
    return html.Div(
        [output],
        className="center-content text-center",
    )


def load_ship_table():
    query = """
    SELECT a.accession_tag, f.familyName, t.species
    FROM joined_ships j
    JOIN taxonomy t ON j.taxid = t.id
    JOIN family_names f ON j.ship_family_id = f.id
    JOIN accessions a ON j.ship_id = a.id
    JOIN ships s on s.accession = a.id
    JOIN gff g ON a.id = g.ship_id
    WHERE s.sequence is NOT NULL AND g.ship_id is NOT NULL AND j.orphan IS NULL
    """
    table_df = pd.read_sql_query(query, engine)
    table_df = table_df.drop_duplicates(subset=["accession_tag"])
    return table_df


def test_pgv(df, ships):
    tmp_pgv = tempfile.NamedTemporaryFile(suffix=".html", delete=False).name
    rows = df[df["accession_tag"].isin(ships)]

    # Call the appropriate function based on number of selections
    if len(ships) > 1:
        if len(ships) > 4:
            raise ValueError("Select up to four Starships to compare.")
        else:
            multi_pgv(gffs, fna, tmp_pgv)
    elif len(ships) == 1:
        gff_file = rows.iloc[0]["gff3"]
        single_pgv(gff_file, tmp_pgv)
    else:
        raise ValueError("No valid selection.")

    # Ensure the file is written and read correctly
    with open(tmp_pgv, "r") as file:
        return html.Div(
            [
                html.Iframe(
                    srcDoc=file.read(),
                    style={
                        "width": "100%",
                        "height": "500px",
                        "border": "none",
                    },
                )
            ]
        )


test_pgv(table_df, ["altals1_s00058"])
