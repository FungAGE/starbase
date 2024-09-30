import warnings

warnings.filterwarnings("ignore")

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import dcc, html, callback
from dash.dependencies import Output, Input, State

import os
import tempfile
import pandas as pd

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Gff
from matplotlib.lines import Line2D
from pygenomeviz.utils import ColorCycler
from pygenomeviz.align import Blast, AlignCoord
from jinja2 import Template

from src.components.sqlite import engine

from src.components.tables import make_ship_table

dash.register_page(__name__)

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
                                html.Div(id="pgv-figure"),
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


def write_tmp(df, type=None):
    tmp = tempfile.NamedTemporaryFile(suffix=f".{type}", delete=False)
    if type == "gff":
        df.to_csv(tmp.name, sep="\t", header=False, index=False)
    elif type == "fa":
        seqid = df["accession_tag"]
        seq = df["ship_sequence"]
        with open(tmp.name, "w") as f:
            f.write(f">{seqid}\n{seq}\n")
    return tmp.name


def fetch_gff(accession):
    query = """
    SELECT g.*
    FROM gff g
    LEFT JOIN accessions a ON g.ship_id = a.id
    LEFT JOIN ships s on s.accession = a.id
    LEFT JOIN joined_ships j ON j.ship_id = a.id
    LEFT JOIN taxonomy t ON j.taxid = t.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    WHERE a.accession_tag = :accession_tag AND j.orphan IS NULL
    """
    df = pd.read_sql_query(query, engine, params={"accession_tag": accession})
    return df


def fetch_fa(accession):
    query = """
    SELECT s.*
    FROM ships s
    LEFT JOIN accessions a ON s.accession = a.id
    WHERE a.accession_tag = :accession_tag
    """
    df = pd.read_sql_query(query, engine, params={"accession_tag": accession})
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

        track = gv.add_feature_track(gff.name, gff.get_seqid2size(), align_label=False)
        for idx, (seqid, features) in enumerate(gff.get_seqid2features("gene").items()):
            segment = track.get_segment(seqid)
            for gene in features:
                add_gene_feature(gene, segment, idx)

    # Run BLAST alignment & filter by user-defined threshold
    align_coords = Blast(seqs, seqtype="protein").run()
    align_coords = AlignCoord.filter(align_coords, length_thr=200, identity_thr=50)

    # Plot BLAST alignment links
    if len(align_coords) > 0:
        min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
        color, inverted_color = "grey", "red"
        for ac in align_coords:
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
    # gv.savefig("tmp/genbank_comparison_by_blast.svg")


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

                if isinstance(selected_rows, list) and all(
                    isinstance(idx, int) for idx in selected_rows
                ):
                    try:
                        rows = table_df.iloc[selected_rows]
                    except IndexError as e:
                        return html.Div("Index out of bounds")

                    tmp_gffs = []
                    # TODO: need to do better filtering/catching for when sequences are missing in the list?
                    tmp_fas = []
                    for index, row in rows.iterrows():
                        accession = row["accession_tag"]
                        gff_df = fetch_gff(accession)
                        tmp_gff = write_tmp(gff_df, "gff")
                        tmp_gffs.append(tmp_gff)

                    if len(selected_rows) > 1:
                        if len(selected_rows) > 4:
                            for index, row in rows.iterrows():
                                fa_df = fetch_fa(accession)
                                tmp_fa = write_tmp(fa_df, "fa")
                                tmp_fas.append(tmp_fa)

                            output = html.P("Select up to four Starships to compare.")
                        else:
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
                        className="auto-resize-750",
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


@callback(
    Output("pgv-table", "children"),
    Input("url", "href"),
)
def load_ship_table(href):
    query = """
    SELECT a.accession_tag, f.familyName, t.species
    FROM joined_ships j
    LEFT JOIN taxonomy t ON j.taxid = t.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    LEFT JOIN accessions a ON j.ship_id = a.id
    LEFT JOIN ships s on s.accession = a.id
    LEFT JOIN gff g ON a.id = g.ship_id
    WHERE s.ship_sequence is NOT NULL AND g.ship_ID is NOT NULL AND j.orphan IS NULL
    """
    table_df = pd.read_sql_query(query, engine)

    if href:
        table = make_ship_table(df=table_df, id="pgv-table", columns=None, pg_sz=15)
        return table


# inject_svg_to_html(
#     "tmp/genbank_comparison_by_blast.svg",
#     "/home/adrian/anaconda3/lib/python3.8/site-packages/pygenomeviz/viewer/pgv-viewer-template.html",
#     "tmp/genbank_comparison_by_blast.html",
# )

# from Bio.SeqFeature import SeqFeature

# def to_stack_features(features: list[SeqFeature]) -> list[list[SeqFeature]]:
#     """Convert feature list to non-overlap stack feature list of lists

#     Parameters
#     ----------
#     features : list[SeqFeature]
#         Features

#     Returns
#     -------
#     stack_features : list[list[SeqFeature]]
#         Stacked features
#     """
#     sorted_features = sorted(features, key=lambda f: int(f.location.start))  # type: ignore

#     def is_overlap(feature1: SeqFeature, feature2: SeqFeature) -> bool:
#         """Check if features overlap each other"""
#         start1, end1 = int(feature1.location.start), int(feature1.location.end)  # type: ignore
#         start2, end2 = int(feature2.location.start), int(feature2.location.end)  # type: ignore
#         return start1 < end2 and start2 < end1

#     stack_features: list[list[SeqFeature]] = []
#     for feature in sorted_features:
#         placed = False
#         for sublist_features in stack_features:
#             if not is_overlap(feature, sublist_features[-1]):
#                 sublist_features.append(feature)
#                 placed = True
#                 break
#         if not placed:
#             stack_features.append([feature])

#     return stack_features


# def add_stacked_features(track, features):
#     """Add stacked features to a GenomeViz track"""
#     stacked_features = to_stack_features(features)

#     # Iterate over each stack (sublists of non-overlapping features)
#     for idx, stack in enumerate(stacked_features):
#         # Adjust y_offset for each stack, with increasing offset for each new row
#         y_offset = idx * 1  # Adjust this value as needed for spacing between stacks
#         for feature in stack:
#             # Get necessary feature details
#             start = feature.location.start
#             end = feature.location.end
#             strand = feature.location.strand

#             # Set color for the feature based on some condition (you can customize this)
#             if "tyr" in feature.qualifiers.get("Alias", [""])[0]:
#                 color = "tomato"
#             else:
#                 color = "grey"

#             # Add feature to the track with a y_offset for stacking
#             track.add_feature(
#                 start=int(start),
#                 end=int(end),
#                 strand=strand,
#                 plotstyle="bigarrow",
#                 color=color,
#                 label_type="gene",
#                 y_offset=y_offset,  # Ensure stacked features have different y offsets
#             )
