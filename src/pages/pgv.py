import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import dcc, html, callback
from dash.dependencies import Output, Input, State

import os
import tempfile
import pandas as pd
import logging

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Gff
from matplotlib.lines import Line2D
from pygenomeviz.utils import ColorCycler
from pygenomeviz.align import Blast, AlignCoord, MMseqs, MUMmer
from Bio import SeqIO
from jinja2 import Template

from src.components.cache import cache
from src.components.mariadb import engine
from src.components.tables import make_ship_table
from src.components.cache_manager import load_from_cache
from src.components.sql_queries import (
    fetch_all_ships,
    fetch_accession_gff,
    fetch_ship_table,
)


logger = logging.getLogger(__name__)

dash.register_page(__name__)

ColorCycler.set_cmap("tab10")

table_columns = [
    {
        "name": "Accession",
        "id": "accession_tag",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
    {
        "name": "Starship Family",
        "id": "familyName",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
    {
        "name": "Species",
        "id": "species",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
]

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
                            ],
                        ),
                    ],
                ),
                dmc.GridCol(
                    span={
                        "sm": 12,
                        "lg": 8,
                    },
                    children=[
                        dmc.Center(
                            [
                                dmc.Button(
                                    "Show Starship(s) in Viewer",
                                    id="update-button",
                                    n_clicks=0,
                                    className="text-custom text-custom-sm text-custom-md",
                                ),
                            ]
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
                                html.Div(id="pgv-message"),
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


def write_tmp(df, seqid, file_type=None, temp_dir=None):
    # Create the file path with your desired naming convention
    file_path = os.path.join(temp_dir, f"{seqid}.{file_type}")

    # Write to the file based on the type
    if file_type == "gff":
        df.to_csv(file_path, sep="\t", header=False, index=False)
    elif file_type == "fa":
        seq = df["sequence"][0]
        with open(file_path, "w") as f:
            f.write(f">{seqid}\n{seq}\n")

    return file_path


def load_gff(accession):

    df = load_from_cache(f"accession_gff_{accession}")
    if df is None:
        df = fetch_accession_gff(accession)

    if df.empty:
        print(f"No GFF records found for accession_tag: {accession}")
    else:
        # Identify the rows where 'end' is less than 'start'
        mask = df["end"] < df["start"]

        # Swap the 'start' and 'end' values where the mask is True
        df.loc[mask, ["start", "end"]] = df.loc[mask, ["end", "start"]].values

    return df


def load_fa(accession):
    df = load_from_cache("all_ships")
    if df is None:
        df = fetch_all_ships()

    if isinstance(accession, str):
        df = df[df["accession_tag"] == accession]
    else:
        df = df[df["accession_tag"].isin(accession)]

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


def is_valid_sequence_file(file_path):
    valid_extensions = (".fa", ".fna", ".fasta", ".gb", ".gbk", ".gbff")
    return file_path.endswith(valid_extensions)


def multi_pgv(gff_files, seqs, tmp_file):
    gff_list = list(map(Gff, gff_files))
    gv = GenomeViz(track_align_type="center", fig_track_height=0.7)
    gv.set_scale_bar()

    for gff in gff_list:
        # Check GFF SeqIDs
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
                return
            if os.path.getsize(seq_file) == 0:
                print(f"Error: {seq_file} is an empty file.")
                return
            if not is_valid_sequence_file(seq_file):
                print(f"Error: {seq_file} does not have a valid extension.")
                return
    # BLAST
    len_thr = 50
    id_thr = 30
    align_coords = Blast(seqs, seqtype="nucleotide").run()

    align_coords = AlignCoord.filter(
        align_coords, length_thr=len_thr, identity_thr=id_thr
    )

    # Run MMseqs RBH search
    # align_coords = MMseqs(seqs, threads=2).run()

    # Run MuMMer
    # align_coords = MUMmer(seqs).run()

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
        message = None
    else:
        message = html.H4(
            f"No alignments found between ships (Length threshold: {len_thr}bp, ID threshold: {id_thr}%)"
        )

    fig = plot_legend(gv)
    gv.savefig_html(tmp_file, fig)

    return message


@cache.memoize()
@callback(
    Output("pgv-table", "children"),
    Input("url", "href"),
)
def load_ship_table(href):
    table_df = load_from_cache("ship_table")
    if table_df is None:
        table_df = fetch_ship_table()
    if href:
        table = make_ship_table(
            df=table_df, id="pgv-table", columns=table_columns, pg_sz=15
        )
        return table


@callback(
    [Output("pgv-figure", "children"), Output("pgv-message", "children")],
    Input("update-button", "n_clicks"),
    [
        State("pgv-table", "derived_virtual_selected_rows"),
        State("pgv-table", "derived_virtual_data"),
    ],
)
def update_pgv(n_clicks, selected_rows, table_data):
    message = None
    if n_clicks > 0:
        tmp_pgv = tempfile.NamedTemporaryFile(suffix=".html", delete=True).name

        if table_data and selected_rows is not None:
            try:
                table_df = pd.DataFrame(table_data)
                # Debug: Check columns and selected_rows
                logger.info(f"Columns in table_df: {table_df.columns.tolist()}")
                logger.info(f"Selected rows: {selected_rows}")

                if isinstance(selected_rows, list) and all(
                    isinstance(idx, int) for idx in selected_rows
                ):
                    try:
                        rows = table_df.iloc[selected_rows]
                    except IndexError as e:
                        return html.Div("Index out of bounds")

                    with tempfile.TemporaryDirectory() as temp_dir:
                        tmp_gffs = []
                        tmp_fas = []
                        for index, row in rows.iterrows():
                            accession = row["accession_tag"]
                            logging.info(f"Fetching GFF for accession: {accession}")
                            gff_df = load_gff(accession)

                            tmp_gff = write_tmp(gff_df, accession, "gff", temp_dir)
                            tmp_gffs.append(tmp_gff)
                            logging.info(f"Fetching FA for accession: {accession}")
                            fa_df = load_fa(accession)
                            tmp_fa = write_tmp(fa_df, accession, "fa", temp_dir)
                            tmp_fas.append(str(tmp_fa))
                            output = html.P("Select up to four Starships to compare.")
                        if len(selected_rows) > 1 and len(selected_rows) <= 4:
                            message = multi_pgv(tmp_gffs, tmp_fas, tmp_pgv)
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
                logging.error(f"Exception: {e}")
                output = html.H4(
                    "Error while comparing ships using BLAST. Try another combination."
                )
        else:
            output = html.H4("Select Starship(s) to visualize.")
    else:
        output = html.H4(
            "Select up to 4 Starships in the table above and click the button to visualize."
        )
    return (
        html.Div(
            [output],
            className="center-content text-center",
        ),
        message,
    )


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
