import dash
import dash_bootstrap_components as dbc
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

from src.components.tables import make_ship_table

dash.register_page(__name__)

ColorCycler.set_cmap("tab10")

layout = dbc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        dbc.Row(
            justify="center",
            align="top",
            children=[
                dbc.Col(
                    style={"padding": "20px"},
                    sm=12,
                    lg=8,
                    children=[
                        html.Div(id="pgv-table"),
                        dbc.Button(
                            "Show Starship(s) in Viewer",
                            id="update-button",
                            n_clicks=0,
                            className="d-grid gap-2 col-4 mx-auto",
                            style={"fontSize": "1rem"},
                        ),
                    ],
                )
            ],
        ),
        dbc.Row(
            justify="center",
            align="top",
            children=[
                dbc.Col(
                    style={"padding": "20px"},
                    sm=12,
                    lg=8,
                    children=[
                        dcc.Loading(
                            id="loading-1",
                            type="default",
                            children=[
                                html.Div(id="pgv-figure"),
                                html.Div(
                                    id="pgv-figure-error-message",
                                    style={"color": "red"},
                                ),
                            ],
                        )
                    ],
                )
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


def multi_pgv(gff_files, fna_files, tmp_file):

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
    align_coords = Blast(fna_files, seqtype="protein").run()
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
        State("pgv-table", "selected_rows"),
        State("pgv-table", "data"),
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

                    if "gff3" in rows.columns and "fna" in rows.columns:
                        gffs = rows["gff3"].dropna().tolist()
                        fna = rows["fna"].dropna().tolist()

                        if len(selected_rows) > 1:
                            if len(selected_rows) > 4:
                                return html.Div(
                                    "Select up to four Starships to compare."
                                )
                            else:
                                multi_pgv(gffs, fna, tmp_pgv)
                        elif len(selected_rows) == 1:
                            gff_file = rows.iloc[0]["gff3"]
                            single_pgv(gff_file, tmp_pgv)
                        else:
                            return html.Div("No valid selection.")

                        try:
                            with open(tmp_pgv, "r") as file:
                                pgv_content = file.read()
                        except IOError:
                            return html.Div("Failed to read the temporary file.")

                        return html.Div(
                            [
                                html.Iframe(
                                    srcDoc=pgv_content,
                                    style={
                                        "width": "100%",
                                        "height": "500px",
                                        "border": "none",
                                    },
                                )
                            ]
                        )
                    else:
                        return html.H4("Required columns are missing from the data.")
                else:
                    return html.H4("Invalid row selection.")
            except Exception as e:
                print("Exception:", e)
                return html.H4("Error processing data.")
        else:
            return html.H4("Select Starship(s) to visualize.")
    else:
        return html.H4(
            "Select the Starship(s) in the table above and click the button to visualize."
        )


def is_valid_file(file_path):
    if isinstance(file_path, str) and os.path.isfile(file_path):
        return os.path.getsize(file_path) > 0
    return False


@callback(
    Output("pgv-table", "children"),
    [
        Input("joined-ships", "data"),
        Input("url", "href"),
    ],
)
def load_ship_table(cached_data, href):
    if href:
        initial_df = pd.DataFrame(cached_data)

        specified_columns = [
            "starshipID",
            "starship_family",
            "genus",
            "species",
        ]

        initial_df["valid_gff3"] = initial_df["gff3"].apply(
            lambda x: x if is_valid_file(x) else None
        )
        initial_df["valid_fna"] = initial_df["fna"].apply(
            lambda x: x if is_valid_file(x) else None
        )

        filtered_df = initial_df.dropna(subset=["valid_gff3"])
        filtered_df = filtered_df.drop(columns=["valid_gff3", "valid_fna"])
        table = make_ship_table(filtered_df, "pgv-table", specified_columns)
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
