import dash
import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html, callback
from dash.exceptions import PreventUpdate
from dash.dependencies import Output, Input, State
import dash_bio as dashbio

import io
from io import StringIO
import os
import re
import tempfile
import base64
from datetime import date
import subprocess
import json
import pandas as pd

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Gff
from matplotlib.lines import Line2D
from pygenomeviz.utils import ColorCycler
from pygenomeviz.parser import Genbank
from pygenomeviz.align import Blast, AlignCoord
from jinja2 import Template

from src.components.shipTable import make_table
from src.data.joined_ships import df

dash.register_page(__name__)

ColorCycler.set_cmap("tab10")

df_sub = df[
    [
        "starshipID",
        "captain_superfamily",
        "starship_family",
        "starship_navis",
        "starship_haplotype",
        "genus",
        "species",
    ]
]

layout = dbc.Container(
    fluid=True,
    children=[
        dbc.Row(
            justify="center",
            align="top",
            children=[
                dbc.Col(
                    style={"padding": "20px"},
                    sm=12,
                    lg=8,
                    children=[make_table(df_sub)],
                )
            ],
        )
    ],
)


def plot_legend(gv):
    fig = gv.plotfig(fast_render=False)

    # Plot legend for groups
    _ = fig.legend(
        handles=[
            Line2D(
                [], [], marker=">", color="tomato", label="captain", ms=12, ls="none"
            ),
            Line2D(
                [], [], marker=">", color="lime", label="auxillary", ms=12, ls="none"
            ),
            Line2D([], [], marker=">", color="grey", label="Others", ms=12, ls="none"),
        ],
        fontsize=12,
        title="Groups",
        title_fontsize=12,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        handlelength=1.0,
    )
    return fig


def add_gene_feature(gene, track):
    # Get gene name in GFF attributes column (e.g. `gene=araD;`)
    gene_name = str(gene.qualifiers.get("Alias", [""])[0])
    print(gene_name)

    start = gene.location.start
    end = gene.location.end
    strand = gene.strand

    # Set user-defined feature color based on gene name
    if "tyr" in gene_name:
        color = "tomato"
    elif any(substring in gene_name for substring in ["fre", "DUF3723", "nlr", "plp"]):
        color = "lime"
    else:
        color = "grey"

    track.add_feature(start, end, strand, color=color, label_type="gene")


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

    print(f"HTML with embedded SVG saved to {output_html_file}")


def single_pgv(gff_file):
    gff = Gff(gff_file)

    gv = GenomeViz()
    gv.set_scale_xticks()
    # gv.set_scale_bar(ymargin=0.5)

    seq_size = gff.get_seqid2size()
    print(seq_size)

    for seqid, size in gff.get_seqid2size().items():
        track = gv.add_feature_track(seqid, size, labelsize=15)
        track.add_sublabel(size=10, color="grey")
        gene_features = gff.get_seqid2features(feature_type="gene")[seqid]
        for gene in gene_features:
            add_gene_feature(gene, track)
    gv.savefig_html("tmp/gff_features.html")


# single_pgv("tmp/gff/altals1_s00058.gff")


def multi_pgv(gff_files, fna_files):

    gff_list = list(map(Gff, gff_files))

    gv = GenomeViz(track_align_type="center", fig_track_height=0.7)
    gv.set_scale_bar()

    for gff in gff_list:
        track = gv.add_feature_track(gff.name, gff.get_seqid2size(), align_label=False)
        for seqid, features in gff.get_seqid2features("gene").items():
            segment = track.get_segment(seqid)
            for gene in features:
                add_gene_feature(gene, segment)

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
    gv.savefig_html("tmp/genbank_comparison_by_blast.html", fig)
    # gv.savefig("tmp/genbank_comparison_by_blast.svg")


@callback(
    Output("pgv-figure", "figure"),
    [
        Input("table", "derived_virtual_data"),
        Input("table", "derived_virtual_selected_rows"),
    ],
)
def update_pgv(table_data, selected_rows):
    if selected_rows:
        table_df = pd.DataFrame(table_data)
        rows = df.iloc[selected_rows]
        gff_files = rows["gff3"].tolist()
        fna_files = 
        if selected_rows > 1:
            multi_pgv(gff_files,fna_files)


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
#             strand = feature.strand

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
