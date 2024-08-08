import dash
from dash import dash_table, dcc, html, callback
from dash.dependencies import Output, Input, State
import dash_bootstrap_components as dbc

from Bio import Phylo
import pandas as pd
import plotly.graph_objs as go
import numpy as np
import re

colors = {
    "superfam01-1": "#8dd3c7",
    "superfam01-2": "#ededa8",
    "superfam01-3": "#adabc4",
    "superfam01-4": "#33a02c",
    "superfam01-5": "#fb8072",
    "superfam02-1": "#80b1d3",
    "superfam02-2": "#b3de69",
    "superfam02-3": "#b0b0b0",
    "superfam03-1": "#fdb45a",
    "superfam03-2": "#fccde5",
    "superfam03-3": "#ffed6f",
    "superfam03-4": "#bc80bd",
    "superfam03-5": "#ccebc5",
}


def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f"rgb({r}, {g}, {b})"


rgb_colors = {key: hex_to_rgb(value) for key, value in colors.items()}


def is_rgb_string(s):
    rgb_pattern = r"^rgb\(\s*(0|[1-9][0-9]?|[1-9][0-9][0-9]?|255)\s*,\s*(0|[1-9][0-9]?|[1-9][0-9][0-9]?|255)\s*,\s*(0|[1-9][0-9]?|[1-9][0-9][0-9]?|255)\s*\)$"
    return bool(re.match(rgb_pattern, s))


def get_x_coordinates(tree):
    """Associates to  each clade an x-coord.
    returns dict {clade: x-coord}
    """
    xcoords = tree.depths()
    # tree.depth() maps tree clades to depths (by branch length).
    # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth
    # is the distance from root to clade

    #  If there are no branch lengths, assign unit branch lengths
    if not max(xcoords.values()):
        xcoords = tree.depths(unit_branch_lengths=True)
    return xcoords


def get_y_coordinates(tree, dist=1):
    """
    returns  dict {clade: y-coord}
    The y-coordinates are  (float) multiple of integers (i*dist below)
    dist depends on the number of tree leafs
    """
    maxheight = tree.count_terminals()  # Counts the number of tree leafs.
    # Rows are defined by the tips/leafs
    ycoords = dict(
        (leaf, maxheight - i * dist)
        for i, leaf in enumerate(reversed(tree.get_terminals()))
    )

    def calc_row(clade):
        for subclade in clade:
            if subclade not in ycoords:
                calc_row(subclade)
        ycoords[clade] = (ycoords[clade.clades[0]] + ycoords[clade.clades[-1]]) / 2

    if tree.root.clades:
        calc_row(tree.root)
    return ycoords


def get_clade_lines(
    orientation="horizontal",
    y_curr=0,
    x_start=0,
    x_curr=0,
    y_bot=0,
    y_top=0,
    line_color=None,
    line_width=0.5,
):
    """define a shape of type 'line', for branch"""
    branch_line = dict(
        type="line", layer="below", line=dict(color=line_color, width=line_width)
    )
    if orientation == "horizontal":
        branch_line.update(x0=x_start, y0=y_curr, x1=x_curr, y1=y_curr)
    elif orientation == "vertical":
        branch_line.update(x0=x_curr, y0=y_bot, x1=x_curr, y1=y_top)
    else:
        raise ValueError("Line type can be 'horizontal' or 'vertical'")

    return branch_line


def draw_clade(
    metadata,
    rgb_colors,
    clade,
    x_start,
    line_shapes,
    line_width=1,
    x_coords=0,
    y_coords=0,
):
    """Recursively draw the tree branches, down from the given clade"""
    x_curr = x_coords[clade]
    y_curr = y_coords[clade]

    if clade in metadata["tip"].tolist():
        df = metadata[metadata["tip"] == clade]
        line_color = df.loc[0, "color"]
    else:
        line_color = "rgb(25,25,25)"

    # Draw a horizontal line from start to here
    branch_line = get_clade_lines(
        orientation="horizontal",
        y_curr=y_curr,
        x_start=x_start,
        x_curr=x_curr,
        line_color=line_color,
        line_width=line_width,
    )

    line_shapes.append(branch_line)

    if clade.clades:
        # Draw a vertical line connecting all children
        y_top = y_coords[clade.clades[0]]
        y_bot = y_coords[clade.clades[-1]]

        line_shapes.append(
            get_clade_lines(
                orientation="vertical",
                x_curr=x_curr,
                y_bot=y_bot,
                y_top=y_top,
                line_color=line_color,
                line_width=line_width,
            )
        )

        # Draw descendants
        for child in clade:
            draw_clade(
                metadata,
                rgb_colors,
                child,
                x_curr,
                line_shapes,
                x_coords=x_coords,
                y_coords=y_coords,
            )


def get_rectangle(
    x_start,
    x_end,
    y_start,
    y_end,
    fill_color="rgba(0,0,0,0.1)",  # Semi-transparent color
    border_color="rgba(0,0,0,0.5)",
    line_width=1,
):
    """Define a rectangle shape"""
    return dict(
        type="rect",
        x0=x_start,
        x1=x_end,
        y0=y_start,
        y1=y_end,
        fillcolor=fill_color,
        line=dict(color=border_color, width=line_width),
    )


def create_tree(tree_file, metadata, highlight_names=None):
    tree = Phylo.read(tree_file, "newick")
    graph_title = "Captain Gene Phylogeny"

    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    draw_clade(
        metadata,
        rgb_colors,
        tree.root,
        0,
        line_shapes,
        line_width=1,
        x_coords=x_coords,
        y_coords=y_coords,
    )

    # Draw rectangles for the specified names
    if highlight_names:
        x_start = (
            min(x_coords[name] for name in highlight_names if name in x_coords) - 0.5
        )  # Adjust as needed
        x_end = max(x_coords.items[highlight_names]) + 0.5  # Adjust as needed
        y_start = min(y_coords.items[highlight_names]) - 0.5  # Adjust as needed
        y_end = max(y_coords.items[highlight_names]) + 0.5  # Adjust as needed
        rectangle = get_rectangle(
            x_start=x_start,
            x_end=x_end,
            y_start=y_start,
            y_end=y_end,
            fill_color="rgba(255,0,0,0.2)",  # Highlight color
            border_color="rgba(255,0,0,0.8)",
            line_width=2,
        )
        line_shapes.append(rectangle)

    my_tree_clades = x_coords.keys()
    X = []
    Y = []
    text = []

    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
        text.append(cl.name)

    layout = dict(
        height=1200,
        width=1000,
        title=graph_title,
        autosize=True,
        automargin=True,
        showlegend=True,
        xaxis=dict(
            range=[0, 6],
            showline=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
        ),
        yaxis=dict(
            range=[0, 1250],
            showline=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
        ),
        hovermode="closest",
        shapes=line_shapes,
        legend={"x": 0, "y": 1},
        font=dict(family="Open Sans"),
    )

    nodes = []

    fig = dict(data=nodes, layout=layout)
    return fig


def plot_tree():
    tree_file = "src/data/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.treefile"
    metadata = pd.read_csv("src/data/superfam-clades.tsv", sep="\t")
    metadata["color"] = metadata["superfam"].map(rgb_colors)

    # Specify names to highlight
    highlight_names = metadata[metadata["superfam"] == "superfam01-1"]["tip"].tolist()

    fig = create_tree(tree_file, metadata, highlight_names=highlight_names)

    layout = html.Div(
        [dcc.Graph(id="phylogeny-graph", className="div-card", figure=fig)]
    )
    return layout


# @callback(
#     Output("phylogeny-graph", "figure"),
#     [
#         Input("table", "derived_virtual_data"),
#         Input("table", "derived_virtual_selected_rows"),

#     ],
# )
# def update_phylogeny_tree(data,selected_row):

#     if tree_file_filtred in tree_fig:
#         fig = tree_fig[tree_file_filtred]
#     else:
#         if ord_by_elt == "Country" or ord_by_elt == "Division" or ord_by_elt == "Date":
#             fig = create_tree(
#                 virus_name, tree_file_filtred, metadata_file_filtred, ord_by_elt
#             )

#         tree_fig[tree_file_filtred] = fig

#     return fig
