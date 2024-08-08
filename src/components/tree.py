import dash
from dash import dash_table, dcc, html, callback
from dash.dependencies import Output, Input, State
import dash_bootstrap_components as dbc

from Bio import Phylo
import pandas as pd
import plotly.graph_objs as go
import numpy as np


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
    line_color="rgb(25,25,25)",
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
    clade,
    x_start,
    line_shapes,
    line_color="rgb(15,15,15)",
    line_width=1,
    x_coords=0,
    y_coords=0,
):
    """Recursively draw the tree branches, down from the given clade"""

    x_curr = x_coords[clade]
    y_curr = y_coords[clade]

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
            draw_clade(child, x_curr, line_shapes, x_coords=x_coords, y_coords=y_coords)


def create_tree(tree_file, metadata, superfamily):
    tree = Phylo.read(tree_file, "newick")
    graph_title = "Captain Gene Phylogeny"

    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    draw_clade(
        tree.root,
        0,
        line_shapes,
        line_color="rgb(25,25,25)",
        line_width=1,
        x_coords=x_coords,
        y_coords=y_coords,
    )
    my_tree_clades = x_coords.keys()
    X = []
    Y = []
    text = []

    axis = dict(
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        title="",  # y title
    )

    label_legend = set(list(metadata[superfamily]))
    nodes = []

    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
        text.append(cl.name)

    layout = dict(
        height=1200,
        width=800,
        title=graph_title,
        # dragmode="select",
        autosize=True,
        automargin=True,
        showlegend=True,
        xaxis=dict(
            showline=True,
            zeroline=False,
            showgrid=True,  # To visualize the vertical lines
            ticklen=4,
            showticklabels=True,
            title="Branch Length",
        ),
        yaxis=axis,
        hovermode="closest",
        shapes=line_shapes,
        legend={"x": 0, "y": 1},
        font=dict(family="Open Sans"),
    )

    fig = dict(data=nodes, layout=layout)
    return fig


def plot_tree():
    tree_file = "src/data/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.treefile"
    metadata = pd.read_csv("src/data/joined_ships.csv")

    fig = create_tree(tree_file, metadata, "superfamily")
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
