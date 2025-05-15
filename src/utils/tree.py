import warnings

from Bio import Phylo
import plotly.graph_objs as go

import tempfile
import subprocess
import os
import logging

from src.utils.seq_utils import load_fasta_to_dict
from src.database.sql_manager import fetch_captain_tree, fetch_sf_data

warnings.filterwarnings("ignore")

logger = logging.getLogger(__name__)


default_highlight_colors = {
    "Phoenix": "#00cc96",
    "Hephaestus": "#ab63fa",
    "Tardis": "#ff6692",
    "Serenity": "#fecb52",
    "Prometheus": "#636efa",
    "Enterprise": "#ef553b",
    "Galactica": "#19d3f3",
    "Moya": "#ff97ff",
    "Arwing": "#ffa15a",
    "Voyager": "#b6e880",
    "Family-11": "#f7a799",
    "superfam03-3": "#bbbbbb",
    "superfam03-2": "#bbbbbb",
}

default_highlight_families = default_highlight_colors.keys()


def hex_to_rgba(hex_color):
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    a = 0.6
    return f"rgba({r}, {g}, {b}, {a})"


rgb_colors = {
    key: hex_to_rgba(value) for key, value in default_highlight_colors.items()
}


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


def extract_coordinates(line_shapes):
    for shape in line_shapes:
        x_coords = (shape["x0"], shape["x1"])
        y_coords = (shape["y0"], shape["y1"])

    return x_coords, y_coords


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
    """Define a shape of type 'line' for a branch."""
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
    clade,
    x_start,
    line_shapes,
    tip_dots,
    line_width=1,
    x_coords=0,
    y_coords=0,
):
    """Recursively draw the tree branches, down from the given clade."""
    x_curr = x_coords[clade]
    y_curr = y_coords[clade]

    # All lines and tips will be black
    line_color = "rgb(25,25,25)"

    # Add a dot at the tip if it's a terminal node
    if clade.name in metadata["tip"].values:
        tip_dot = dict(
            x=x_curr,
            y=y_curr,
            mode="markers",
            marker=dict(size=5, color=line_color),
            name=clade.name,
        )
        tip_dots.append(tip_dot)

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

        # Vertical line connecting the children
        vertical_line = get_clade_lines(
            orientation="vertical",
            x_curr=x_curr,
            y_bot=y_bot,
            y_top=y_top,
            line_color=line_color,
            line_width=line_width,
        )

        line_shapes.append(vertical_line)

        # Draw descendants
        for child in clade:
            draw_clade(
                metadata,
                child,
                x_curr,
                line_shapes,
                tip_dots,
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


def get_text_label(
    x,
    y,
    text,
    font_size=12,
    font_color="black",
    x_anchor="center",
    y_anchor="middle",
):
    """Define a text annotation"""
    return dict(
        x=x,
        y=y,
        text=text,
        showarrow=False,
        font=dict(size=font_size, color=font_color),
        xanchor=x_anchor,
        yanchor=y_anchor,
    )


def superfam_highlight(
    metadata,
    highlights,
    x_coords=None,
    y_coords=None,
):
    superfam_df = metadata[metadata["familyName"] == highlights]
    rectangle = None
    scatter = None
    text_label = None
    if not superfam_df.empty:
        highlight_names = superfam_df["tip"].tolist()

        color = (
            superfam_df.iloc[0]["color"]
            if "color" in superfam_df.columns and not superfam_df["color"].isna().all()
            else "rgba(25, 25, 25, 0.6)"
        )
        if x_coords is not None and y_coords is not None:
            x_coord_list = [
                x_coords[clade] for clade in x_coords if clade.name in highlight_names
            ]
            y_coord_list = [
                y_coords[clade] for clade in y_coords if clade.name in highlight_names
            ]

            x_start, x_end = min(x_coord_list) - 1, max(x_coord_list)
            y_start, y_end = min(y_coord_list) - 0.5, max(y_coord_list) + 0.5

            rectangle = get_rectangle(
                x_start,
                x_end,
                y_start,
                y_end,
                fill_color=color,
                border_color=color,
                line_width=2,
            )

            scatter = go.Scatter(
                x=[(x_start + x_end) / 2],
                y=[(y_start + y_end) / 2],
                mode="markers",
                marker=dict(size=10, color="rgba(255,255,255,0)"),
                text=[highlights],
                hoverinfo="text",
                name=highlights,
            )

            text_label = get_text_label(
                # x=(x_end + x_start) / 2,
                x=7,
                y=(y_end + y_start) / 2,
                text=highlights,
                font_size=24,
                font_color="black",
                x_anchor="center",
                y_anchor="middle",
            )
    return rectangle, scatter, text_label


def add_tip_labels(
    tip,
    x_coords=0,
    y_coords=0,
):
    text_label = None

    if x_coords is not None and y_coords is not None and tip:
        x_coord_list = [
            x_coords[clade]
            for clade in x_coords
            if clade.name is not None and clade.name in tip
        ]
        y_coord_list = [
            y_coords[clade]
            for clade in y_coords
            if clade.name is not None and clade.name in tip
        ]

        if x_coord_list and y_coord_list:
            # x_start, x_end = min(x_coord_list) - 1, max(x_coord_list)
            y_start, y_end = min(y_coord_list) - 0.5, max(y_coord_list) + 0.5

            text_label = get_text_label(
                x=4,
                y=(y_end + y_start) / 2,
                text=tip,
                font_size=24,
                font_color="black",
                x_anchor="center",
                y_anchor="middle",
            )

    return text_label


def plot_tree(highlight_families=None, tips=None):
    tree_string = fetch_captain_tree()
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as temp_file:
        temp_file.write(tree_string)
        tree_file = temp_file.name

        tree = Phylo.read(tree_file, "newick")

    metadata = fetch_sf_data()

    metadata["color"] = metadata["familyName"].map(rgb_colors)

    # Add debug logging
    logger.debug(f"Metadata columns: {metadata.columns.tolist()}")
    logger.debug(f"Metadata head: \n{metadata.head()}")
    logger.debug(f"Highlight families: {highlight_families}")
    logger.debug(f"Tips: {tips}")

    # graph_title = "Captain Gene Phylogeny"
    graph_title = None

    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    text_labels = []
    tip_dots = []
    centroids = []
    nodes = []

    draw_clade(
        metadata,
        tree.root,
        0,
        line_shapes,
        tip_dots,
        line_width=1,
        x_coords=x_coords,
        y_coords=y_coords,
    )

    # Always show rectangles and tooltips for all families
    highlights = default_highlight_families
    for highlight in highlights:
        rectangle, scatter, text_label = superfam_highlight(
            metadata,
            highlight,
            x_coords=x_coords,
            y_coords=y_coords,
        )
        line_shapes.append(rectangle)
        centroids.append(scatter)
        # Only add text labels if highlight_families is specified
        if highlight_families is not None and (
            highlight_families == "all" or highlight == highlight_families
        ):
            text_labels.append(text_label)

    if tips:
        for tip in tips:
            tips_labels = add_tip_labels(tip=tip, x_coords=x_coords, y_coords=y_coords)
            text_labels.append(tips_labels)

    nodes = centroids

    layout = dict(
        height=1200,
        # width=1000,
        title=graph_title,
        autosize=True,
        showlegend=False,
        xaxis=dict(
            range=[0, 8],
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
        shapes=line_shapes,
    )

    if highlight_families is not None:
        # legend = {"x": 0, "y": 1}
        # Filter out None elements from text_labels
        annotations = [label for label in text_labels if label is not None]
        font = dict(family="Open Sans")
        layout["annotations"] = annotations
        layout["font"] = font

    fig = go.Figure(data=nodes, layout=layout)

    return fig


def run_mafft(query, ref_msa):
    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    # MODEL = "PROTGTR+G+F"
    fasta_dict = load_fasta_to_dict(query)
    tmp_headers = list(fasta_dict.keys())

    mafft_cmd = f"mafft --thread 2 --auto --addfragments {query} --keeplength {ref_msa} | seqkit grep -n -f {tmp_headers} > {tmp_fasta}"
    subprocess.run(
        mafft_cmd,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return tmp_headers, tmp_fasta


def add_to_tree(query_msa, tree, ref_msa, model, tmp_dir):
    if query_msa:
        # Construct the epa-ng command
        epa_cmd = f"epa-ng -T 2 --redo --ref-msa {ref_msa} --tree {tree} --query {query_msa} --model {model} --out-dir {tmp_dir}"

        # Run the command
        result = subprocess.run(
            epa_cmd,
            shell=True,
            stdout=subprocess.DEVNULL,  # Uncomment to suppress stdout
            stderr=subprocess.DEVNULL,  # Uncomment to suppress stderr
        )

        # Check for errors in command execution
        if result.returncode != 0:
            raise RuntimeError(
                f"epa-ng command failed with return code {result.returncode}"
            )

        tmp_tree = os.path.join(tmp_dir, "epa_result.jplace")
        return tmp_tree


def gappa(tmp_tree):
    out_dir = os.path.dirname(tmp_tree)
    out_tree = os.path.join(out_dir, "epa_result.newick")
    gappa_cmd = (
        f"gappa examine graft --out-dir {out_dir} --threads 2 --jplace-path {tmp_tree}"
    )
    subprocess.run(
        gappa_cmd,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    logging.debug(f"Tree output: {out_tree}")

    return out_tree
