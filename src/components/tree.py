from Bio import Phylo
import pandas as pd
import plotly.graph_objs as go

tree_file = "src/data/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.treefile"
metadata = pd.read_csv("src/data/superfam-clades.tsv", sep="\t")

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

default_highlight_clades = [
    "superfam01-1",
    "superfam01-2",
    "superfam01-3",
    "superfam01-4",
    "superfam01-5",
    "superfam02-1",
    "superfam02-2",
    "superfam02-3",
    "superfam03-1",
    "superfam03-2",
    "superfam03-3",
    "superfam03-4",
    "superfam03-5",
]


def hex_to_rgba(hex_color):
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    a = 0.6
    return f"rgba({r}, {g}, {b}, {a})"


rgb_colors = {key: hex_to_rgba(value) for key, value in colors.items()}

metadata["color"] = metadata["superfam"].map(rgb_colors)


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

    if clade.name in metadata["tip"].values:
        line_color = metadata.loc[metadata["tip"] == clade.name, "color"].values[0]
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
    superfam_clade,
    x_coords=None,
    y_coords=None,
):
    df = metadata[metadata["superfam"] == superfam_clade]

    if not df.empty:
        highlight_names = df["tip"].tolist()

        color = (
            df.iloc[0]["color"]
            if "color" in df.columns and not df["color"].isna().all()
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
            text=[superfam_clade],
            hoverinfo="text",
            name=superfam_clade,
        )

        text_label = get_text_label(
            # x=(x_end + x_start) / 2,
            x=7,
            y=(y_end + y_start) / 2,
            text=superfam_clade,
            font_size=24,
            font_color="black",
            x_anchor="center",
            y_anchor="middle",
        )
        return rectangle, scatter, text_label


def plot_tree(tree_file, metadata, highlight_clades=default_highlight_clades):
    tree = Phylo.read(tree_file, "newick")

    graph_title = "Captain Gene Phylogeny"

    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    text_labels = []
    scatter_points = []
    nodes = []

    draw_clade(
        metadata,
        tree.root,
        0,
        line_shapes,
        line_width=1,
        x_coords=x_coords,
        y_coords=y_coords,
    )

    if highlight_clades is not None:
        for superfam_clade in highlight_clades:
            rectangle, scatter, text_label = superfam_highlight(
                metadata,
                superfam_clade,
                x_coords=x_coords,
                y_coords=y_coords,
            )
            line_shapes.append(rectangle)
            scatter_points.append(scatter)
            text_labels.append(text_label)
        nodes = scatter_points

    layout = dict(
        height=1200,
        width=1000,
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

    if highlight_clades is not None:
        legend = {"x": 0, "y": 1}
        annotations = text_labels
        font = dict(family="Open Sans")
        # layout["legend"] = legend
        layout["annotations"] = annotations
        layout["font"] = font

    fig = go.Figure(data=nodes, layout=layout)

    return fig
