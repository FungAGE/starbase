from dash import html
from ete3 import Tree, TreeStyle, NodeStyle

tree_file = "src/data/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.treefile.headers"
nodes = "CryptonF1RhMu"


def generate_tree(tree_file, nodes):
    colors = {
        "1834": {"superfam": "fam1-1", "color": "#8dd3c7"},
        "1811": {"superfam": "fam1-2", "color": "#ededa8"},
        "2032": {"superfam": "fam1-3", "color": "#adabc4"},
        "2087": {"superfam": "fam1-4", "color": "#33a02c"},
        "2153": {"superfam": "fam1-5", "color": "#fb8072"},
        "1235": {"superfam": "fam2-1", "color": "#80b1d3"},
        "1384": {"superfam": "fam2-2", "color": "#b3de69"},
        "1497": {"superfam": "fam2-3", "color": "#b0b0b0"},
        "1643": {"superfam": "fam3-1", "color": "#fdb45a"},
        "1748": {"superfam": "fam3-2", "color": "#fccde5"},
        "1778": {"superfam": "fam3-3", "color": "#ffed6f"},
        "1531": {"superfam": "fam3-4", "color": "#bc80bd"},
        "1606": {"superfam": "fam3-5", "color": "#ccebc5"},
    }

    t = Tree(tree_file)

    ancestor = t.get_common_ancestor(nodes)
    tip = t.search_nodes(name=nodes)

    # Let's now add some custom features to our nodes. add_features can be
    # used to add many features at the same time.
    tip.add_features(vowel=False, confidence=1.0)
    ancestor.add_features(color=colors.get(leaf.name, "none"))

    # Define a node style
    nstyle = NodeStyle()
    nstyle["fgcolor"] = "blue"
    nstyle["size"] = 10

    # Apply the node style to all nodes
    for n in t.traverse():
        n.set_style(nstyle)

    # Create a tree style
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True

    # Save the tree as an HTML file
    t.render("src/data/tree.html", tree_style=ts)

    return html.Iframe(
        id="tree",
        srcDoc=open("src/data/tree.png", "r").read(),
        width="100%",
        height="600px",
    )
