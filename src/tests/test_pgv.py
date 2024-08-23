import os
import tempfile
import pandas as pd

from dash import html
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Gff
from matplotlib.lines import Line2D
from pygenomeviz.utils import ColorCycler
from pygenomeviz.align import Blast, AlignCoord
from jinja2 import Template

ColorCycler.set_cmap("tab10")

df = pd.read_csv("src/data/joined_ships.csv")


def is_valid_file(file_path):
    if isinstance(file_path, str) and os.path.isfile(file_path):
        return os.path.getsize(file_path) > 0
    return False


df["valid_gff3"] = df["gff3"].apply(lambda x: x if is_valid_file(x) else None)
df["valid_fna"] = df["fna"].apply(lambda x: x if is_valid_file(x) else None)

filtered_df = df.dropna(subset=["valid_gff3"])
filtered_df = filtered_df.drop(columns=["valid_gff3", "valid_fna"])


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
    attributes = gene.qualifiers
    gene_name = str(gene.qualifiers.get("Alias", [""])[0])

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
    align_coords = AlignCoord.filter(align_coords, length_thr=500, identity_thr=30)

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


def test_pgv(df, ships):
    tmp_pgv = tempfile.NamedTemporaryFile(suffix=".html", delete=False).name
    rows = df[df["starshipID"].isin(ships)]

    # Ensure the columns exist in the DataFrame
    if "gff3" in rows.columns and "fna" in rows.columns:
        gffs = rows["gff3"].dropna().tolist()
        fna = rows["fna"].dropna().tolist()

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
    else:
        raise ValueError("Required columns are missing from the data.")


test_pgv(df, ["altals1_s00058"])
