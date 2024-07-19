from pygenomeviz import GenomeViz
from pygenomeviz.parser import Gff
from pygenomeviz.utils import ColorCycler

ColorCycler.set_cmap("tab10")

gff = Gff("tmp/gff/altals1_s00058.gff")

gv = GenomeViz(fig_track_height=0.8, feature_track_ratio=0.4)
gv.set_scale_xticks()
# gv.set_scale_bar(ymargin=0.5)

# color = ColorCycler()

track = gv.add_feature_track(gff.name, gff.get_seqid2size())

gene_features = gff.extract_features("gene")

for feature in gene_features:
    # Get gene name in GFF attributes column (e.g. `gene=araD;`)
    gene_name = str(feature.qualifiers.get("Alias", [""])[0])
    print(gene_name)

    # Set user-defined feature color based on gene name
    if "tyr" in gene_name:
        color = "tomato"
    elif any(substring in gene_name for substring in ["fre", "DUF3723", "nlr", "plp"]):
        color = "lime"
    else:
        color = "grey"

    track.add_features(feature, plotstyle="bigarrow", color=color, label_type="gene")

# fig = gv.plotfig()

# # Plot legend for groups
# _ = fig.legend(
#     handles=[
#         Line2D([], [], marker=">", color="tomato", label="captain", ms=12, ls="none"),
#         Line2D([], [], marker=">", color="lime", label="auxillary", ms=12, ls="none"),
#         Line2D([], [], marker=">", color="grey", label="Others", ms=12, ls="none"),
#     ],
#     fontsize=12,
#     title="Groups",
#     title_fontsize=12,
#     loc="center left",
#     bbox_to_anchor=(1.02, 0.5),
#     handlelength=1.0,
# )
# fig.savefig("result.png")

gv.savefig_html("gff_features.html")
