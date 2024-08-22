import warnings

warnings.filterwarnings("ignore")
import io
import os
import base64
import tempfile
import matplotlib.pyplot as plt
import plotly.express as px

import logomaker as lm

from Bio.Align.Applications import ClustalwCommandline


def agg_df(df, groups):
    agg = df.groupby(groups).starshipID.agg(
        count="count", nunique="nunique", duplicates=lambda x: x.size - x.nunique()
    )

    agg = agg.reset_index()

    return agg


def create_sunburst_plot(df, groups, title):
    selection = agg_df(df, groups)

    pie = px.sunburst(
        selection,
        path=groups,
        values="count",
    )

    pie.update_layout(
        autosize=True,
        title_font=dict(size=24),
        title={
            "text": title,
            "y": 1,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        margin=dict(t=50, l=0, r=0, b=0),
    )
    return pie


def are_all_strings_same_length(strings):
    return len(set(len(s) for s in strings)) == 1


def make_logo(seqs):

    if not seqs:  # If all sequences are empty, return None
        return None

    temp_in_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False)
    temp_out_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False)

    # Write sequences to the temporary input file
    with open(temp_in_file.name, "w") as file:
        for idx, seq in enumerate(seqs):
            if seq != ".":
                header = f">seq{idx + 1}"
                file.write(f"{header}\n{seq}\n")

    if os.path.exists(temp_in_file.name) and os.path.getsize(temp_in_file.name) > 0:

        # Run ClustalW for alignment
        clustalw_cline = ClustalwCommandline(
            "clustalw",
            infile=temp_in_file.name,
            outfile=temp_out_file.name,
            output="FASTA",
        )
        clustalw_cline()

        # Read the aligned sequences
        with open(temp_out_file.name, "r") as f:
            lines = f.readlines()

        aln_seqs = [
            seq.strip().upper()
            for seq in lines
            if not seq.startswith("#") and not seq.startswith(">")
        ]

        # Check if all sequences are the same length
        if are_all_strings_same_length(aln_seqs):
            fig, ax = plt.subplots()
            counts_mat = lm.alignment_to_matrix(
                sequences=aln_seqs, to_type="counts", characters_to_ignore=".-X"
            )
            lm.Logo(counts_mat, ax=ax)

            # Save the figure to a BytesIO object
            buf = io.BytesIO()
            fig.savefig(buf, format="png")
            buf.seek(0)

            # Encode the image to base64
            img_base64 = base64.b64encode(buf.getvalue()).decode("utf-8")
            plt.close(fig)  # Close the figure to free up resources

            # Clean up temporary files
            os.remove(temp_in_file.name)
            os.remove(temp_out_file.name)
            return img_base64
        else:
            # Clean up temporary files
            os.remove(temp_in_file.name)
            os.remove(temp_out_file.name)
            return None
    else:
        return None
