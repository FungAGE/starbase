import warnings

warnings.filterwarnings("ignore")
import io
import os
import base64
import tempfile
import matplotlib.pyplot as plt
import plotly.express as px
import pickle

from Bio.Seq import Seq
import logomaker as lm
from Bio.Align.Applications import ClustalwCommandline

from src.utils.seq_utils import clean_sequence


def agg_df(df, groups):
    if all(group in df.columns for group in groups):
        agg = df.groupby(groups).accession_tag.agg(
            count="count", nunique="nunique", duplicates=lambda x: x.size - x.nunique()
        )

        agg = agg.reset_index()

        return agg
    else:
        raise ValueError(
            f"One or more of the groups {groups} do not exist in the DataFrame"
        )


def create_sunburst_plot(df, type, title_switch=True, cache_bust=None):
    # Define color schemes and settings based on type
    settings = {
        "ship": {
            "groups": ["familyName"],
            "title": "Starships by Family/Navis",
            "color_sequence": px.colors.qualitative.Set3,  # Discrete colors for families
            "hover_data": ["count", "nunique", "duplicates"]
        },
        "tax": {
            "groups": ["order", "family"],
            "title": "Starships by Order/Family",
            "color_sequence": px.colors.qualitative.Pastel,  # Discrete colors for taxonomy
            "hover_data": ["count", "nunique"]
        }
    }
    
    if type not in settings:
        raise ValueError(f"Unknown plot type: {type}")
        
    config = settings[type]
    selection = agg_df(df, config["groups"])
    
    # Create enhanced sunburst plot
    fig = px.sunburst(
        selection,
        path=config["groups"],
        values="count",
        color=config["groups"][0],  # Color by the first grouping level
        color_discrete_sequence=config["color_sequence"],
        custom_data=config["hover_data"],
        branchvalues="total",
        maxdepth=2,
    )
    
    # Enhanced styling
    fig.update_layout(
        template="plotly_white",
        font_family="Arial, sans-serif",
        autosize=True,
        showlegend=True,
        title={
            "text": config["title"] if title_switch else None,
            "y": 0.95,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
            "font": {"size": 24}
        },
        margin=dict(
            t=30,
            l=10,
            r=10,
            b=10,
            pad=4
        ),
        transition={
            "duration": 500,
            "easing": "cubic-in-out"
        },
        hoverlabel=dict(
            bgcolor="white",
            font_size=14,
            font_family="Arial, sans-serif"
        ),
        # Add uirevision to control caching behavior
        uirevision=str(cache_bust) if cache_bust is not None else True
    )
    
    # Customize hover information
    hover_template = (
        "<b>%{label}</b><br><br>" +
        "Count: %{value}<br>" +  # Changed from customdata to value
        "Percentage: %{percentParent:.1%}<br>" +
        "<extra></extra>"
    )
    
    fig.update_traces(
        hovertemplate=hover_template,
        textinfo="label+percent parent",
        insidetextorientation="radial",
        selector=dict(type="sunburst"),
        marker=dict(
            line=dict(color="white", width=1)
        ),
    )
    
    return fig


def are_all_strings_same_length(strings):
    return len(set(len(s) for s in strings)) == 1


def make_logo(seqs, fig_name=None, type=None):

    if not seqs:  # If all sequences are empty, return None
        return None

    temp_in_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False)
    temp_out_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False)

    # Write sequences to the temporary input file
    with open(temp_in_file.name, "w") as file:
        clean_seqs = [str(Seq(seq)) for seq in seqs if seq and clean_sequence(seq)]

        for idx, cs in enumerate(clean_seqs):
            if cs != ".":
                header = f">seq{idx + 1}"
                file.write(f"{header}\n{cs}\n")

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
            seq.strip()
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

            if fig_name:
                fig.savefig(fig_name, format="png")
                return None
            else:
                # Save the figure to a BytesIO object
                buf = io.BytesIO()
                fig.savefig(buf, format="png")
                buf.seek(0)

                img_base64 = base64.b64encode(buf.getvalue()).decode("utf-8")
                output = f"data:image/png;base64,{img_base64}"
                plt.close(fig)

                os.remove(temp_in_file.name)
                os.remove(temp_out_file.name)
                return output

        else:
            os.remove(temp_in_file.name)
            os.remove(temp_out_file.name)
            return None
    else:
        return None
