import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import dcc, html, callback, no_update

from dash.dependencies import Output, Input, State

import os
import tempfile
import pandas as pd
import logging

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Gff
from matplotlib.lines import Line2D
from pygenomeviz.utils import ColorCycler
from pygenomeviz.align import Blast, AlignCoord, MMseqs, MUMmer

from src.components.callbacks import create_modal_callback
from src.components.error_boundary import handle_callback_error, create_error_boundary



logger = logging.getLogger(__name__)

dash.register_page(__name__)

ColorCycler.set_cmap("tab10")

table_columns = [
    {
        "id": "accession_tag",
        "name": "Accession",
        "selectable": True,
    },
    {
        "id": "familyName",
        "name": "Starship Family",
        "selectable": True,
    },
    {
        "id": "species",
        "name": "Species",
        "selectable": True,
    },
]

modal = dmc.Modal(
    id="pgv-modal",
    opened=False,
    centered=True,
    overlayProps={"blur": 3},
    size="lg",
    children=[
        dmc.Title(id="pgv-modal-title", order=3),
        dmc.Space(h="md"),
        html.Div(id="pgv-modal-content"),
    ],
)

layout = create_error_boundary(
    dmc.Container(
        fluid=True,
        children=[
            dcc.Location(id="url", refresh=False),
        
        # Header Section
        dmc.Paper(
            children=[
                dmc.Title("Starship Genome Viewer", order=1, mb="md"),
                dmc.Text(
                    "Compare and visualize up to 4 Starship sequences",
                    size="lg",
                    c="dimmed",
                ),
            ],
            p="xl",
            radius="md",
            withBorder=False,
            mb="xl",
        ),
        
        # Main content
        dmc.Grid(
            children=[
                dmc.GridCol(
                    span={"base": 12, "md": 6},
                    children=dmc.Paper(
                                    children=[
                                        dmc.Title("Select Starships", order=2),
                                        dmc.Button(
                                            dmc.Text("Show Selected Starship(s)",size="lg"),
                                            id="update-button",
                                            variant="gradient",
                                            gradient={"from": "indigo", "to": "cyan"},
                                            leftSection=html.I(className="bi bi-eye"),
                                        ),
                                    ],
                                    p="xl",
                                    radius="md",
                                    withBorder=True,
                                    style={"position": "relative"}
                        ),
                        id="pgv-controls-container"
                ),
                
                # Right Column - Visualization Section
                dmc.GridCol(
                    span={"base": 12, "md": 6},  # Full width on mobile, half on medium+ screens
                    children=[
                        dmc.Paper(
                            children=dmc.Stack([
                                # Message Area
                                html.Div(
                                    id="pgv-message",
                                    style={"textAlign": "center"},
                                ),
                                # Visualization Area
                                dcc.Loading(
                                    id="loading-1",
                                    type="circle",
                                    children=html.Div(
                                        id="pgv-figure",
                                        style={
                                            "height": "800px",
                                            "width": "100%",
                                            "overflow": "auto",
                                            "backgroundColor": "#f8f9fa",
                                            "border": "1px solid #dee2e6",
                                            "borderRadius": "4px",
                                        },
                                    ),
                                ),
                            ], gap="md"),
                            p="xl",
                            radius="md",
                            withBorder=True,
                            h="100%",  # Make paper fill the height
                        ),
                    ],
                ),
            ],
            gutter="xl",
        ),
        modal,
    ],
    py="xl",
    )
)

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
    gene_name = str(gene.qualifiers.get("Alias", [""])[0])
    attributes = gene.qualifiers

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


def write_tmp(data, seqid, file_type, temp_dir):
    """Write sequence or GFF data to temporary files.
    
    Args:
        data: String for sequence data, DataFrame for GFF data
        seqid: Accession ID for the file name
        file_type: Either 'fa' or 'gff'
        temp_dir: Directory to write the temporary files
    """
    logger.debug(
        f"Entering write_tmp with seqid={seqid}, file_type={file_type}, temp_dir={temp_dir}"
    )

    file_path = os.path.join(temp_dir, f"{seqid}.{file_type}")
    logger.debug(f"File path set to {file_path}")

    try:
        if file_type == "fa":
            logger.debug("Writing sequence data to FASTA file")
            with open(file_path, "w") as f:
                f.write(f">{seqid}\n{data}\n")
        elif file_type == "gff":
            logger.debug("Writing GFF data")
            data.to_csv(file_path, sep="\t", header=False, index=False)
        else:
            logger.error(f"Unsupported file_type: {file_type}")
            raise ValueError(f"Unsupported file_type: {file_type}")

    except Exception as e:
        logger.error(f"An error occurred while writing the file: {e}")
        raise

    logger.debug(f"Successfully wrote file: {file_path}")
    return file_path


# def load_gff(accession):

#     df = cache.get(f"accession_gff_{accession}")
#     if df is None:
#         df = fetch_accession_gff(accession)

#     if df.empty:
#         logger.error(f"No GFF records found for accession: {accession}")
#     else:
#         # Identify the rows where 'end' is less than 'start'
#         mask = df["end"] < df["start"]

#         # Swap the 'start' and 'end' values where the mask is True
#         df.loc[mask, ["start", "end"]] = df.loc[mask, ["end", "start"]].values

#         # TODO: find a way to implement this
#         # # flip coordinates to favour captain order in starship
#         # df["priority"] = df["attributes"].str.contains("tyr").astype(int)

#         # df_sorted = df.sort_values(
#         #     by=["priority", "start", "end"], ascending=[False, True, True]
#         # ).drop(columns="priority")

#     return df


# def load_fa(accession):
#     df = cache.get("all_ships")
#     if df is None:
#         df = fetch_all_ships()

#     if isinstance(accession, str):
#         df = df[df["accession_tag"] == accession]
#     else:
#         df = df[df["accession_tag"].isin(accession)]

#     if df.empty:
#         logger.error(f"No fasta records found for accession: {accession}")

#     return df


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


def is_valid_sequence_file(file_path):
    valid_extensions = (".fa", ".fna", ".fasta", ".gb", ".gbk", ".gbff")
    return file_path.endswith(valid_extensions)


def multi_pgv(gff_files, seqs, tmp_file, len_thr=50, id_thr=30):
    gff_list = list(map(Gff, gff_files))
    gv = GenomeViz(track_align_type="center", fig_track_height=0.7)
    gv.set_scale_bar()

    # Store sequence sizes in a dictionary for validation
    seq_sizes = {}
    for gff in gff_list:
        for seqid, size in gff.get_seqid2size().items():
            seq_sizes[seqid] = size

    for gff in gff_list:
        for seqid, features in gff.get_seqid2features("gene").items():
            logger.info(f"Processing seqid: {seqid}")

            if seqid not in gff.get_seqid2size():
                logger.error(f"Error: SeqID {seqid} not found in GFF sizes")
                continue

            track = gv.add_feature_track(seqid, gff.get_seqid2size(), align_label=False)
            segment = track.get_segment(seqid)

            for idx, gene in enumerate(features):
                start = int(gene.location.start)
                end = int(gene.location.end)
                seg_size = gff.get_seqid2size()[seqid]

                # Validate coordinates
                if start < 0 or end > seg_size:
                    logger.error(f"Invalid coordinates: start={start}, end={end}, seg_size={seg_size}")
                    continue

                add_gene_feature(gene, segment, idx)

        for seq_file in seqs:
            if not os.path.isfile(seq_file):
                logger.error(f"Error: {seq_file} is not a valid file path.")
                return
            if os.path.getsize(seq_file) == 0:
                logger.error(f"Error: {seq_file} is an empty file.")
                return
            if not is_valid_sequence_file(seq_file):
                logger.error(f"Error: {seq_file} does not have a valid extension.")
                return

    # BLAST
    align_coords = Blast(seqs, seqtype="nucleotide").run()
    align_coords = AlignCoord.filter(
        align_coords, length_thr=len_thr, identity_thr=id_thr
    )

    # Run MMseqs RBH search
    # align_coords = MMseqs(seqs, threads=2).run()

    # Run MuMMer
    # align_coords = MUMmer(seqs).run()

    if len(align_coords) > 0:
        min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
        logger.info(f"Minimum identity: {min_ident}")

        color, inverted_color = "blue", "orange"
        for ac in align_coords:
            # Validate coordinates before adding link
            query_start, query_end = ac.query_link[2], ac.query_link[3]
            ref_start, ref_end = ac.ref_link[2], ac.ref_link[3]
            
            query_seqid = ac.query_link[1]
            ref_seqid = ac.ref_link[1]
            
            # Skip if coordinates are out of bounds
            if (query_start < 0 or query_end > seq_sizes[query_seqid] or 
                ref_start < 0 or ref_end > seq_sizes[ref_seqid]):
                logger.warning(
                    f"Skipping alignment with invalid coordinates: "
                    f"Query({query_start}, {query_end}) vs size {seq_sizes[query_seqid]}, "
                    f"Ref({ref_start}, {ref_end}) vs size {seq_sizes[ref_seqid]}"
                )
                continue
                
            logger.info(f"Adding link between {ac.query_link} and {ac.ref_link}")
            gv.add_link(
                ac.query_link,
                ac.ref_link,
                color=color,
                inverted_color=inverted_color,
                v=ac.identity,
                vmin=min_ident,
            )
        gv.set_colorbar([color, inverted_color], vmin=min_ident)
        message = None
    else:
        message = html.H4(
            f"No alignments found between ships (Length threshold: {len_thr}bp, ID threshold: {id_thr}%)"
        )

    fig = plot_legend(gv)
    gv.savefig_html(tmp_file, fig)

    return message


@callback(
    Output("pgv-table", "children"),
    Input("url", "href"),
    prevent_initial_call=False
)
@handle_callback_error
def load_ship_table(href):
    from src.database.sql_manager import fetch_ship_table
    from src.components.tables import make_pgv_table
    table_df = fetch_ship_table(curated=True)
    
    if table_df is not None and not table_df.empty:
        table = make_pgv_table(
            df=table_df,
            columns=table_columns,
            id="pgv-table",
            select_rows=True,
            pg_sz=10
        )
        logger.info("Table created successfully")
        return table
    
    logger.error("No data available or empty DataFrame")
    return html.Div("No data available")

@callback(
    [Output("pgv-figure", "children"), Output("pgv-message", "children")],
    Input("update-button", "n_clicks"),
    [
        State("pgv-table", "selected_rows"),
        State("pgv-table", "data"),
        State("length-threshold", "value"),
        State("identity-threshold", "value"),
    ],
)
@handle_callback_error
def update_pgv(n_clicks, selected_rows, table_data, len_thr, id_thr):
    from src.database.sql_manager import fetch_accession_ship
    message = None
    if not n_clicks:
        return no_update, "Select Starships from the table and click 'Show Selected Starships'"

    if n_clicks > 0:
        tmp_pgv = tempfile.NamedTemporaryFile(suffix=".html", delete=True).name
        if table_data and selected_rows is not None:
            try:
                if isinstance(selected_rows, list) and len(selected_rows) > 0:
                    with tempfile.TemporaryDirectory() as temp_dir:
                        tmp_gffs = []
                        tmp_fas = []
                        for idx in selected_rows:
                            row = table_data[idx]
                            accession = row["accession_tag"]
                            ship_data = fetch_accession_ship(accession)
                            fa_df = ship_data["sequence"]
                            tmp_fa = write_tmp(fa_df, accession, "fa", temp_dir)
                            tmp_fas.append(str(tmp_fa))

                            gff_df = ship_data["gff"]
                            tmp_gff = write_tmp(gff_df, accession, "gff", temp_dir)
                            tmp_gffs.append(tmp_gff)

                        if len(selected_rows) > 1 and len(selected_rows) <= 4:
                            message = multi_pgv(tmp_gffs, tmp_fas, tmp_pgv, len_thr, id_thr)
                        elif len(selected_rows) == 1:
                            single_pgv(tmp_gffs[0], tmp_pgv)
                        else:
                            output = html.P("Please select between 1 and 4 Starships.")
                            return output, None
                        try:
                            with open(tmp_pgv, "r") as file:
                                pgv_content = file.read()
                        except IOError:
                            output = html.P("Failed to read the temporary file.")

                        output = html.Iframe(
                            srcDoc=pgv_content,
                            style={
                                "width": "100%",
                                "height": "800px",
                                "border": "none",
                            },
                        )
                else:
                    output = html.H4("Please select at least one Starship.")
            except Exception as e:
                logger.error(f"Exception: {e}")
                output = html.H4(
                    "Error while comparing ships using BLAST. Try another combination."
                )
        else:
            output = html.H4("Select Starship(s) to visualize.")
    else:
        output = html.H4(
            "Select up to 4 Starships in the table above and click the button to visualize."
        )
    return (
        html.Div(
            [output],
            className="center-content text-center",
        ),
        message,
    )


toggle_modal = create_modal_callback(
    "pgv-table",
    "pgv-modal",
    "pgv-modal-content",
    "pgv-modal-title"
)