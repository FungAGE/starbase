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
from Bio import SeqIO
from jinja2 import Template

from src.components.cache import cache
from src.components.sql_engine import starbase_engine
from src.components.tables import make_ship_table
from src.components.cache_manager import load_from_cache
from src.components.sql_queries import (
    fetch_all_ships,
    fetch_accession_gff,
    fetch_ship_table,
)
from src.components.callbacks import create_accession_modal


logger = logging.getLogger(__name__)

dash.register_page(__name__)

ColorCycler.set_cmap("tab10")

table_columns = [
    {
        "name": "Accession",
        "id": "accession_tag",
        "deletable": False,
        "selectable": True,
    },
    {
        "name": "Starship Family",
        "id": "familyName",
        "deletable": False,
        "selectable": True,
    },
    {
        "name": "Species",
        "id": "species",
        "deletable": False,
        "selectable": True,
    },
]

modal =     dbc.Modal(
            [
                dbc.ModalHeader(dbc.ModalTitle(id="pgv-modal-title")),
                dbc.ModalBody(id="pgv-modal-content"),
            ],
            id="pgv-modal",
            is_open=False,
)

layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        modal,
        dmc.Grid(
            justify="center",
            align="top",
            children=[
                dmc.GridCol(
                    style={"padding": "40px"},
                    span={
                        "sm": 12,
                        "lg": 8,
                    },
                    children=[
                        dcc.Loading(
                            id="loading",
                            type="circle",
                            children=[
                                html.Div(id="pgv-table", className="center-content"),
                            ],
                        ),
                    ],
                ),
                dmc.GridCol(
                    span={
                        "sm": 12,
                        "lg": 8,
                    },
                    children=[
                        dmc.Center(
                            [
                                dmc.Button(
                                    "Show Starship(s) in Viewer",
                                    id="update-button",
                                    n_clicks=0,
                                    className="text-custom text-custom-sm text-custom-md",
                                ),
                            ]
                        ),
                    ],
                ),
                dmc.GridCol(
                    style={"padding": "40px"},
                    span={
                        "sm": 12,
                        "lg": 8,
                    },
                    children=[
                        dcc.Loading(
                            id="loading-1",
                            type="circle",
                            children=[
                                html.Div(id="pgv-message"),
                                html.Div(
                                    id="pgv-figure",
                                    style={
                                        "height": "800px",
                                        "width": "100%",
                                        "overflow": "auto",
                                        "border": "1px solid #ccc",
                                    },
                                ),
                            ],
                        )
                    ],
                ),
            ],
        ),
    ],
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


def write_tmp(df, seqid, file_type=None, temp_dir=None):
    logger.debug(
        "Entering write_tmp with seqid=%s, file_type=%s, temp_dir=%s",
        seqid,
        file_type,
        temp_dir,
    )

    if temp_dir is None:
        logger.warning("temp_dir is None; this may cause errors in file path creation.")

    file_path = os.path.join(temp_dir, f"{seqid}.{file_type}")
    logger.debug("File path set to %s", file_path)

    try:
        if file_type == "gff":
            logger.debug("file_type is 'gff'. Attempting to write dataframe to file.")
            df.to_csv(file_path, sep="\t", header=False, index=False)
            logger.debug("Dataframe written to %s", file_path)
        elif file_type == "fa":
            logger.debug(
                "file_type is 'fa'. Checking for 'sequence' column and first entry."
            )

            if "sequence" not in df.columns:
                logger.error("'sequence' column is missing from DataFrame.")
                raise KeyError("'sequence' column is missing from DataFrame.")

            if df.empty or df["sequence"].empty:
                logger.error("The 'sequence' column is empty.")
                raise ValueError("The 'sequence' column is empty.")

            seq = df["sequence"].iloc[0]  # Use .iloc[0] to avoid KeyError
            with open(file_path, "w") as f:
                f.write(f">{seqid}\n{seq}\n")
            logger.debug("Sequence data written to %s", file_path)
        else:
            logger.error("Unsupported file_type provided: %s", file_type)
            raise ValueError("Unsupported file_type: %s" % file_type)

    except Exception as e:
        logger.exception("An error occurred while writing the file: %s", e)
        raise

    logger.debug("Exiting write_tmp with file_path=%s", file_path)
    return file_path


def load_gff(accession):

    df = load_from_cache(f"accession_gff_{accession}")
    if df is None:
        df = fetch_accession_gff(accession)

    if df.empty:
        logger.error(f"No GFF records found for accession: {accession}")
    else:
        # Identify the rows where 'end' is less than 'start'
        mask = df["end"] < df["start"]

        # Swap the 'start' and 'end' values where the mask is True
        df.loc[mask, ["start", "end"]] = df.loc[mask, ["end", "start"]].values

        # TODO: find a way to implement this
        # # flip coordinates to favour captain order in starship
        # df["priority"] = df["attributes"].str.contains("tyr").astype(int)

        # df_sorted = df.sort_values(
        #     by=["priority", "start", "end"], ascending=[False, True, True]
        # ).drop(columns="priority")

    return df


def load_fa(accession):
    df = load_from_cache("all_ships")
    if df is None:
        df = fetch_all_ships()

    if isinstance(accession, str):
        df = df[df["accession_tag"] == accession]
    else:
        df = df[df["accession_tag"].isin(accession)]

    if df.empty:
        logger.error(f"No fasta records found for accession: {accession}")

    return df


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


def multi_pgv(gff_files, seqs, tmp_file):
    gff_list = list(map(Gff, gff_files))
    gv = GenomeViz(track_align_type="center", fig_track_height=0.7)
    gv.set_scale_bar()

    for gff in gff_list:
        # Check GFF SeqIDs
        for seqid, features in gff.get_seqid2features("gene").items():
            logger.info(f"Processing seqid: {seqid}")

            if seqid not in gff.get_seqid2size():
                logger.error(f"Error: SeqID {seqid} not found in GFF sizes")
                continue

            track = gv.add_feature_track(seqid, gff.get_seqid2size(), align_label=False)
            segment = track.get_segment(seqid)

            for idx, gene in enumerate(features):
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
    len_thr = 50
    id_thr = 30
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
)
def load_ship_table(href):
    table_df = load_from_cache("ship_table")
    if table_df is None:
        table_df = fetch_ship_table()
    if href:
        table = make_ship_table(
            df=table_df, columns=table_columns, id="pgv-table", pg_sz=15
        )
        return table


@callback(
    [Output("pgv-figure", "children"), Output("pgv-message", "children")],
    Input("update-button", "n_clicks"),
    [
        State("pgv-table", "derived_virtual_selected_rows"),
        State("pgv-table", "derived_virtual_data"),
    ],
)
def update_pgv(n_clicks, selected_rows, table_data):
    message = None
    if n_clicks > 0:
        tmp_pgv = tempfile.NamedTemporaryFile(suffix=".html", delete=True).name
        if table_data and selected_rows is not None:
            try:
                table_df = pd.DataFrame(table_data)
                logger.info(f"Columns in table_df: {table_df.columns.tolist()}")
                logger.info(f"Selected rows: {selected_rows}")

                if isinstance(selected_rows, list) and all(
                    isinstance(idx, int) for idx in selected_rows
                ):
                    try:
                        rows = table_df.iloc[selected_rows]
                    except IndexError as e:
                        return html.Div("Index out of bounds")

                    with tempfile.TemporaryDirectory() as temp_dir:
                        tmp_gffs = []
                        tmp_fas = []
                        for index, row in rows.iterrows():
                            accession = row["accession_tag"]

                            logger.info(f"Fetching FA for accession: {accession}")
                            fa_df = load_fa(accession)
                            tmp_fa = write_tmp(fa_df, accession, "fa", temp_dir)
                            tmp_fas.append(str(tmp_fa))

                            logger.info(f"Fetching GFF for accession: {accession}")
                            gff_df = load_gff(accession)

                            tmp_gff = write_tmp(gff_df, accession, "gff", temp_dir)
                            tmp_gffs.append(tmp_gff)

                            output = html.P("Select up to four Starships to compare.")
                        if len(selected_rows) > 1 and len(selected_rows) <= 4:
                            message = multi_pgv(tmp_gffs, tmp_fas, tmp_pgv)
                        elif len(selected_rows) == 1:
                            single_pgv(tmp_gffs[0], tmp_pgv)
                        else:
                            output = html.P("No valid selection.")
                        try:
                            with open(tmp_pgv, "r") as file:
                                pgv_content = file.read()
                        except IOError:
                            output = html.P("Failed to read the temporary file.")

                        output = html.Iframe(
                            srcDoc=pgv_content,
                            style={
                                "width": "100%",
                                "height": "100%",
                                "border": "none",
                            },
                        )
                else:
                    output = html.H4("Invalid row selection.")
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


@callback(
    Output("pgv-modal", "is_open"),
    Output("pgv-modal-content", "children"),
    Output("pgv-modal-title", "children"),
    Output("pgv-table", "active_cell"),
    Input("pgv-table", "active_cell"),
    State("pgv-modal", "is_open"),
    State("pgv-table", "data"),
)
def toggle_modal(cell_clicked, is_open, table_data):    
    # If no cell was clicked, keep modal closed
    if cell_clicked is None:
        return False, no_update, no_update, no_update
        
    if table_data:
        try:
            row = cell_clicked["row"]
            row_data = table_data[row]
            accession = row_data.get("accession_tag")
            if accession:
                modal_content, modal_title = create_accession_modal(accession)
                return True, modal_content, modal_title, None
            else:
                return False, "No accession data found", "Error", None
                
        except Exception as e:
            logger.error(f"Error in toggle_modal: {str(e)}")
            return False, "Error loading modal", "Error", None
            
    return is_open, no_update, no_update, no_update
