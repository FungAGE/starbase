import dash
from dash import dcc, html, callback
from dash.dependencies import Output, Input
import dash_bio as dashbio

import os
import pandas as pd

dash.register_page(__name__)

df = pd.read_csv("src/assets/joined_ships.csv")


def create_genome_file_dict(directory):
    genome_file_dict = {}
    for genome_dir in os.listdir(directory):
        genome_path = os.path.join(directory, genome_dir)
        if os.path.isdir(genome_path):
            fasta_files = [f for f in os.listdir(genome_path) if f.endswith(".fasta")]
            fasta_index_files = [
                f for f in os.listdir(genome_path) if f.endswith(".fasta.fai")
            ]
            gff_files = [f for f in os.listdir(genome_path) if f.endswith(".gff")]

            if fasta_files and fasta_index_files and gff_files:
                genome_file_dict[genome_dir] = {
                    "fasta": os.path.join(genome_path, fasta_files[0]),
                    "fasta_index": os.path.join(genome_path, fasta_index_files[0]),
                    "gff": os.path.join(genome_path, gff_files[0]),
                }
    return genome_file_dict


# Replace 'directory_path' with the path to the directory containing genome subdirectories
directory_path = "data/Starships"
genome_file_dict = create_genome_file_dict(directory_path)
print(genome_file_dict)


HOSTED_GENOME_DICT = [
    {"value": "mm10", "label": "Mouse (GRCm38/mm10)"},
    {"value": "rn6", "label": "Rat (RGCS 6.0/rn6)"},
    {"value": "gorGor4", "label": "Gorilla (gorGor4.1/gorGor4)"},
    {"value": "panTro4", "label": "Chimp (SAC 2.1.4/panTro4)"},
    {"value": "panPan2", "label": "Bonobo (MPI-EVA panpan1.1/panPan2)"},
    {"value": "canFam3", "label": "Dog (Broad CanFam3.1/canFam3)"},
    {"value": "ce11", "label": "C. elegans (ce11)"},
]

layout = html.Div(
    [
        dcc.Loading(id="default-igv-container"),
        html.Hr(),
        html.P("Select the genome to display below."),
        dcc.Dropdown(
            id="default-igv-genome-select", options=HOSTED_GENOME_DICT, value="ce11"
        ),
    ]
)


# Return the IGV component with the selected genome.
@callback(
    Output("default-igv-container", "children"),
    Input("default-igv-genome-select", "value"),
)
def return_igv(genome):
    return html.Div(
        [
            dashbio.Igv(
                id="reference-igv",
                reference={
                    "id": "ASM985889v3",
                    "name": "Sars-CoV-2 (ASM985889v3)",
                    "fastaURL": "https://s3.amazonaws.com/igv.org.genomes/covid_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna",
                    "indexURL": "https://s3.amazonaws.com/igv.org.genomes/covid_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.fai",
                    "order": 1000000,
                    "tracks": [
                        {
                            "name": "Annotations",
                            "url": "https://s3.amazonaws.com/igv.org.genomes/covid_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz",
                            "displayMode": "EXPANDED",
                            "nameField": "gene",
                            "height": 150,
                            "color": "rgb(176,141,87)",
                        }
                    ],
                },
            )
        ]
    )


if __name__ == "__main__":
    app.run(debug=True)
