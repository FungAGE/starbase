import dash
from dash import html
import dash_bio

dash.register_page(__name__)

HOSTED_GENOME_DICT = [{"value": "altalt8_s00067", "label": "altalt8_s00067"}]

# layout = html.Div(
#     [
#         html.P("Select the Starship to display below."),
#         dcc.Dropdown(
#             id="default-starship-select",
#             options=HOSTED_GENOME_DICT,
#             value="altalt8_s00067",
#         ),
#         html.Hr(),
#         dcc.Loading(id="default-container"),
#     ]
# )


# # Return the IGV component with the selected starship.
# @callback(
#     Output("default-container", "children"), Input("default-starship-select", "value")
# )
# def return_jbrowse(starship):
#     with open("src/tmp/" + starship + ".json", "r") as file:
#         data = json.load(file)
#     my_assembly = data["assembly"]
#     # print(my_assembly)
#     my_tracks = data["tracks"]
#     my_location = data["location"]
#     my_default_session = data["defaultSession"]
#     return html.Div(
#         [
#             dash_jbrowse.LinearGenomeView(
#                 id="lgv",
#                 assembly=my_assembly,
#                 tracks=my_tracks,
#                 defaultSession=my_default_session,
#                 location=my_location,
#             ),
#         ]
#     )


text_style = {"color": "#506784", "font-family": "Open Sans"}

_COMPONENT_ID = "igv-chart"


def header_colors():
    return {
        "bg_color": "#0F5BA7",
        "font_color": "white",
    }


layout = html.Div(
    [
        # dash_bio.Igv(
        #     id="reference-igv",
        #     reference={
        #         "id": "ASM985889v3",
        #         "name": "Sars-CoV-2 (ASM985889v3)",
        #         "fastaURL": "https://s3.amazonaws.com/igv.org.genomes/covid_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna",
        #         "indexURL": "https://s3.amazonaws.com/igv.org.genomes/covid_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.fai",
        #         "order": 1000000,
        #         "tracks": [
        #             {
        #                 "name": "Annotations",
        #                 "url": "https://s3.amazonaws.com/igv.org.genomes/covid_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz",
        #                 "displayMode": "EXPANDED",
        #                 "nameField": "gene",
        #                 "height": 150,
        #                 "color": "rgb(176,141,87)",
        #             }
        #         ],
        #     },
        # )
        dash_bio.Igv(
            id="igv-chart",
            reference={
                "id": "altalt8_s00067",
                "name": "altalt8_s00067",
                "fastaURL": "src/tmp/altalt8_s00067.fa",
                "indexURL": "src/tmp/altalt8_s00067.fa.fai",
                # "order": 1000000,
                "tracks": [
                    {
                        "name": "Annotations",
                        "url": "src/tmp/altalt8_s00067.gff.gz",
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
