import dash
import dash_jbrowse
import json
import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html, callback
from dash.dependencies import Output, Input, State

dash.register_page(__name__)

HOSTED_GENOME_DICT = [
    {"value": "altalt8_s00067", "label": "altalt8_s00067"}
]

layout = html.Div(
    [
        html.P("Select the Starship to display below."),
        dcc.Dropdown(
            id="default-starship-select", options=HOSTED_GENOME_DICT, value="altalt8_s00067"
        ),
        html.Hr(),
        dcc.Loading(id="default-container"),
    ]
)


# Return the IGV component with the selected starship.
@callback(
    Output("default-container", "children"), Input("default-starship-select", "value")
)
def return_jbrowse(starship):
    with open('src/assets/' + starship + '.json', "r") as file:
            data = json.load(file)
    my_assembly = data["assembly"]
    # print(my_assembly)
    my_tracks = data["tracks"]
    my_location = data["location"]
    my_default_session = data["defaultSession"]
    return html.Div(
        [
            dash_jbrowse.LinearGenomeView(
                id='lgv',
                assembly=my_assembly,
                tracks=my_tracks,
                defaultSession=my_default_session,
                location=my_location,
            ),
        ]
    )

# text_style = {
#     'color': "#506784",
#     'font-family': 'Open Sans'
# }

# _COMPONENT_ID = 'igv-chart'

# def description():
#     return 'A high-performance genomics viewer with an interactive UI and support for a ' \
#            'wide variety of data types and features.'


# def header_colors():
#     return {
#         'bg_color': '#0F5BA7',
#         'font_color': 'white',
#     }


# layout = html.Div(id='igv-body', className='app-body', children=[
#         html.Div(id='igv-control-tabs', className='control-tabs', children=[
#             dcc.Tabs(
#                 id='igv-tabs',
#                 value='what-is',
#                 children=[
#                     dcc.Tab(
#                         label='About',
#                         value='what-is',
#                         children=html.Div(className='control-tab', children=[
#                             html.H4(className='what-is', children='What is IGV?'),
#                             dcc.Markdown(
#                                 """
#                                 The Dash IGV component is a high-performance genomics
#                                 data visualization component developed originally by the [IGV
#                                 Team](https://igv.org/) based at UC San Diego and the Broad
#                                 Institute. It offers
#                                 support for array-based and next-generation sequencing data,
#                                 and a smooth, interactive UI for real-time exploration of large
#                                 scale genomic data. This includes visualizing alignments,
#                                 copy number,
#                                 genome-wide interactions, gene expression and methylation, and more
#                                 data types. Data tracks, interactions, and analysis can be
#                                 controlled
#                                 by integrating with a Dash app to create a complete dynamic
#                                 workflow.


#                                 Read more about the component here:
#                                 https://github.com/igvteam/igv.js/
#                                 """
#                             )
#                         ])
#                     ),
#                     dcc.Tab(
#                         label='Data',
#                         value='data',
#                         children=html.Div(className='control-tab', children=[
#                             html.Div(className='app-controls-block', children=[
#                                 html.Div(
#                                     className='fullwidth-app-controls-name',
#                                     children="Select a Genome"
#                                 ),
#                                 dcc.Dropdown(
#                                     id='genome-dropdown',
#                                     options=HOSTED_GENOME_DICT,
#                                     value='rn6',
#                                 ),
#                                 html.Div(
#                                     className='app-controls-desc',
#                                     children='Select a Genome Identifier to display the remotely '
#                                              'hosted '
#                                              'genome.'
#                                 ),
#                             ]),
#                             html.Hr(
#                                 className='igv-separator'
#                             ),
#                             html.Div(
#                                 className='app-controls-block',
#                                 children=[
#                                     html.Div(className='app-controls-name',
#                                              children='Minimum Window Size'),
#                                     dcc.Slider(
#                                         className='control-slider',
#                                         id='minimum-bases',
#                                         value=100,
#                                         min=10,
#                                         max=200,
#                                         step=10,
#                                         marks=dict((i, str(i)) for i in range(10, 190, 30))
#                                     ),

#                                     html.Div(
#                                         className='app-controls-desc',
#                                         children='Minimum window size in base pairs when zooming '
#                                                  'in.'
#                                     ),
#                                 ],
#                             ),
#                         ])
#                     )
#                 ]
#             )
#         ]),
#         dcc.Loading(parent_className='dashbio-loading', id='igv-output'),
#     ])


# # Return the IGV component with the selected genome and base length
# @callback(
#     Output('igv-output', 'children'),
#     [Input('genome-dropdown', 'value'),
#       Input('minimum-bases', 'value')]
# )
# def return_igv(genome, bases):
#     return (
#         html.Div([
#             dash_bio.Igv(
#                 id=_COMPONENT_ID,
#                 reference={
#                 'id': 'altalt8_s00067',
#                 'name': 'altalt8_s00067',
#                 'fastaURL': 'src/assets/altalt8_s00067.fa',
#                 'indexURL': 'src/assets/altalt8_s00067.fa.fai',
#                 'tracks': [
#                     {
#                         'name': 'Annotations',
#                         'url': 'src/assets/altalt8_s00067.gff.gz',
#                         'displayMode': 'EXPANDED',
#                         'nameField': 'gene',
#                         'height': 150,
#                         'color': 'rgb(176,141,87)'
#                     }
#                 ]
#             },
#                 minimumBases=bases,
#             )
#         ])
#     )