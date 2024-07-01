import dash
from dash import html

dash.register_page(__name__, title="Home", name="Home", path="/")

index_string = """
<!DOCTYPE html>
<html>
    <head>
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <meta charset="UTF-8">
        <title>Dash App</title>
        <link href="/assets/output.css" rel="stylesheet">
    </head>
    <body class="bg-gray-100">
        <div id="react-entry-point">
            <script type="text/javascript">
                window.PlotlyConfig = {MathJaxConfig: 'local'};
            </script>
        </div>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
"""


def mod_home_ui():
    working = [
        "Catalogue/Wiki of Starship Metadata",
        "Submission of new Starship sequences",
        "BLAST/HMMER searches",
    ]
    working_ul = html.Ul(
        [html.Li(item) for item in working],
        style={
            "fontSize": "0.6vw",
        },
    )

    not_working = [
        "Synteny/Genome Browser",
        html.P(
            [
                html.Span(
                    "starfish",
                    className="logo-text",
                ),
                " webserver",
            ],
        ),
    ]
    not_working_ul = html.Ul(
        [html.Li(item) for item in not_working],
        style={
            "fontSize": "0.6vw",
        },
    )

    return html.Div(
        className="container mx-auto p-4",
        children=[
            html.Div(
                className="flex flex-col items-start md:items-center",
                children=[
                    html.Div(
                        className="w-full md:w-auto",
                        children=[
                            html.H1(
                                children=[
                                    html.Span(
                                        "starbase",
                                        className="logo-text",
                                        style={"fontSize": "7vw"},
                                    ),
                                    html.Br(),
                                    "A database and toolkit for exploring large eukaryotic transposable elements in Fungi",
                                ],
                                className="text-left md:text-center text-2xl md:text-4xl",
                            ),
                            html.Div(
                                className="flex flex-col md:flex-row",
                                children=[
                                    html.Div(
                                        className="w-full md:w-1/2 p-2",
                                        children=[
                                            html.Div(
                                                className="bg-white shadow rounded-lg",
                                                children=[
                                                    html.Div(
                                                        className="bg-gray-200 p-4 rounded-t-lg",
                                                        children=[
                                                            html.H4(
                                                                children=[
                                                                    "What can I currently use ",
                                                                    html.Span(
                                                                        "starbase",
                                                                        className="logo-text",
                                                                    ),
                                                                    " for?",
                                                                ],
                                                                className="text-left md:text-center text-sm md:text-base",
                                                            )
                                                        ],
                                                    ),
                                                    html.Div(
                                                        className="p-4",
                                                        children=[working_ul],
                                                    ),
                                                ],
                                            )
                                        ],
                                    ),
                                    html.Div(
                                        className="w-full md:w-1/2 p-2",
                                        children=[
                                            html.Div(
                                                className="bg-white shadow rounded-lg",
                                                children=[
                                                    html.Div(
                                                        className="bg-gray-200 p-4 rounded-t-lg",
                                                        children=[
                                                            html.H4(
                                                                children=[
                                                                    "Functions of ",
                                                                    html.Span(
                                                                        "starbase",
                                                                        className="logo-text",
                                                                    ),
                                                                    " under active development:",
                                                                ],
                                                                className="text-left md:text-center text-sm md:text-base",
                                                            )
                                                        ],
                                                    ),
                                                    html.Div(
                                                        className="p-4",
                                                        children=[not_working_ul],
                                                    ),
                                                ],
                                            )
                                        ],
                                    ),
                                ],
                            ),
                            html.Div(
                                className="flex justify-center my-4",
                                children=[
                                    html.Img(
                                        src="assets/images/starbase-map.png",
                                        className="w-full",
                                    )
                                ],
                            ),
                            html.Div(
                                className="bg-white shadow rounded-lg",
                                children=[
                                    html.Div(
                                        className="bg-gray-200 p-4 rounded-t-lg",
                                        children=[
                                            html.H4(
                                                "Data Availability",
                                                className="text-left md:text-center text-lg md:text-2xl",
                                            )
                                        ],
                                    ),
                                    html.Div(
                                        className="p-4",
                                        children=[
                                            html.P(
                                                children=[
                                                    "We have been maintaining ",
                                                    html.Span(
                                                        "starbase",
                                                        className="logo-text",
                                                    ),
                                                    " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export. In the mean time, you can retrieve all Starship sequences, annotations, and more, in a single .zip file (size ~100Mb)",
                                                ],
                                                className="text-left md:text-center text-xs md:text-sm",
                                            ),
                                            html.Div(
                                                className="text-center mt-4",
                                                children=[
                                                    html.Button(
                                                        children=[
                                                            html.P(
                                                                children=[
                                                                    "Download the latest version of ",
                                                                    html.Span(
                                                                        "starbase",
                                                                        className="logo-text",
                                                                    ),
                                                                    ".",
                                                                ],
                                                                className="text-left md:text-center text-xs md:text-sm",
                                                            )
                                                        ],
                                                        id="dl_package",
                                                        className="bg-blue-500 hover:bg-blue-700 text-white font-bold py-2 px-4 rounded",
                                                    )
                                                ],
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                        ],
                    ),
                ],
            ),
        ],
    )


layout = mod_home_ui
