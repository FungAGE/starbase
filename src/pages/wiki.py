import dash
import dash_bootstrap_components as dbc
from dash import html


dash.register_page(__name__)

layout = html.Div(
    [
        html.Div(
            style={
                "display": "flex",
                "justify-content": "center",
                "align-items": "center",
            },
            children=[
                html.Div(
                    [
                        html.Div(
                            [
                                dbc.Card(
                                    dbc.CardBody(
                                        [
                                            html.H3("What is a Starship?"),
                                            html.P(
                                                "Starships are novel family of class II DNA transposons, endemic to Pezizomycotina. Starships can be extremely large (~20-700kb), making up to 2% of fungal genomes. These elements replicate within the host genome via tyrosine recombinases (captain genes) [2]. They can also pick up and carry relevant genetic 'cargo', including genes for metal resistance in Paecilomyces, cheese making in Penicillium, and enable the reansfer of formaldehyde resistance in Aspergillus nidulans and Penicillium chrysogenum."
                                            ),
                                            html.Br(),
                                            html.Div(
                                                style={
                                                    "display": "flex",
                                                    "justify-content": "center",
                                                    "align-items": "center",
                                                    "backgroundColor": "white",
                                                },
                                                children=[
                                                    html.Img(
                                                        src="assets/images/starship-model.png",
                                                        width="85%",
                                                    )
                                                ],
                                                className="box-body",
                                            ),
                                        ]
                                    ),
                                    color="primary",
                                    inverse=True,
                                )
                            ]
                        ),
                    ],
                    className="box",
                    style={"width": "75%"},
                )
            ],
            className="row",
        )
    ]
)
