import dash
import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html, callback
from dash.dependencies import Output, Input, State

import dash_ag_grid as dag
import base64

import pandas as pd

df = pd.read_csv('src/assets/joined_ships.csv')
ship_count = df[["starshipID"]].nunique()
species = df['genus'] + '-' + df['species']
species_count = species.nunique()

dash.register_page(__name__)

layout = html.Div([
    html.Div([
        html.Div([
            html.Div([
                html.H3("What is a Starship?"),
                html.P("Starships are novel family of class II DNA transposons, endemic to Pezizomycotina. Starships can be extremely large (~20-700kb), making up to 2% of fungal genomes. These elements replicate within the host genome via tyrosine recombinases (captain genes) [2]. They can also pick up and carry relevant genetic 'cargo', including genes for metal resistance in Paecilomyces, cheese making in Penicillium, and enable the reansfer of formaldehyde resistance in Aspergillus nidulans and Penicillium chrysogenum.")
            ], className="box-body"),
            html.Div([
                html.Img(src="assets/images/starship-model.png", width="85%")
            ], className="box-body")
        ], className="box", style={"width": "100%"})
    ], className="row"),
    html.Br(),
    dbc.Row([
        dbc.Col(
            dbc.Card(
                dbc.CardBody([
                    html.H4("Total number of Starships in starbase", className="card-title"),
                    html.P(ship_count, className="card-text"),
                ]),
                className="shadow-sm"
            ),
            width=4
        ),
        dbc.Col(
            dbc.Card(
                dbc.CardBody([
                    html.H4("Fungal species with Starships", className="card-title"),
                    html.P(species_count, className="card-text"),
                ]),
                className="shadow-sm"
            ),
            width=4
        )
    ])
])