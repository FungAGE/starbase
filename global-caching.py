import copy
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
from flask_caching import Cache
import numpy as np
import os
import pandas as pd
import time

app = dash.Dash(__name__)

CACHE_CONFIG = {
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': '/tmp/dash_cache',
    'CACHE_DEFAULT_TIMEOUT': 300  # optional, but good to specify
}

# Initialize the cache with the Flask server instance
cache = Cache(app.server, config=CACHE_CONFIG)

N = 100

df = pd.DataFrame({
    'category': (
        (['apples'] * 5 * N) +
        (['oranges'] * 10 * N) +
        (['figs'] * 20 * N) +
        (['pineapples'] * 15 * N)
    )
})
df['x'] = np.random.randn(len(df['category']))
df['y'] = np.random.randn(len(df['category']))

app.layout = html.Div([
    html.Link(href="https://codepen.io/chriddyp/pen/bWLwgP.css", rel="stylesheet"),
    html.Link(href="https://codepen.io/chriddyp/pen/brPBPO.css", rel="stylesheet"),
    dcc.Dropdown(
        id='dropdown',
        options=[{'label': i, 'value': i} for i in df['category'].unique()],
        value='apples'
    ),
    html.Div([
        html.Div(dcc.Graph(id='graph-1'), className="six columns"),
        html.Div(dcc.Graph(id='graph-2'), className="six columns"),
    ], className="row"),
    html.Div([
        html.Div(dcc.Graph(id='graph-3'), className="six columns"),
        html.Div(dcc.Graph(id='graph-4'), className="six columns"),
    ], className="row"),

    # hidden signal value
    html.Div(id='signal', style={'display': 'none'})
])


@cache.memoize()
def global_store(value):
    print('Computing value with {}'.format(value))
    time.sleep(5)
    return df[df['category'] == value]


def generate_figure(value, figure):
    if value is None:
        return {}

    filtered_dataframe = global_store(value)
    figure['data'][0]['x'] = filtered_dataframe['x']
    figure['data'][0]['y'] = filtered_dataframe['y']
    figure['layout'] = {'margin': {'l': 20, 'r': 10, 'b': 20, 't': 10}}
    return figure


@app.callback(Output('signal', 'children'), [Input('dropdown', 'value')])
def compute_value(value):
    global_store(value)
    return value


@app.callback(Output('graph-1', 'figure'), [Input('signal', 'children')])
def update_graph_1(value):
    return generate_figure(value, {
        'data': [{
            'type': 'scatter',
            'mode': 'markers',
            'marker': {
                'opacity': 0.5,
                'size': 14,
                'line': {'border': 'thin darkgrey solid'}
            }
        }]
    })


@app.callback(Output('graph-2', 'figure'), [Input('signal', 'children')])
def update_graph_2(value):
    return generate_figure(value, {
        'data': [{
            'type': 'scatter',
            'mode': 'lines',
            'line': {'shape': 'spline', 'width': 0.5},
        }]
    })


@app.callback(Output('graph-3', 'figure'), [Input('signal', 'children')])
def update_graph_3(value):
    return generate_figure(value, {
        'data': [{
            'type': 'histogram2d',
        }]
    })


@app.callback(Output('graph-4', 'figure'), [Input('signal', 'children')])
def update_graph_4(value):
    return generate_figure(value, {
        'data': [{
            'type': 'histogram2dcontour',
        }]
    })


if __name__ == '__main__':
    app.run_server(debug=True, threaded=True, host='localhost', port=8005)
