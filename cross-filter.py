import dash
from dash.dependencies import Input, Output
from dash import dcc, html
import numpy as np
import pandas as pd

external_stylesheets = [
    "https://codepen.io/chriddyp/pen/bWLwgP.css",
    ]

app = dash.Dash(__name__,
    external_stylesheets=external_stylesheets)

np.random.seed(0)
df = pd.DataFrame({
    'Column {}'.format(i): np.random.rand(30) + i * 10
    for i in range(6)})

app.layout = html.Div([
    html.Div(
        dcc.Graph(
            id='g1',
            config={'displayModeBar': False}
        ), className='four columns'
    ),
    html.Div(
        dcc.Graph(
            id='g2',
            config={'displayModeBar': False}
        ), className='four columns'),
    html.Div(
        dcc.Graph(
            id='g3',
            config={'displayModeBar': False}
        ), className='four columns')
], className='row')


def highlight(x, y):
    def callback(*selectedDatas):
        # Initialize selected points as the full dataframe index
        selected_points = df.index

        # Cross-filter: for each selected data (from the graphs)
        for selected_data in selectedDatas:
            if selected_data and selected_data['points']:
                selected_index = [
                    p['customdata'] for p in selected_data['points']
                ]
                if selected_index:
                    # Use intersection to only keep selected points across graphs
                    selected_points = np.intersect1d(selected_points, selected_index)
        
        # Create figure structure
        figure = {
            'data': [
                {
                    'x': df[x],
                    'y': df[y],
                    'text': df.index,
                    'customdata': df.index,
                    'type': 'scatter',
                    'mode': 'markers+text',
                    'marker': {
                        'color': 'rgba(0, 116, 217, 0.7)',
                        'size': 12,
                        'line': {
                            'color': 'rgb(0, 116, 217)',
                            'width': 0.5
                        }
                    },
                    'selected_points': selected_points,  # This ensures cross-filtering
                    'unselected': {
                        'marker': {
                            'opacity': 0.3,
                        },
                        'textfont': {
                            'color': 'rgba(0, 0, 0, 0)'  # Transparent text for unselected points
                        }
                    }
                },
            ],
            'layout': {
                'margin': {'l': 15, 'r': 0, 'b': 15, 't': 5},
                'dragmode': 'select',
                'hovermode': 'closest',
                'showlegend': False
            }
        }

        # Highlight the selected region with a rectangle
        if selectedDatas[0] and selectedDatas[0].get('range'):
            shape = {
                'type': 'rect',
                'line': {
                    'width': 1,
                    'dash': 'dot',
                    'color': 'darkgrey'
                },
                'x0': selectedDatas[0]['range']['x'][0],
                'x1': selectedDatas[0]['range']['x'][1],
                'y0': selectedDatas[0]['range']['y'][0],
                'y1': selectedDatas[0]['range']['y'][1]
            }
            figure['layout']['shapes'] = [shape]
        else:
            figure['layout']['shapes'] = []

        return figure

    return callback


# Callbacks for updating figures with cross-filtering
app.callback(
    Output('g1', 'figure'),
    [Input('g1', 'selectedData'),
     Input('g2', 'selectedData'),
     Input('g3', 'selectedData')]
)(highlight('Column 0', 'Column 1'))

app.callback(
    Output('g2', 'figure'),
    [Input('g2', 'selectedData'),
     Input('g1', 'selectedData'),
     Input('g3', 'selectedData')]
)(highlight('Column 2', 'Column 3'))

app.callback(
    Output('g3', 'figure'),
    [Input('g3', 'selectedData'),
     Input('g1', 'selectedData'),
     Input('g2', 'selectedData')]
)(highlight('Column 4', 'Column 5'))

if __name__ == '__main__':
    app.run_server(debug=True, threaded=True, host='localhost', port=8005)
