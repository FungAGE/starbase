from datetime import datetime

from dash import dash
from dash import html, dcc
from dash.dependencies import Input, Output

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, suppress_callback_exceptions=True, external_stylesheets=external_stylesheets)

IDS = [1, 2, 3]


def divs_list():
    return [html.Div([
        dcc.Markdown(
            '',
            id='model-{}-markdown'.format(id)
        ),
        html.P(
            '',
            id='model-{}-p'.format(id)
        ),
        html.Button(
            'Delete',
            id='model-{}-delete-button'.format(id),
            style={'width': '49%'}
        ),
        html.Button(
            'Start/Stop',
            id='model-{}-toggle-button'.format(id),
            style={'marginLeft': '2%', 'width': '49%'}
        ),

        html.Hr()
    ], id='model-k2-{}'.format(id)) for id in IDS]


def create_callbacks():
    def delete(i):
        def callback(*data):
            return """# Delete
            {}""".format(datetime.now().time().isoformat())

        return callback

    def toggle(i):
        def callback(*data):
            return "Toggle {}".format(datetime.now().time().isoformat())

        return callback

    for model_id in IDS:
        app.callback(
            Output('model-{}-markdown'.format(model_id), 'children'),
            [Input('model-{}-delete-button'.format(model_id), 'n_clicks')]
        )(delete(model_id))

        app.callback(
            Output('model-{}-p'.format(model_id), 'children'),
            [Input('model-{}-toggle-button'.format(model_id), 'n_clicks')]
        )(toggle(model_id))


def layout():
    divs = divs_list()

    return html.Div([  # Wrap the return value in html.Div
        html.Div([
            html.H1('Example'),
            html.Div(divs, id='model-k2-divs'),
        ], style={'maxWidth': '720px', 'margin': '0 auto'}, id='example')
    ])


app.layout = html.Div(children=[
    html.Meta(
        name='viewport',
        content='width=device-width, initial-scale=1.0'
    ),

    html.Div([
        html.Div([
            html.Div(
                html.Div(id="page", className="content"),
                className="content-container"
            ),
        ], className="container-width")
    ], className="background"),

    dcc.Location(id='location', refresh=False),
], style={'width': '70%', 'margin': '0 auto'})


toc = html.Div([
        html.Ul([
            html.Li(dcc.Link('Home', href='/')),
            html.Li(dcc.Link('Example', href='/example')),
        ])
    ])

pages = {
    'index': {
        'url': '/',
        'content': html.Div([
            html.H1('Table of Content'),
            toc,
        ], className="toc", style={'textAlign': 'center'})
    },
    'example': {
        'url': '/example',
        'content': layout
    }
}


@app.callback(
    Output('page', 'children'),
    [Input('location', 'pathname')])
def display_content(pathname):
    if pathname is None:
        return html.Div()
    matched = [c for c in list(pages.keys())
               if pages[c]['url'] == pathname]

    if matched and matched[0] == 'example':
        create_callbacks()
        c = pages[matched[0]]['content']
        page_content = c()
        content = html.Div(page_content)
    else:
        content = pages['index']['content']

    return content

if __name__ == '__main__':
    create_callbacks()
    app.run_server(debug=True, threaded=True, host='0.0.0.0', port=5005)