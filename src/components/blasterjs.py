import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc

# Initialize the Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div([
    html.H1("BLAST Results Viewer"),
    
    # Load the HTML template via iframe
    html.Iframe(
        src="/assets/html/blast.html",  # path to your HTML file
        style={"width": "100%", "height": "800px", "border": "none"}
    )
])

# Run the app
if __name__ == "__main__":
    app.run_server(debug=True, threaded=True, host="localhost", port=8005)
