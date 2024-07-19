import dash
from dash import html

dash.register_page(__name__)

# Define the layout
layout = html.Div(
    [html.Iframe(src="tmp/mummer/result.html", width="100%", height="500")]
)

# Read the HTML file
with open("tmp/mummer/result.html", "r") as file:
    html_content = file.read()

# Define the layout
layout = html.Div([html.Div(html_content, style={"whiteSpace": "pre-line"})])
