import dash
from dash import dcc, html
import pandas as pd
from dash.dependencies import Input, Output
from sqlalchemy import create_engine

# Create a SQLite engine and initialize the table
if __name__ == '__main__':
    engine = create_engine('sqlite:///database_folder/starbase.sqlite')
    query = "SELECT name FROM sqlite_master WHERE type='table'"
    sql_tbls = pd.read_sql_query(query, engine)

# def load_ship_metadata(app):
#     @app.callback(Output("joined-ships", "data"), Input("url", "href"))
def ship_metadata():
    query = """
    SELECT j.*, t."order", t.family, f.longFamilyID, f.familyName, a.accession_tag
    FROM joined_ships j
    JOIN taxonomy t ON j.taxid = t.id
    JOIN family_names f ON j.ship_family_id = f.id
    JOIN accessions a ON j.ship_id = a.id
    """
    df = pd.read_sql_query(query, engine)
    # for saving in cache:
    # data = df.to_dict(orient="records")
    return df


# def load_ship_papers(app):
#     @app.callback(
#         Output("paper-cache", "data"),
#         Input("url", "href"),
#     )

def ship_papers():
    query = """
    SELECT p.Title, p.Author, p.PublicationYear, p.DOI, p.Url, p.shortCitation, f.familyName, f.type_element_reference
    FROM papers p
    JOIN family_names f ON p.shortCitation = f.type_element_reference
    """
    df = pd.read_sql_query(query, engine)
    return df

# Callback to handle FASTA download
# def dl_fa(app):
#     @app.callback(
#         Output("download-fasta", "data"),
#         Input("download-all-btn", "n_clicks"),
#         Input("download-selected-btn", "n_clicks"),
#         [
#             State("download-table", "derived_virtual_data"),
#             State("download-table", "derived_virtual_selected_rows"),
#         ],
#         prevent_initial_call=True,
#     )
def fasta_table(ship_names=None):
    query = f"SELECT ship_name, ship_sequence FROM ships"
    if ship_names is not None and len(ship_names) > 1:
        placeholders = ",".join(["?"] * len(ship_names))
        query += f"WHERE ship_name IN ({placeholders})"
        df = pd.read_sql_query(query, engine, params=ship_names)
    else:
        df = pd.read_sql_query(query, engine)

    return df

    # # Create FASTA content
    # fasta_content = [
    #     f">{row['ship_name']}\n{row['ship_sequence']}"
    #     for _, row in df.iterrows()
    # ]
    # fasta_str = "\n".join(fasta_content)

    # # Send the FASTA file for download
    # return dcc.send_string(fasta_str, filename="starships.fasta")