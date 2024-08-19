from dash import dash_table, html


def make_table(df, columns=None):
    if columns is None:
        columns = df.columns

    table = html.Div(
        [
            dash_table.DataTable(
                id="table",
                columns=[
                    {
                        "name": i,
                        "id": i,
                        "deletable": False,
                        "selectable": True,
                    }
                    for i in columns
                ],
                data=df.to_dict("records"),
                editable=False,
                filter_action="native",
                sort_action="native",
                sort_mode="multi",
                row_selectable="multi",
                row_deletable=False,
                selected_columns=[],
                selected_rows=[],
                page_action="native",
                page_current=0,
                page_size=10,
                style_table={
                    "overflowX": "auto",
                    "maxWidth": "95%",
                },
                style_cell={
                    "minWidth": "150px",
                    "width": "150px",
                    "maxWidth": "150px",
                    "whiteSpace": "normal",
                },
            ),
            html.Div(id="table-container"),
        ]
    )
    return table
