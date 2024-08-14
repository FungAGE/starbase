import plotly.express as px


def agg_df(df, groups):
    agg = df.groupby(groups).starshipID.agg(
        count="count", nunique="nunique", duplicates=lambda x: x.size - x.nunique()
    )

    agg = agg.reset_index()

    return agg


def create_sunburst_plot(df, path, title):
    selection = agg_df(df, path)

    pie = px.sunburst(
        selection,
        path=path,
        values="count",
    )

    pie.update_layout(
        autosize=True,
        title_font=dict(size=24),
        title={
            "text": title,
            "y": 1,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        margin=dict(t=50, l=0, r=0, b=0),
    )
    return pie
