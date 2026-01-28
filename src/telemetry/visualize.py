import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from src.config.settings import PAGE_MAPPING
from src.config.logging import get_logger

logger = get_logger(__name__)

def create_time_series_figure(time_series_data):
    """Create a time series visualization."""
    if time_series_data:
        df = pd.DataFrame(time_series_data, columns=["date", "count"])
        df["date"] = pd.to_datetime(df["date"])

        # Calculate moving average for trend
        df["moving_avg"] = df["count"].rolling(window=3, center=True).mean()

        fig = go.Figure()

        # Add main line with gradient
        fig.add_trace(
            go.Scatter(
                x=df["date"],
                y=df["count"],
                mode="lines+markers",
                name="Daily Visitors",
                line=dict(color="#6366f1", width=3, shape="spline"),
                marker=dict(size=8, color="#6366f1", line=dict(width=2, color="white")),
                fill="tonexty",
                fillcolor="rgba(99, 102, 241, 0.1)",
            )
        )

        # Add trend line
        fig.add_trace(
            go.Scatter(
                x=df["date"],
                y=df["moving_avg"],
                mode="lines",
                name="Trend (3-day avg)",
                line=dict(color="#ef4444", width=2, dash="dash"),
                opacity=0.7,
            )
        )

        # Calculate growth rate
        if len(df) > 1:
            current = df["count"].iloc[-1]
            previous = df["count"].iloc[-2]
            growth_rate = ((current - previous) / previous * 100) if previous > 0 else 0
            growth_text = (
                f"📈 +{growth_rate:.1f}%"
                if growth_rate > 0
                else f"📉 {growth_rate:.1f}%"
            )
        else:
            growth_text = "📊 New data"

        fig.update_layout(
            title=dict(
                text=f"Visitor Trends {growth_text}",
                x=0.5,
                font=dict(size=18, color="#1f2937"),
            ),
            xaxis=dict(
                title="Date",
                showgrid=True,
                gridcolor="rgba(0,0,0,0.05)",
                zeroline=False,
                tickformat="%b %d",
            ),
            yaxis=dict(
                title="Unique Visitors",
                showgrid=True,
                gridcolor="rgba(0,0,0,0.05)",
                zeroline=False,
            ),
            plot_bgcolor="white",
            paper_bgcolor="white",
            height=400,
            margin=dict(l=50, r=20, t=60, b=50),
            hovermode="x unified",
            hoverlabel=dict(bgcolor="white", font_size=12, font_family="Inter"),
            showlegend=True,
            legend=dict(
                orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1
            ),
        )

        # Focus on last 14 days
        if not df.empty:
            end_date = df["date"].max()
            start_date = end_date - pd.Timedelta(days=14)
            fig.update_layout(xaxis_range=[start_date, end_date])

    else:
        fig = go.Figure()
        fig.update_layout(
            title="No visitor data available",
            xaxis_title="Date",
            yaxis_title="Unique Visitors",
            plot_bgcolor="white",
            paper_bgcolor="white",
            height=400,
        )

    return fig


def create_endpoints_figure(endpoints_data):
    """Create a modern horizontal bar chart for page visits."""
    if endpoints_data:
        filtered_data = [
            (PAGE_MAPPING.get(row[0], row[0]), row[1])
            for row in endpoints_data
            if row[0] in PAGE_MAPPING
        ]

        if filtered_data:
            endpoints = [row[0] for row in filtered_data]
            counts = [row[1] for row in filtered_data]

            # Create color gradient
            colors = px.colors.qualitative.Set3[: len(endpoints)]

            fig = go.Figure(
                go.Bar(
                    y=endpoints,
                    x=counts,
                    orientation="h",
                    marker=dict(color=colors, line=dict(color="white", width=1)),
                    text=counts,
                    textposition="auto",
                    texttemplate="%{text}",
                    textfont=dict(color="white", size=12),
                )
            )

            fig.update_layout(
                title=dict(
                    text="Page Popularity", x=0.5, font=dict(size=18, color="#1f2937")
                ),
                xaxis=dict(
                    title="Unique Visitors",
                    showgrid=True,
                    gridcolor="rgba(0,0,0,0.05)",
                    zeroline=False,
                ),
                yaxis=dict(showgrid=False, zeroline=False),
                plot_bgcolor="white",
                paper_bgcolor="white",
                height=400,
                margin=dict(l=50, r=20, t=60, b=50),
                bargap=0.3,
                showlegend=False,
            )

        else:
            fig = go.Figure()
            fig.update_layout(
                title="No page visit data available",
                plot_bgcolor="white",
                paper_bgcolor="white",
                height=400,
            )
    else:
        fig = go.Figure()
        fig.update_layout(
            title="No page visit data available",
            plot_bgcolor="white",
            paper_bgcolor="white",
            height=400,
        )

    return fig


def create_map_figure(locations):
    """Create a modern map visualization."""
    fig = go.Figure()

    fig.update_layout(
        mapbox=dict(
            style="carto-positron",
            zoom=1,
            center=dict(lat=20, lon=0),
        ),
        margin={"r": 0, "t": 30, "l": 0, "b": 0},
        height=400,
        title=dict(
            text="Global Visitor Distribution",
            x=0.5,
            font=dict(size=18, color="#1f2937"),
        ),
        paper_bgcolor="white",
    )

    if locations and len(locations) > 0:
        df = pd.DataFrame(locations)

        # Count visitors per location
        location_counts = (
            df.groupby(["lat", "lon", "city", "country"])
            .size()
            .reset_index(name="visits")
        )

        # Create size gradient
        min_visits = location_counts["visits"].min()
        max_visits = location_counts["visits"].max()
        location_counts["marker_size"] = 8 + (
            location_counts["visits"] - min_visits
        ) * (25 / (max_visits - min_visits if max_visits > min_visits else 1))

        fig.add_trace(
            go.Scattermapbox(
                lat=location_counts["lat"],
                lon=location_counts["lon"],
                mode="markers",
                marker=dict(
                    size=location_counts["marker_size"],
                    color=location_counts["visits"],
                    colorscale="Viridis",
                    opacity=0.8,
                    sizemode="diameter",
                    colorbar=dict(title="Visits", x=0.95, len=0.8),
                ),
                text=location_counts.apply(
                    lambda row: f"<b>{row['city']}, {row['country']}</b><br>👥 {int(row['visits'])} visitors",
                    axis=1,
                ),
                hoverinfo="text",
                hovertemplate="%{text}<extra></extra>",
            )
        )

        # Adjust center and zoom
        center_lat = location_counts["lat"].mean()
        center_lon = location_counts["lon"].mean()
        fig.update_layout(
            mapbox=dict(center=dict(lat=center_lat, lon=center_lon), zoom=1.5)
        )

    return fig