import logging
from src.components.sql_queries import (
    fetch_meta_data,
    fetch_paper_data,
    cache_sunburst_plot,
    fetch_download_data,
    fetch_all_ships,
    fetch_ship_table,
    fetch_all_captains,
    fetch_captain_tree,
    fetch_sf_data,
)

logger = logging.getLogger(__name__)

families = [
    "Phoenix",
    "Hephaestus",
    "Tardis",
    "Serenity",
    "Prometheus",
    "Enterprise",
    "Galactica",
    "Moya",
    "Arwing",
    "Voyager",
    "Family-11",
]


def precompute_all():
    """Precompute and cache all necessary data and figures."""

    logger.info("Starting precomputation of all data and figures.")

    try:
        logger.info("Precomputing meta data...")
        meta_data = fetch_meta_data(curated=True)
        logger.info("Meta data precomputed and cached successfully.")
    except Exception as e:
        logger.error(f"Failed to precompute meta data: {str(e)}")

    try:
        logger.info("Precomputing paper data...")
        paper_data = fetch_paper_data()
        logger.info("Paper data precomputed and cached successfully.")
    except Exception as e:
        logger.error(f"Failed to precompute paper data: {str(e)}")

    try:
        logger.info("Creating sunburst figures...")
        for family in families:
            sunburst_figure = cache_sunburst_plot(family, meta_data)
        logger.info("Sunburst figures created and cached successfully.")
    except Exception as e:
        logger.error(f"Failed to create sunburst figures: {str(e)}")

    try:
        logger.info("Precomputing fetch_download_data...")
        download_data = fetch_download_data()
    except Exception as e:
        logger.error(f"Failed to fetch_download_data: {str(e)} ")
    try:
        logger.info("Precomputing fetch_all_ships...")
        all_ships = fetch_all_ships()
    except Exception as e:
        logger.error(f"Failed to fetch_all_ships: {str(e)} ")
    try:
        logger.info("Precomputing fetch_ship_table...")
        ship_table = fetch_ship_table()
    except Exception as e:
        logger.error(f"Failed to fetch_ship_table: {str(e)} ")
    try:
        logger.info("Precomputing fetch_all_captains...")
        all_captains = fetch_all_captains()
    except Exception as e:
        logger.error(f"Failed to fetch_all_captains: {str(e)} ")
    try:
        logger.info("Precomputing fetch_captain_tree...")
        captain_tree = fetch_captain_tree()
    except Exception as e:
        logger.error(f"Failed to fetch_captain_tree: {str(e)} ")
    try:
        logger.info("Precomputing fetch_sf_data...")
        sf_data = fetch_sf_data()
    except Exception as e:
        logger.error(f"Failed to fetch_sf_data: {str(e)} ")

    logger.info("Precomputation complete.")
