from src.database.sql_manager import (
    fetch_meta_data,
    fetch_paper_data,
    fetch_ship_table,
    fetch_all_ships,
    get_database_stats,
    fetch_download_data
)

precompute_tasks = {
    "meta_data": lambda: fetch_meta_data(curated=True),
    "paper_data": fetch_paper_data,
    "ship_table": fetch_ship_table,
    "all_ships": fetch_all_ships,
    "database_stats": get_database_stats,
    "download_data_curated_true_derep_false": lambda: fetch_download_data(curated=True, dereplicate=False),
    "download_data_curated_true_derep_true": lambda: fetch_download_data(curated=True, dereplicate=True),
    "download_data_curated_false_derep_false": lambda: fetch_download_data(curated=False, dereplicate=False),
    "download_data_curated_false_derep_true": lambda: fetch_download_data(curated=False, dereplicate=True),
}