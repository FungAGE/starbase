"""Extract data from Django MySQL for validation and load."""


def extract_from_mysql(mysql_config, tables=None):
    """Read configured tables from Django MySQL into a records structure (e.g. dict of DataFrames or list of dicts).

    Args:
        mysql_config: Dict with connection info (e.g. from Django DATABASES['default']):
            ENGINE, NAME, USER, PASSWORD, HOST, PORT.
        tables: Optional list of table names to extract; default None means all relevant tables.

    Returns:
        records: Structure suitable for validators/transformers (format TBD; e.g. dict[str, pd.DataFrame]).
    """
    # TODO: Use sqlalchemy + pymysql to connect and read tables into pandas DataFrames or similar.
    # from sqlalchemy import create_engine
    # engine = create_engine(...)
    # return {table: pd.read_sql(f"SELECT * FROM {table}", engine) for table in tables or DEFAULT_TABLES}
    return {}
