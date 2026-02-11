"""Load validated records into a SQLite database."""


def load_to_sqlite(records, output_path, schema=None):
    """Write records to SQLite at output_path. Optionally apply schema (tables DDL).

    Args:
        records: Structure produced by extract/transform (e.g. dict[str, pd.DataFrame]).
        output_path: Path for the .db file.
        schema: Optional DDL or schema module for table creation.

    Returns:
        output_path (str).
    """
    # TODO: Create SQLite DB, create tables from schema, insert from records.
    # from sqlalchemy import create_engine
    # engine = create_engine(f"sqlite:///{output_path}")
    # for table, df in records.items(): df.to_sql(table, engine, if_exists="replace", index=False)
    return output_path

