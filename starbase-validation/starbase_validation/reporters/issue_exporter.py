"""Export validation issues back to Django for curator review."""


def report_to_django(issues, mysql_config, options=None):
    """Write validation issues to Django MySQL (e.g. validation_issues table) or CSV for manual import.

    Args:
        issues: List of issue dicts (type, category, table_name, record_id, details, etc.).
        mysql_config: Django-style DB config for MySQL.
        options: Optional dict with export_format ('table' | 'csv'), csv_path, etc.
    """
    # TODO: Create engine, write to validation_issues table or issues_df.to_csv(csv_path).
    # TODO: Optional send_email_to_curators(...).
    pass

