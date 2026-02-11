"""Generate quality and validation reports."""


def generate_report(validation_results, version=None, output_path=None):
    """Write a human-readable or machine report for validation run and optional DB version.

    Args:
        validation_results: ValidationResults from validate_all.
        version: Optional version string of created database.
        output_path: Optional file path to write report (e.g. JSON or Markdown).
    """
    # TODO: Summarize blocking_issues, warnings, score; optionally write to file.
    pass

