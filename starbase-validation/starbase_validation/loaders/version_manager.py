"""Manage versioned SQLite database creation and metadata."""


class DatabaseVersionManager:
    """Create versioned SQLite databases with metadata (version, quality score, etc.)."""

    def __init__(self, version_registry_path=None):
        """Optional version_registry_path: file or DB storing last version for increment."""
        self.version_registry_path = version_registry_path

    def get_last_version(self):
        """Return last published version string (e.g. '1.2.3') or None."""
        # TODO: Read from registry or metadata table in latest DB.
        return None

    def calculate_version(self, records, last_version=None):
        """Compute next semantic version (MAJOR.MINOR.PATCH) from last_version and change heuristics."""
        last = last_version or self.get_last_version() or "0.0.0"
        # TODO: has_breaking_changes -> major; has_new_records -> minor; else patch.
        major, minor, patch = (int(x) for x in last.split("."))
        return f"{major}.{minor}.{patch + 1}"

    def create_version(self, records, output_dir, validation_results=None):
        """Create new versioned DB and metadata table. Returns version string."""
        version = self.calculate_version(records)
        # TODO: Call sqlite_loader, then add database_version table with version, created_at,
        #       validation_passed, total_records, quality_score, warnings, source.
        return version
