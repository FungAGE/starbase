"""Validation results container and threshold evaluation."""


class ValidationResults:
    """Aggregated results from all validators with blocking/warning classification."""

    def __init__(self):
        self.issues = []  # list of dicts: {type, category, table, record_id, details, severity}
        self.blocking_issues = []
        self.warnings = []
        self.score = 1.0  # 0.0–1.0 quality score

    def add(self, other):
        """Merge another ValidationResults or list of issues into this one."""
        if isinstance(other, ValidationResults):
            self.issues.extend(other.issues)
            self.blocking_issues.extend(other.blocking_issues)
            self.warnings.extend(other.warnings)
        elif isinstance(other, (list, tuple)):
            self.issues.extend(other)
        return self

    def evaluate_thresholds(self, quality_rules):
        """Classify issues as blocking vs warning using quality_rules and compute score."""
        blocking_types = set(quality_rules.get("blocking_issues", []))
        warning_types = set(quality_rules.get("warning_issues", []))
        self.blocking_issues = [i for i in self.issues if i.get("type") in blocking_types]
        self.warnings = [i for i in self.issues if i.get("type") in warning_types]
        # Simple score: 1.0 if no blocking, reduced by blocking count
        if self.blocking_issues:
            self.score = max(0.0, 1.0 - 0.2 * len(self.blocking_issues))
        else:
            self.score = 1.0 - 0.02 * min(len(self.warnings), 25)
        return self

    def has_blocking_issues(self):
        return len(self.blocking_issues) > 0

    def is_acceptable(self):
        return not self.has_blocking_issues()
