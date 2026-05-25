#!/usr/bin/env python3
"""
Fill missing starship_features fields from a starfish stats TSV.

Usage:
    python -m src.database.migrations.fill_features_from_stats_tsv \\
        /path/to/mycodb.final.starships.named.stats

    # Preview changes without writing:
    python -m src.database.migrations.fill_features_from_stats_tsv \\
        /path/to/mycodb.final.starships.named.stats --dry-run

The TSV is expected to have (tab-separated) these columns:
    elementID  elementCaptainID  elementContigID  elementBegin  elementEnd
    elementLength  elementStrand  emptySiteID  emptyContigID  emptyBegin
    emptyEnd  emptyLength  emptyStrand  emptySiteSeq  quality  [warnings]

Lines starting with '#' are treated as comments / alternate headers and skipped.

Matching strategy
-----------------
Primary match:   starshipID == elementID  AND  captainID == elementCaptainID
Fallback match:  starshipID == elementID  only (used when captainID is NULL in the DB)

Only NULL / empty fields are updated — existing values are never overwritten.
"""

import argparse
import csv
import sys
import os
from typing import Dict, List, Optional, Tuple

# Allow both `python -m ...` and direct execution
_HERE = os.path.dirname(__file__)
_PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from src.config.database import StarbaseSession
from src.config.logging import get_logger
from src.database.models.schema import StarshipFeatures

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# TSV column → ORM attribute mapping
# ---------------------------------------------------------------------------
TSV_TO_DB: Dict[str, str] = {
    "elementContigID": "contigID",
    "elementBegin": "elementBegin",
    "elementEnd": "elementEnd",
    "elementLength": "elementLength",
    "elementStrand": "strand",
    "emptySiteID": "emptySiteID",
    "emptyContigID": "emptyContig",
    "emptyBegin": "emptyBegin",
    "emptyEnd": "emptyEnd",
    "emptySiteSeq": "emptySeq",
    "quality": "boundaryType",
}

# Fields that should be cast to int before storage
INT_FIELDS = {"elementBegin", "elementEnd", "elementLength", "emptyBegin", "emptyEnd"}

# Sentinel values that count as "missing" in the TSV
MISSING_VALUES = {".", "", "NA", "N/A", "None", "NULL", "null"}


def _is_missing(value: Optional[str]) -> bool:
    return value is None or value.strip() in MISSING_VALUES


def _cast(db_attr: str, raw: str):
    """Convert a raw TSV string to the appropriate Python type."""
    if db_attr in INT_FIELDS:
        try:
            return int(raw)
        except (ValueError, TypeError):
            return None
    return raw.strip() or None


def load_tsv(path: str) -> List[Dict[str, str]]:
    """
    Parse the starfish stats TSV, skipping comment/alternate header lines.

    Returns a list of dicts keyed by the column names from the first non-comment
    header line found.
    """
    rows: List[Dict[str, str]] = []
    header: Optional[List[str]] = None

    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for line in reader:
            if not line:
                continue
            # Lines that start with '#' are comments or duplicate headers
            if line[0].startswith("#"):
                # Use the last such line as the real header (strip leading '#')
                candidate = [c.strip().lstrip("#").strip() for c in line]
                if "elementID" in candidate:
                    header = candidate
                continue
            if header is None:
                # No '#' header found yet — treat first data line as header
                header = [c.strip() for c in line]
                continue
            row = dict(zip(header, [c.strip() for c in line]))
            rows.append(row)

    logger.info(f"Loaded {len(rows)} rows from {path}")
    return rows


def fill_features(tsv_path: str, dry_run: bool = True) -> Dict:
    """
    Fill missing starship_features fields from the TSV.

    Returns a summary dict with counts.
    """
    rows = load_tsv(tsv_path)

    # Build a lookup: (starshipID, captainID) -> row, and (starshipID,) -> [rows]
    by_both: Dict[Tuple[str, str], Dict] = {}
    by_ship: Dict[str, List[Dict]] = {}
    for row in rows:
        sid = row.get("elementID", "").strip()
        cid = row.get("elementCaptainID", "").strip()
        if not sid:
            continue
        by_both[(sid, cid)] = row
        by_ship.setdefault(sid, []).append(row)

    session = StarbaseSession()
    summary = {
        "tsv_rows": len(rows),
        "db_features": 0,
        "matched": 0,
        "unmatched": 0,
        "fields_filled": 0,
        "features_updated": 0,
        "dry_run": dry_run,
    }

    try:
        features: List[StarshipFeatures] = session.query(StarshipFeatures).all()
        summary["db_features"] = len(features)
        logger.info(f"Queried {len(features)} starship_features rows from DB")

        for feat in features:
            sid = feat.starshipID or ""
            cid = feat.captainID or ""

            # Try specific match first, then fall back to starshipID only
            tsv_row = by_both.get((sid, cid))
            if tsv_row is None and sid in by_ship:
                candidates = by_ship[sid]
                tsv_row = candidates[0] if len(candidates) == 1 else None
                if tsv_row is None:
                    logger.debug(
                        f"Ambiguous: {len(candidates)} TSV rows for starshipID={sid!r}; "
                        "skipping (supply captainID in DB to disambiguate)"
                    )

            if tsv_row is None:
                summary["unmatched"] += 1
                continue

            summary["matched"] += 1
            fields_changed = 0

            for tsv_col, db_attr in TSV_TO_DB.items():
                raw = tsv_row.get(tsv_col)
                if _is_missing(raw):
                    continue  # TSV has nothing useful for this field

                current = getattr(feat, db_attr, None)
                # Only fill if the DB value is truly empty
                if current is not None and str(current).strip() not in MISSING_VALUES:
                    continue

                value = _cast(db_attr, raw)
                if value is None:
                    continue

                if not dry_run:
                    setattr(feat, db_attr, value)
                fields_changed += 1
                logger.debug(
                    f"[{'DRY' if dry_run else 'SET'}] "
                    f"starshipID={sid!r} {db_attr}: {current!r} → {value!r}"
                )

            if fields_changed:
                summary["fields_filled"] += fields_changed
                summary["features_updated"] += 1

        if not dry_run:
            session.commit()
            logger.info(
                f"Committed: {summary['features_updated']} features updated, "
                f"{summary['fields_filled']} fields filled"
            )
        else:
            logger.info(
                f"[DRY RUN] Would update {summary['features_updated']} features, "
                f"filling {summary['fields_filled']} fields"
            )

    except Exception as exc:
        logger.error(f"Error filling features: {exc}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()

    return summary


def print_summary(summary: Dict) -> None:
    mode = "DRY RUN" if summary["dry_run"] else "APPLIED"
    print(f"\n{'=' * 60}")
    print(f"fill_features_from_stats_tsv  [{mode}]")
    print(f"{'=' * 60}")
    print(f"  TSV rows loaded:          {summary['tsv_rows']}")
    print(f"  DB features queried:      {summary['db_features']}")
    print(f"  Matched to TSV:           {summary['matched']}")
    print(f"  No TSV match:             {summary['unmatched']}")
    print(f"  Features to update:       {summary['features_updated']}")
    print(f"  Fields to fill:           {summary['fields_filled']}")
    if summary["dry_run"]:
        print("\n  Re-run with --apply to write changes.")
    print(f"{'=' * 60}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Fill missing starship_features fields from a starfish stats TSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "tsv",
        help="Path to the starfish stats TSV (e.g. mycodb.final.starships.named.stats)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        default=True,
        help="Preview changes without writing to the DB (default: True)",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Write changes to the DB (overrides --dry-run)",
    )
    args = parser.parse_args()

    dry_run = not args.apply

    if not os.path.isfile(args.tsv):
        parser.error(f"TSV file not found: {args.tsv}")

    summary = fill_features(tsv_path=args.tsv, dry_run=dry_run)
    print_summary(summary)


if __name__ == "__main__":
    main()
