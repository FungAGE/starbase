"""
Update joined_ships for source = 'starbase_fillin' using a fill-in metadata TSV.

Matches rows by joined_ships.starshipID = TSV starshipID, then:
- Updates joined_ships: genome_id (from genomes.ome = TSV 'ome'), tax_id, optional ship_family_id/ship_navis_id.
- Creates new starship_features rows using TSV contigID, element_start (-> elementBegin), element_end (-> elementEnd), element_length (-> elementLength).

Usage:
    python -m src.database.migrations.update_joined_ships_from_fillin_metadata /path/to/metadata.tsv
    python -m src.database.migrations.update_joined_ships_from_fillin_metadata /path/to/metadata.tsv --dry-run
"""

import argparse
import csv
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

from sqlalchemy import text
from src.config.logging import get_logger
from src.database.sql_engine import get_starbase_session

logger = get_logger(__name__)

FILLIN_SOURCE = "starbase_fillin"


def load_tsv(path: Path):
    """Load metadata TSV; return list of dicts with keys: starshipID, ome, ship_family, ship_navis, etc."""
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def run(path: Path, dry_run: bool = False):
    rows = load_tsv(path)
    if not rows:
        logger.warning("TSV is empty")
        return

    with get_starbase_session() as session:
        # Resolve genome_id and tax_id by ome
        ome_to_genome = {}
        genome_id_to_tax_id = {}
        result = session.execute(
            text("SELECT id, ome, taxonomy_id FROM genomes WHERE ome IS NOT NULL AND trim(ome) != ''")
        )
        for r in result:
            ome_to_genome[r.ome.strip()] = r.id
            genome_id_to_tax_id[r.id] = r.taxonomy_id

        # Optional: resolve family and navis by name
        name_to_family_id = {}
        fam_result = session.execute(text("SELECT id, familyName FROM family_names"))
        for r in fam_result:
            if r.familyName:
                name_to_family_id[r.familyName.strip()] = r.id

        name_to_navis_id = {}
        nav_result = session.execute(text("SELECT id, navis_name FROM navis_names"))
        for r in nav_result:
            if r.navis_name:
                name_to_navis_id[r.navis_name.strip()] = r.id

        updated = 0
        features_created = 0
        skipped_no_ome = 0
        skipped_no_genome = 0
        skipped_no_joined = 0

        for row in rows:
            starship_id = (row.get("starshipID") or "").strip()
            ome = (row.get("ome") or "").strip()
            if not starship_id:
                continue
            if not ome:
                skipped_no_ome += 1
                continue

            genome_id = ome_to_genome.get(ome)
            if genome_id is None:
                skipped_no_genome += 1
                logger.debug("No genome for ome=%r", ome)
                continue

            tax_id = genome_id_to_tax_id.get(genome_id)

            # Find joined_ships row: starshipID + source = starbase_fillin (need ship_id, accession_id for starship_features)
            j = session.execute(
                text(
                    "SELECT id, genome_id, tax_id, ship_family_id, ship_navis_id, ship_id, accession_id FROM joined_ships "
                    "WHERE starshipID = :sid AND source = :src"
                ),
                {"sid": starship_id, "src": FILLIN_SOURCE},
            ).fetchone()

            if not j:
                skipped_no_joined += 1
                continue

            if j.ship_id is None:
                skipped_no_joined += 1
                logger.debug("joined_ships id=%s has no ship_id, skip", j.id)
                continue

            # Build update: at least genome_id (and tax_id when available)
            set_parts = ["genome_id = :gid"]
            params = {"gid": genome_id, "jid": j.id}

            if tax_id is not None:
                set_parts.append("tax_id = :tid")
                params["tid"] = tax_id

            ship_family = (row.get("ship_family") or "").strip()
            ship_navis = (row.get("ship_navis") or "").strip()
            if ship_family and ship_family in name_to_family_id:
                set_parts.append("ship_family_id = :fid")
                params["fid"] = name_to_family_id[ship_family]
            if ship_navis and ship_navis in name_to_navis_id:
                set_parts.append("ship_navis_id = :nid")
                params["nid"] = name_to_navis_id[ship_navis]

            if dry_run:
                logger.info(
                    "Would update joined_ships id=%s starshipID=%s -> genome_id=%s tax_id=%s",
                    j.id,
                    starship_id,
                    genome_id,
                    tax_id,
                )
            else:
                session.execute(
                    text(
                        "UPDATE joined_ships SET " + ", ".join(set_parts) + " WHERE id = :jid"
                    ),
                    params,
                )
            updated += 1

            # Create new starship_features row from TSV contigID, element_start, element_end, element_length
            contig_id = (row.get("contigID") or "").strip() or None
            try:
                elem_start = int(row["element_start"]) if row.get("element_start") else None
            except (ValueError, TypeError):
                elem_start = None
            try:
                elem_end = int(row["element_end"]) if row.get("element_end") else None
            except (ValueError, TypeError):
                elem_end = None
            try:
                elem_len = int(row["element_length"]) if row.get("element_length") else None
            except (ValueError, TypeError):
                elem_len = None

            has_feature_data = (
                contig_id is not None or elem_start is not None or elem_end is not None or elem_len is not None
            )
            if not has_feature_data:
                continue

            if dry_run:
                logger.info(
                    "Would insert starship_features ship_id=%s accession_id=%s contigID=%s elementBegin=%s elementEnd=%s elementLength=%s",
                    j.ship_id,
                    j.accession_id,
                    contig_id,
                    elem_start,
                    elem_end,
                    elem_len,
                )
                features_created += 1
            else:
                session.execute(
                    text("""
                        INSERT INTO starship_features
                        (accession_id, ship_id, contigID, starshipID, elementBegin, elementEnd, elementLength)
                        VALUES (:aid, :sid, :contig, :starship_id, :elem_begin, :elem_end, :elem_len)
                    """),
                    {
                        "aid": j.accession_id,
                        "sid": j.ship_id,
                        "contig": contig_id,
                        "starship_id": starship_id,
                        "elem_begin": elem_start,
                        "elem_end": elem_end,
                        "elem_len": elem_len,
                    },
                )
                features_created += 1

        logger.info(
            "Updated %s joined_ships rows, created %s starship_features rows (dry_run=%s). Skipped: no_ome=%s, no_genome=%s, no_joined_row=%s",
            updated,
            features_created,
            dry_run,
            skipped_no_ome,
            skipped_no_genome,
            skipped_no_joined,
        )


def main():
    parser = argparse.ArgumentParser(
        description="Update joined_ships (source=starbase_fillin) from fill-in metadata TSV"
    )
    parser.add_argument(
        "metadata_tsv",
        type=Path,
        help="Path to metadata TSV (columns: starshipID, ome, ship_family, ship_navis, ...)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Log updates without writing to the database",
    )
    args = parser.parse_args()

    if not args.metadata_tsv.exists():
        logger.error("File not found: %s", args.metadata_tsv)
        sys.exit(1)

    run(args.metadata_tsv, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
