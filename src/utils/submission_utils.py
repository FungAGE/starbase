#!/usr/bin/env python3
"""
Comprehensive submission utility for adding new Starship sequences to the database.

This module handles:
1. Input validation and schema enforcement
2. Sequence uniqueness checking (MD5-based)
3. Handling duplicate sequences from different genomes/taxa
4. Database insertion with proper foreign key relationships
5. Accession assignment using the standard workflow

Usage:
    from src.utils.submission_utils import SubmissionProcessor

    processor = SubmissionProcessor()

    # Single submission
    result = processor.process_submission({
        'sequence': 'ATGC...',
        'starshipID': 'SS-1.1',
        'genus': 'Fusarium',
        'species': 'oxysporum',
        'evidence': 'starfish',
        'source': 'publication_2024',
        # ... other fields
    })

    # Batch submissions
    results = processor.process_batch([submission1, submission2, ...])
"""

import pandas as pd
from typing import Dict, List, Optional, Tuple, Any
from datetime import datetime

from src.config.logging import get_logger
from src.config.database import StarbaseSession
from src.database.models.schema import (
    Ships,
    Accessions,
    ShipAccessions,
    JoinedShips,
    Taxonomy,
    FamilyNames,
    Navis,
    Haplotype,
    Genome,
    Gff,
    StarshipFeatures,
)
from src.utils.classification_utils import (
    assign_accession,
    ensure_ship_has_ssb,
    generate_md5_hash,
)
from src.utils.seq_utils import clean_sequence, revcomp
from src.database.sql_manager import fetch_ships

logger = get_logger(__name__)


# ============================================================================
# SCHEMA DEFINITION
# ============================================================================

REQUIRED_FIELDS = {
    "sequence": str,
    "starshipID": str,
    "evidence": str,
    "source": str,
}

OPTIONAL_FIELDS = {
    # Taxonomy
    "genus": str,
    "species": str,
    "strain": str,
    "tax_name": str,  # Full taxonomy name
    "tax_order": str,
    "tax_family": str,
    # Classification
    "ship_family": str,
    "ship_navis": str,
    "ship_haplotype": str,
    # Genome info
    "assembly_accession": str,
    "ome": str,
    "genome_contig": str,
    "contig_length": int,
    "contig_id": str,
    # Location
    "element_start": int,
    "element_end": int,
    "element_strand": str,
    # Captain
    "captain_id": int,
    # Curation
    "curated_status": str,  # 'curated', 'uncurated', 'needs_review'
    "curator": str,
    "notes": str,
    # GFF annotations (list of dicts)
    "gff_entries": list,
}


# ============================================================================
# VALIDATION
# ============================================================================


class ValidationError(Exception):
    """Custom exception for validation errors."""

    pass


def validate_submission(data: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Validate submission data against schema.

    Args:
        data: Submission data dictionary

    Returns:
        Tuple of (is_valid, list_of_errors)
    """
    errors = []

    # Check required fields
    for field, expected_type in REQUIRED_FIELDS.items():
        if field not in data or data[field] is None:
            errors.append(f"Missing required field: {field}")
        elif not isinstance(data[field], expected_type):
            errors.append(f"Field '{field}' must be of type {expected_type.__name__}")

    # Validate optional fields if present
    for field, expected_type in OPTIONAL_FIELDS.items():
        if field in data and data[field] is not None:
            if not isinstance(data[field], expected_type):
                errors.append(
                    f"Field '{field}' must be of type {expected_type.__name__}"
                )

    # Validate sequence
    if "sequence" in data and data["sequence"]:
        seq = data["sequence"]
        if not seq.strip():
            errors.append("Sequence cannot be empty")
        elif not all(c in "ATGCNatgcn" for c in seq.replace("\n", "").replace(" ", "")):
            errors.append("Sequence contains invalid characters (only ATGCN allowed)")

        # Check minimum length
        if len(seq.replace("\n", "").replace(" ", "")) < 100:
            errors.append("Sequence too short (minimum 100 bp)")

    # Validate curated_status
    if "curated_status" in data and data["curated_status"]:
        valid_statuses = ["curated", "uncurated", "needs_review", "pending"]
        if data["curated_status"] not in valid_statuses:
            errors.append(f"curated_status must be one of: {valid_statuses}")

    # Validate strand
    if "element_strand" in data and data["element_strand"]:
        if data["element_strand"] not in ["+", "-", "."]:
            errors.append("element_strand must be '+', '-', or '.'")

    # Validate and normalize coordinates
    # Note: start can be > end for reverse strand features
    # We'll auto-detect strand and normalize coordinates
    if "element_start" in data and "element_end" in data:
        start = data.get("element_start")
        end = data.get("element_end")
        if start is not None and end is not None:
            # Auto-detect strand if not provided
            if start > end:
                # Reverse strand: swap coordinates and set strand to '-'
                if not data.get("element_strand"):
                    data["element_strand"] = "-"
                # Swap so start < end for storage
                data["element_start"], data["element_end"] = end, start
            elif start < end:
                # Forward strand: set default strand to '+'
                if not data.get("element_strand"):
                    data["element_strand"] = "+"
            else:
                # start == end is invalid
                errors.append("element_start cannot equal element_end")

    return len(errors) == 0, errors


# ============================================================================
# DUPLICATE DETECTION
# ============================================================================


class DuplicateInfo:
    """Information about duplicate sequences."""

    def __init__(
        self,
        is_duplicate: bool,
        existing_ship_id: Optional[int] = None,
        existing_accession: Optional[str] = None,
        different_taxon: bool = False,
        match_type: str = None,
    ):
        self.is_duplicate = is_duplicate
        self.existing_ship_id = existing_ship_id
        self.existing_accession = existing_accession
        self.different_taxon = different_taxon
        self.match_type = match_type  # 'exact', 'similar', 'contained'


def check_sequence_duplicate(
    sequence: str, genus: Optional[str] = None, species: Optional[str] = None
) -> DuplicateInfo:
    """
    Check if sequence already exists in database.

    Returns:
        DuplicateInfo object with details about any duplicates found
    """
    session = StarbaseSession()

    try:
        # Clean and hash the sequence
        clean_seq = clean_sequence(sequence)
        md5_hash = generate_md5_hash(clean_seq)
        revcomp_seq = revcomp(clean_seq)
        md5_revcomp = generate_md5_hash(revcomp_seq)

        # Check for exact MD5 match
        matching_ships = (
            session.query(Ships)
            .filter(
                (Ships.md5 == md5_hash)
                | (Ships.md5 == md5_revcomp)
                | (Ships.rev_comp_md5 == md5_hash)
                | (Ships.rev_comp_md5 == md5_revcomp)
            )
            .all()
        )

        if not matching_ships:
            return DuplicateInfo(is_duplicate=False)

        # Found exact match - check if it's from a different taxon
        for ship in matching_ships:
            # Get taxonomy info for this ship
            joined_ship = (
                session.query(JoinedShips)
                .filter(JoinedShips.ship_id == ship.id)
                .first()
            )

            if joined_ship and joined_ship.taxonomy:
                existing_genus = joined_ship.taxonomy.genus
                existing_species = joined_ship.taxonomy.species

                # Check if taxonomy differs
                different_taxon = (
                    genus
                    and species
                    and (existing_genus != genus or existing_species != species)
                )

                # Get accession
                accession_tag = None
                if ship.accession_id:
                    accession = (
                        session.query(Accessions)
                        .filter(Accessions.id == ship.accession_id)
                        .first()
                    )
                    if accession:
                        accession_tag = accession.accession_tag

                return DuplicateInfo(
                    is_duplicate=True,
                    existing_ship_id=ship.id,
                    existing_accession=accession_tag,
                    different_taxon=different_taxon,
                    match_type="exact",
                )

        # Found match but no taxonomy info
        return DuplicateInfo(
            is_duplicate=True,
            existing_ship_id=matching_ships[0].id,
            existing_accession=None,
            different_taxon=bool(genus and species),
            match_type="exact",
        )

    finally:
        session.close()


# ============================================================================
# DATABASE INSERTION
# ============================================================================


class SubmissionProcessor:
    """Main class for processing submissions."""

    def __init__(self, dry_run: bool = False):
        self.dry_run = dry_run
        self.stats = {
            "total": 0,
            "validated": 0,
            "inserted": 0,
            "duplicates_same_taxon": 0,
            "duplicates_diff_taxon": 0,
            "errors": [],
        }

    def process_submission(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process a single submission.

        Args:
            data: Submission data dictionary

        Returns:
            Result dictionary with status and details
        """
        self.stats["total"] += 1
        result = {
            "success": False,
            "ship_id": None,
            "accession": None,
            "errors": [],
            "warnings": [],
            "duplicate_info": None,
        }

        try:
            # Step 1: Validate
            is_valid, validation_errors = validate_submission(data)
            if not is_valid:
                result["errors"] = validation_errors
                self.stats["errors"].append(
                    {
                        "starshipID": data.get("starshipID", "unknown"),
                        "errors": validation_errors,
                    }
                )
                return result

            self.stats["validated"] += 1

            # Step 2: Check for duplicates
            duplicate_info = check_sequence_duplicate(
                data["sequence"], data.get("genus"), data.get("species")
            )
            result["duplicate_info"] = duplicate_info

            if duplicate_info.is_duplicate and not duplicate_info.different_taxon:
                # Exact duplicate in same taxon - skip
                result["warnings"].append(
                    f"Duplicate sequence found (ship_id={duplicate_info.existing_ship_id}, "
                    f"accession={duplicate_info.existing_accession})"
                )
                self.stats["duplicates_same_taxon"] += 1
                logger.warning(f"Skipping duplicate: {data.get('starshipID')}")
                return result

            if duplicate_info.is_duplicate and duplicate_info.different_taxon:
                result["warnings"].append(
                    "Same sequence in different taxon - will create new joined_ships entry"
                )
                self.stats["duplicates_diff_taxon"] += 1

            # Step 3: Process insertion
            if not self.dry_run:
                ship_id, accession = self._insert_to_database(data, duplicate_info)
                result["ship_id"] = ship_id
                result["accession"] = accession
                result["success"] = True
                self.stats["inserted"] += 1
            else:
                result["success"] = True
                result["warnings"].append("DRY RUN: No actual database changes made")

        except Exception as e:
            logger.error(f"Error processing submission: {str(e)}")
            result["errors"].append(str(e))
            self.stats["errors"].append(
                {"starshipID": data.get("starshipID", "unknown"), "error": str(e)}
            )

        return result

    def _insert_to_database(
        self, data: Dict[str, Any], duplicate_info: DuplicateInfo
    ) -> Tuple[int, str]:
        """
        Insert submission data into database with proper foreign key relationships.

        1. Determine if we need to create a new ship or reuse existing
        2. Assign accession
        3. Create Accessions record
        4. Create Ships record
        5. Ensure this ship has an SSB (ship) accession for display/search
        6. Get or create foreign key records
        7. Create JoinedShips record
        8. Add GFF entries if provided
        9. Create StarshipFeatures if location data provided
        10. Commit all changes

        Returns:
            Tuple of (ship_id, accession_tag)
        """
        session = StarbaseSession()

        try:
            clean_seq = clean_sequence(data["sequence"])

            if duplicate_info.is_duplicate and duplicate_info.different_taxon:
                ship_id = duplicate_info.existing_ship_id
                ship = session.query(Ships).filter(Ships.id == ship_id).first()

                if ship.accession_id:
                    accession = (
                        session.query(Accessions)
                        .filter(Accessions.id == ship.accession_id)
                        .first()
                    )
                    accession_tag = accession.accession_tag if accession else None
                else:
                    existing_ships = fetch_ships(with_sequence=True)
                    accession_tag, _ = assign_accession(clean_seq, existing_ships)

                    new_accession = Accessions(
                        accession_tag=accession_tag, version_tag=1
                    )
                    session.add(new_accession)
                    session.flush()

                    ship.accession_id = new_accession.id
                    accession_tag = new_accession.accession_tag
            else:
                # Create new ship
                existing_ships = fetch_ships(with_sequence=True)
                accession_tag, needs_review = assign_accession(
                    clean_seq, existing_ships
                )

                ship_accession = ShipAccessions(ship_accession_tag=accession_tag)
                session.add(ship_accession)
                session.flush()

                accession = Accessions(accession_tag=accession_tag, version_tag=1)
                session.add(accession)
                session.flush()

                ship = Ships(
                    sequence=clean_seq,
                    sequence_length=len(clean_seq),
                    md5=generate_md5_hash(clean_seq),
                    rev_comp_md5=generate_md5_hash(revcomp(clean_seq)),
                    accession_id=accession.id,
                )
                session.add(ship)
                session.flush()
                ship_id = ship.id

            ensure_ship_has_ssb(session, ship_id)

            taxonomy_id = self._get_or_create_taxonomy(session, data)
            family_id = self._get_or_create_family(session, data)
            navis_id = self._get_or_create_navis(session, data, family_id)
            haplotype_id = self._get_or_create_haplotype(
                session, data, navis_id, family_id
            )
            genome_id = self._get_or_create_genome(session, data, taxonomy_id)

            joined_ship = JoinedShips(
                starshipID=data["starshipID"],
                ship_id=ship_id,
                accession_id=accession.id
                if not duplicate_info.different_taxon
                else ship.accession_id,
                evidence=data["evidence"],
                source=data["source"],
                curated_status=data.get("curated_status", "uncurated"),
                tax_id=taxonomy_id,
                ship_family_id=family_id,
                ship_navis_id=navis_id,
                ship_haplotype_id=haplotype_id,
                genome_id=genome_id,
                captain_id=data.get("captain_id"),
                created_at=datetime.now(),
                updated_at=datetime.now(),
            )
            session.add(joined_ship)
            session.flush()

            if data.get("gff_entries"):
                self._add_gff_entries(
                    session, data["gff_entries"], accession.id, ship_id
                )

            self._create_starship_features(
                session, data, ship_id, data.get("captain_id")
            )

            session.commit()
            logger.info(
                f"Successfully inserted {data['starshipID']} (ship_id={ship_id}, accession={accession_tag})"
            )

            return ship_id, accession_tag

        except Exception as e:
            session.rollback()
            logger.error(f"Database error: {str(e)}")
            raise
        finally:
            session.close()

    def _get_or_create_taxonomy(self, session, data: Dict) -> Optional[int]:
        """Get or create taxonomy record."""
        genus = data.get("genus")
        species = data.get("species")

        if not genus or not species:
            return None

        # Try to find existing
        taxonomy = (
            session.query(Taxonomy)
            .filter(Taxonomy.genus == genus, Taxonomy.species == species)
            .first()
        )

        if taxonomy:
            return taxonomy.id

        # Create new
        taxonomy = Taxonomy(
            name=data.get("tax_name", f"{genus} {species}"),
            genus=genus,
            species=species,
            strain=data.get("strain"),
            order=data.get("tax_order"),
            family=data.get("tax_family"),
        )
        session.add(taxonomy)
        session.flush()
        logger.info(f"Created new taxonomy: {genus} {species}")
        return taxonomy.id

    def _get_or_create_family(self, session, data: Dict) -> Optional[int]:
        """Get or create family name record."""
        family_name = data.get("ship_family")
        if not family_name:
            return None

        family = (
            session.query(FamilyNames)
            .filter(FamilyNames.familyName == family_name)
            .first()
        )

        if family:
            return family.id

        family = FamilyNames(familyName=family_name)
        session.add(family)
        session.flush()
        logger.info(f"Created new family: {family_name}")
        return family.id

    def _get_or_create_navis(
        self, session, data: Dict, family_id: Optional[int]
    ) -> Optional[int]:
        """Get or create navis name record."""
        navis_name = data.get("ship_navis")
        if not navis_name:
            return None

        navis = session.query(Navis).filter(Navis.navis_name == navis_name).first()

        if navis:
            return navis.id

        navis = Navis(navis_name=navis_name, ship_family_id=family_id)
        session.add(navis)
        session.flush()
        logger.info(f"Created new navis: {navis_name}")
        return navis.id

    def _get_or_create_haplotype(
        self, session, data: Dict, navis_id: Optional[int], family_id: Optional[int]
    ) -> Optional[int]:
        """Get or create haplotype name record."""
        haplotype_name = data.get("ship_haplotype")
        if not haplotype_name:
            return None

        haplotype = (
            session.query(Haplotype)
            .filter(Haplotype.haplotype_name == haplotype_name)
            .first()
        )

        if haplotype:
            return haplotype.id

        haplotype = Haplotype(
            haplotype_name=haplotype_name, navis_id=navis_id, ship_family_id=family_id
        )
        session.add(haplotype)
        session.flush()
        logger.info(f"Created new haplotype: {haplotype_name}")
        return haplotype.id

    def _get_or_create_genome(
        self, session, data: Dict, taxonomy_id: Optional[int]
    ) -> Optional[int]:
        """Get or create genome record."""
        assembly_accession = data.get("assembly_accession")
        if not assembly_accession:
            return None

        genome = (
            session.query(Genome)
            .filter(Genome.assembly_accession == assembly_accession)
            .first()
        )

        if genome:
            return genome.id

        genome = Genome(
            assembly_accession=assembly_accession,
            taxonomy_id=taxonomy_id,
            ome=data.get("ome"),
            genomeSource=data.get("genome_source"),
            biosample=data.get("biosample"),
        )
        session.add(genome)
        session.flush()
        logger.info(f"Created new genome: {assembly_accession}")
        return genome.id

    def _add_gff_entries(
        self, session, gff_entries: List[Dict], accession_id: int, ship_id: int
    ):
        """Add GFF annotation entries."""
        for entry in gff_entries:
            gff = Gff(
                seqid=entry.get("seqid"),
                source=entry.get("source"),
                type=entry.get("type"),
                start=entry.get("start"),
                end=entry.get("end"),
                score=entry.get("score"),
                strand=entry.get("strand"),
                phase=entry.get("phase"),
                attributes=entry.get("attributes"),
                accession_id=accession_id,
                ship_id=ship_id,
            )
            session.add(gff)

    def _create_starship_features(
        self, session, data: Dict, ship_id: int, captain_id: Optional[int]
    ):
        """Create StarshipFeatures record if location data is provided."""
        contig_id = data.get("contig_id")
        element_start = data.get("element_start")
        element_end = data.get("element_end")

        # Only create if we have the minimum required fields
        if not (contig_id and element_start and element_end):
            return

        # elementLength must equal ABS(elementEnd - elementBegin) + 1 (always compute from coordinates)
        # abs() handles negative-strand features where element_start > element_end
        element_length = abs(element_end - element_start) + 1

        features = StarshipFeatures(
            contigID=contig_id,
            starshipID=data.get("starshipID"),
            captainID=data.get("captain_id"),  # Note: this might need to be fetched
            elementBegin=element_start,
            elementEnd=element_end,
            elementLength=element_length,
            strand=data.get("element_strand"),
            ship_id=ship_id,
            captain_id=captain_id,
        )
        session.add(features)
        logger.info(f"Created StarshipFeatures for {data.get('starshipID')}")

    def process_batch(self, submissions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Process multiple submissions.

        Args:
            submissions: List of submission data dictionaries

        Returns:
            List of result dictionaries
        """
        results = []
        total = len(submissions)

        logger.info(f"Processing batch of {total} submissions...")

        for idx, submission in enumerate(submissions, 1):
            logger.info(
                f"Processing {idx}/{total}: {submission.get('starshipID', 'unknown')}"
            )
            result = self.process_submission(submission)
            results.append(result)

        self._print_summary()
        return results

    def _print_summary(self):
        """Print processing summary."""
        logger.info("=" * 80)
        logger.info("SUBMISSION PROCESSING SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Total submissions: {self.stats['total']}")
        logger.info(f"Validated: {self.stats['validated']}")
        logger.info(f"Inserted: {self.stats['inserted']}")
        logger.info(f"Duplicates (same taxon): {self.stats['duplicates_same_taxon']}")
        logger.info(
            f"Duplicates (different taxon): {self.stats['duplicates_diff_taxon']}"
        )
        logger.info(f"Errors: {len(self.stats['errors'])}")

        if self.stats["errors"]:
            logger.info("\nErrors encountered:")
            for error in self.stats["errors"][:10]:  # Show first 10
                logger.info(f"  - {error}")
            if len(self.stats["errors"]) > 10:
                logger.info(f"  ... and {len(self.stats['errors']) - 10} more")


# ============================================================================
# CONVENIENCE FUNCTIONS
# ============================================================================


def parse_fasta_file(fasta_path: str) -> Dict[str, str]:
    """
    Parse FASTA file and return dict of header -> sequence.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dict mapping header (without '>') to sequence

    Note:
        Only uses the first word of the header (before any whitespace)
        as the key, which matches standard FASTA ID format.
    """
    import screed

    logger.info(f"Reading FASTA file: {fasta_path}")
    sequences = {}

    with screed.open(fasta_path) as fasta:
        for record in fasta:
            # Extract just the ID (first word before whitespace)
            full_header = record.name
            header_id = full_header.split()[0] if full_header else full_header
            sequence = str(record.sequence)
            sequences[header_id] = sequence

    logger.info(f"Loaded {len(sequences)} sequences from FASTA")
    return sequences


def merge_fasta_with_tsv(
    fasta_path: str, tsv_path: str, header_column: str = "starshipID"
) -> List[Dict[str, Any]]:
    """
    Merge FASTA sequences with TSV metadata.

    Args:
        fasta_path: Path to FASTA file with sequences
        tsv_path: Path to TSV file with metadata (no sequence column)
        header_column: Name of TSV column that matches FASTA headers

    Returns:
        List of merged submission dictionaries
    """
    # Parse FASTA
    sequences = parse_fasta_file(fasta_path)

    # Read TSV
    logger.info(f"Reading metadata from {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t")

    if header_column not in df.columns:
        raise ValueError(
            f"Column '{header_column}' not found in TSV. Available columns: {list(df.columns)}"
        )

    # Merge sequences with metadata
    submissions = []
    unmatched_headers = []

    for _, row in df.iterrows():
        row_dict = row.to_dict()
        header = str(row_dict[header_column])

        if header in sequences:
            row_dict["sequence"] = sequences[header]

            # Convert pandas NaN to None for proper validation
            # NaN values fail isinstance() checks
            for key, value in row_dict.items():
                if pd.isna(value):
                    row_dict[key] = None

            submissions.append(row_dict)
        else:
            unmatched_headers.append(header)
            logger.warning(f"No sequence found for header: {header}")

    logger.info(f"Matched {len(submissions)} sequences with metadata")

    if unmatched_headers:
        logger.warning(f"Unmatched headers: {len(unmatched_headers)}")
        if len(unmatched_headers) <= 10:
            for header in unmatched_headers:
                logger.warning(f"  - {header}")
        else:
            logger.warning(f"  (showing first 10: {unmatched_headers[:10]})")

    # Check for extra sequences in FASTA
    tsv_headers = set(df[header_column].astype(str))
    fasta_headers = set(sequences.keys())
    extra_in_fasta = fasta_headers - tsv_headers

    if extra_in_fasta:
        logger.warning(
            f"Found {len(extra_in_fasta)} sequences in FASTA without metadata in TSV"
        )
        if len(extra_in_fasta) <= 10:
            for header in extra_in_fasta:
                logger.warning(f"  - {header}")

    return submissions


def process_from_tsv(tsv_path: str, dry_run: bool = False) -> List[Dict[str, Any]]:
    """
    Process submissions from TSV file (with embedded sequences).

    Expected columns:
    - Required: sequence, starshipID, evidence, source
    - Optional: genus, species, ship_family, etc.

    Args:
        tsv_path: Path to TSV file
        dry_run: If True, don't actually insert to database

    Returns:
        List of result dictionaries
    """
    logger.info(f"Reading submissions from {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t")

    # Convert DataFrame to list of dicts
    submissions = df.to_dict("records")

    processor = SubmissionProcessor(dry_run=dry_run)
    results = processor.process_batch(submissions)

    return results


def process_from_fasta_and_tsv(
    fasta_path: str,
    tsv_path: str,
    dry_run: bool = False,
    header_column: str = "starshipID",
) -> List[Dict[str, Any]]:
    """
    Process submissions from separate FASTA and TSV files.

    This is the RECOMMENDED method for batch submissions.

    Args:
        fasta_path: Path to FASTA file with sequences
        tsv_path: Path to TSV file with metadata (without sequence column)
        dry_run: If True, don't actually insert to database
        header_column: TSV column name that matches FASTA headers (default: 'starshipID')

    Returns:
        List of result dictionaries

    Example:
        # FASTA file (sequences.fasta):
        >SS-1.1
        ATGCATGC...
        >SS-1.2
        GCTAGCTA...

        # TSV file (metadata.tsv):
        starshipID  evidence  source      genus      species
        SS-1.1      starfish  pub2024     Fusarium   oxysporum
        SS-1.2      starfish  pub2024     Fusarium   solani

        # Process:
        results = process_from_fasta_and_tsv('sequences.fasta', 'metadata.tsv')
    """
    # Merge FASTA with TSV
    submissions = merge_fasta_with_tsv(fasta_path, tsv_path, header_column)

    # Process submissions
    processor = SubmissionProcessor(dry_run=dry_run)
    results = processor.process_batch(submissions)

    return results


def process_single(
    sequence: str, starshipID: str, evidence: str, source: str, **kwargs
) -> Dict[str, Any]:
    """
    Process a single submission (convenience function).

    Args:
        sequence: Sequence string
        starshipID: Starship identifier
        evidence: Evidence method
        source: Data source
        **kwargs: Additional optional fields

    Returns:
        Result dictionary
    """
    data = {
        "sequence": sequence,
        "starshipID": starshipID,
        "evidence": evidence,
        "source": source,
        **kwargs,
    }

    processor = SubmissionProcessor(dry_run=False)
    return processor.process_submission(data)
