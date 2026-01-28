import subprocess
import os
import fcntl
import pandas as pd
from sourmash import (
    MinHash,
    SourmashSignature,
    save_signatures,
    load_file_as_signatures,
)
import screed

from src.config.cache import cache
from src.config.settings import BLAST_DB_PATHS
from src.utils.seq_utils import create_ncbi_style_header, write_fasta

from src.config.logging import get_logger

logger = get_logger(__name__)


class FileLock:
    """Simple filesystem lock using fcntl (Linux/Unix)."""

    def __init__(self, lock_path: str):
        self.lock_path = lock_path
        self._fh = None

    def __enter__(self):
        os.makedirs(os.path.dirname(self.lock_path), exist_ok=True)
        self._fh = open(self.lock_path, "w")
        fcntl.flock(self._fh.fileno(), fcntl.LOCK_EX)
        return self

    def __exit__(self, exc_type, exc, tb):
        try:
            if self._fh is not None:
                fcntl.flock(self._fh.fileno(), fcntl.LOCK_UN)
                self._fh.close()
        finally:
            self._fh = None


def create_blast_database(fasta_path, dbtype):
    """
    Create a BLAST database from the provided FASTA file.

    - Uses the FASTA basename (without extension) as the DB prefix
    - Proactively removes any stale BLAST index/taxonomy files to avoid taxonomy LMDB conflicts
    """
    # Use basename without extension as the BLAST DB prefix
    out_prefix, _ = os.path.splitext(fasta_path)

    # Best-effort cleanup of old index/taxonomy files (prevents taxonomy-related crashes)
    stale_suffixes = [
        ".nhr",
        ".nin",
        ".nsq",
        ".ndb",
        ".not",
        ".ntf",
        ".nto",
        ".njs",
    ]
    for suffix in stale_suffixes:
        stale_path = f"{out_prefix}{suffix}"
        try:
            if os.path.exists(stale_path):
                os.remove(stale_path)
        except Exception:
            # Ignore cleanup failures
            pass

    cmd = [
        "makeblastdb",
        "-in",
        fasta_path,
        "-input_type",
        "fasta",
        "-dbtype",
        dbtype,
        "-out",
        out_prefix,
    ]

    logger.debug(f"Creating BLAST database with command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.debug("BLAST database created successfully.")
    except subprocess.CalledProcessError as e:
        stderr = e.stderr or ""
        stdout = e.stdout or ""
        logger.error(
            "Failed to create BLAST database (exit %s).\nSTDOUT:\n%s\nSTDERR:\n%s",
            e.returncode,
            stdout.strip(),
            stderr.strip(),
        )
        raise


def blast_db_exists(blastdb):
    """
    Check if the BLAST database exists by looking for required file extensions.
    """
    # Determine the expected BLAST DB prefix (basename without extension)
    db_prefix, _ = os.path.splitext(blastdb)

    # Required core files for BLAST DB existence check
    required_suffixes = [".nhr", ".nin", ".nsq"]

    for suffix in required_suffixes:
        expected_path = f"{db_prefix}{suffix}"
        if not os.path.exists(expected_path):
            return False

    return True


def create_sourmash_signatures(fasta_path, seq_type="nucl"):
    """
    Create sourmash signatures from a FASTA file and save to disk.

    Args:
        fasta_path: Path to input FASTA file
        seq_type: 'nucl' for nucleotide or 'prot' for protein

    Returns:
        Path to the created signature file (.sig)
    """
    sig_path = fasta_path + ".sig"

    # Set parameters based on sequence type
    if seq_type == "nucl":
        k = 21
        scaled = 1000
        is_protein = False
    else:
        k = 7
        scaled = 100
        is_protein = True

    logger.debug(f"Creating sourmash signatures for {fasta_path} (type={seq_type})")

    try:
        signatures = []
        seq_count = 0

        # Read FASTA and create signatures
        for record in screed.open(fasta_path):
            mh = MinHash(n=0, ksize=k, scaled=scaled, is_protein=is_protein)

            # Add sequence data
            if is_protein:
                mh.add_protein(record.sequence)
            else:
                mh.add_sequence(record.sequence, force=True)

            # Create signature with sequence name
            sig = SourmashSignature(mh, name=record.name)
            signatures.append(sig)
            seq_count += 1

        # Save all signatures to file
        with open(sig_path, "w") as f:
            save_signatures(signatures, f)

        logger.info(f"Created {seq_count} sourmash signatures: {sig_path}")
        return sig_path

    except Exception as e:
        logger.error(f"Failed to create sourmash signatures: {e}")
        raise


def sourmash_sig_exists(fasta_path):
    """
    Check if sourmash signature file exists for a given FASTA file.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        bool: True if .sig file exists
    """
    sig_path = fasta_path + ".sig"
    return os.path.exists(sig_path)


def load_sourmash_signatures(fasta_path):
    """
    Load pre-computed sourmash signatures from disk.

    Args:
        fasta_path: Path to FASTA file (signature file is fasta_path + '.sig')

    Returns:
        List of (seq_name, signature) tuples or None if file doesn't exist
    """
    sig_path = fasta_path + ".sig"

    if not os.path.exists(sig_path):
        logger.warning(f"Signature file not found: {sig_path}")
        return None

    try:
        signatures = []
        for sig in load_file_as_signatures(sig_path):
            signatures.append((sig.name, sig))

        logger.debug(f"Loaded {len(signatures)} signatures from {sig_path}")
        return signatures

    except Exception as e:
        logger.error(f"Failed to load signatures from {sig_path}: {e}")
        return None


def append_sourmash_signature(fasta_path, sequence, sequence_name, seq_type="nucl"):
    """
    Append a single signature to an existing signature file.

    This is more efficient than rebuilding the entire signature file
    when adding individual sequences.

    Args:
        fasta_path: Base path to FASTA file (signature is fasta_path + '.sig')
        sequence: DNA/protein sequence string
        sequence_name: Name/ID for the sequence (e.g., accession tag)
        seq_type: 'nucl' or 'prot'

    Returns:
        bool: True if successful, False otherwise
    """
    sig_path = fasta_path + ".sig"

    # Set parameters based on sequence type
    if seq_type == "nucl":
        k = 21
        scaled = 1000
        is_protein = False
    else:
        k = 7
        scaled = 100
        is_protein = True

    try:
        # Create signature for new sequence
        mh = MinHash(n=0, ksize=k, scaled=scaled, is_protein=is_protein)

        if is_protein:
            mh.add_protein(sequence)
        else:
            mh.add_sequence(sequence, force=True)

        new_sig = SourmashSignature(mh, name=sequence_name)

        # Load existing signatures if file exists
        existing_signatures = []
        if os.path.exists(sig_path):
            try:
                for sig in load_file_as_signatures(sig_path):
                    # Skip if this sequence name already exists (update case)
                    if sig.name != sequence_name:
                        existing_signatures.append(sig)
            except Exception as e:
                logger.warning(
                    f"Could not load existing signatures from {sig_path}: {e}"
                )

        # Combine existing + new signature
        all_signatures = existing_signatures + [new_sig]

        # Save all signatures back to file
        with open(sig_path, "w") as f:
            save_signatures(all_signatures, f)

        logger.debug(
            f"Appended signature for {sequence_name} to {sig_path} (total: {len(all_signatures)})"
        )
        return True

    except Exception as e:
        logger.error(f"Failed to append signature for {sequence_name}: {e}")
        return False


def append_sourmash_signatures_batch(fasta_path, sequences, seq_type="nucl"):
    """
    Append multiple signatures to an existing signature file in one operation.

    More efficient than calling append_sourmash_signature() repeatedly.

    Args:
        fasta_path: Base path to FASTA file (signature is fasta_path + '.sig')
        sequences: List of (sequence_name, sequence_string) tuples
        seq_type: 'nucl' or 'prot'

    Returns:
        int: Number of signatures successfully appended
    """
    sig_path = fasta_path + ".sig"

    # Set parameters
    if seq_type == "nucl":
        k = 21
        scaled = 1000
        is_protein = False
    else:
        k = 7
        scaled = 100
        is_protein = True

    try:
        # Create signatures for new sequences
        new_signatures = []
        logger.debug(
            f"Creating signatures for {len(sequences)} sequences: {[name for name, _ in sequences]}"
        )
        for seq_name, sequence in sequences:
            mh = MinHash(n=0, ksize=k, scaled=scaled, is_protein=is_protein)

            if is_protein:
                mh.add_protein(sequence)
            else:
                mh.add_sequence(sequence, force=True)

            sig = SourmashSignature(mh, name=seq_name)
            new_signatures.append(sig)

        # Load existing signatures
        existing_signatures = []
        replaced_signatures = []
        new_seq_names = {name for name, _ in sequences}

        logger.debug(f"Loading existing signatures from {sig_path}")

        original_count = 0
        if os.path.exists(sig_path):
            try:
                for sig in load_file_as_signatures(sig_path):
                    original_count += 1
                    # Skip if updating (new batch contains this name)
                    if sig.name not in new_seq_names:
                        existing_signatures.append(sig)
                    else:
                        replaced_signatures.append(sig.name)
                logger.debug(
                    f"Loaded {original_count} existing signatures, {len(replaced_signatures)} will be replaced"
                )
            except Exception as e:
                logger.warning(f"Could not load existing signatures: {e}")
        else:
            logger.debug(f"No existing signature file found at {sig_path}")

        # Combine and save
        all_signatures = existing_signatures + new_signatures
        logger.debug(
            f"Writing {len(all_signatures)} total signatures ({len(existing_signatures)} kept + {len(new_signatures)} new/updated)"
        )

        with open(sig_path, "w") as f:
            save_signatures(all_signatures, f)

        # Calculate truly new vs replaced
        num_replaced = len(replaced_signatures)
        num_truly_new = len(new_signatures) - num_replaced

        if num_replaced > 0 and num_truly_new > 0:
            logger.info(
                f"Updated {sig_path}: added {num_truly_new} new, replaced {num_replaced} existing (total: {len(all_signatures)})"
            )
        elif num_replaced > 0:
            logger.info(
                f"Updated {sig_path}: replaced {num_replaced} existing signatures (total: {len(all_signatures)})"
            )
        else:
            logger.info(
                f"Appended {num_truly_new} signatures to {sig_path} (total: {len(all_signatures)})"
            )

        return num_truly_new  # Return only truly new signatures

    except Exception as e:
        logger.error(f"Failed to append batch signatures: {e}")
        return 0


def check_signature_file_status(fasta_path):
    """
    Check the status of a signature file and report statistics.

    Args:
        fasta_path: Base path to FASTA file

    Returns:
        dict: Statistics about the signature file
    """
    sig_path = fasta_path + ".sig"

    stats = {"exists": False, "num_signatures": 0, "accessions": [], "path": sig_path}

    if os.path.exists(sig_path):
        stats["exists"] = True
        try:
            from sourmash import load_file_as_signatures

            for sig in load_file_as_signatures(sig_path):
                stats["num_signatures"] += 1
                if sig.name:
                    stats["accessions"].append(sig.name)
            logger.info(
                f"Signature file {sig_path} contains {stats['num_signatures']} signatures"
            )
        except Exception as e:
            logger.error(f"Error reading signature file {sig_path}: {e}")
            stats["error"] = str(e)
    else:
        logger.warning(f"Signature file not found: {sig_path}")

    return stats


def update_signature_file_from_database(dataset_name="all", seq_type="nucl"):
    """
    Rebuild signature file from current database state.

    Use this as a fallback if incremental updates fail or to ensure
    signatures are in sync with the database.

    Args:
        dataset_name: 'all' or 'curated'
        seq_type: 'nucl' or 'prot'

    Returns:
        bool: True if successful, False otherwise
    """
    from src.database.sql_manager import fetch_ships
    import tempfile
    import shutil

    try:
        fasta_path = BLAST_DB_PATHS["ship"][dataset_name]["nucl"]

        # Fetch current ships
        curated = dataset_name == "curated"
        ships = fetch_ships(dereplicate=True, with_sequence=True, curated=curated)

        if ships.empty:
            logger.warning(f"No ships found for {dataset_name} dataset")
            return False

        # Write current ships to temp FASTA
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
            for _, row in ships.iterrows():
                accession = row.get("accession_display", row.get("accession_tag"))
                sequence = row["sequence"]
                if accession and sequence:
                    tmp.write(f">{accession}\n{sequence}\n")
            tmp_path = tmp.name

        # Recreate signatures from temp FASTA
        create_sourmash_signatures(tmp_path, seq_type=seq_type)

        # Move temp .sig to actual location
        shutil.move(tmp_path + ".sig", fasta_path + ".sig")
        os.unlink(tmp_path)

        logger.info(f"Rebuilt signature file for {dataset_name} dataset")
        return True

    except Exception as e:
        logger.error(f"Failed to update signature file: {e}")
        return False


def create_reference_fasta_for_minimap2(fasta_path):
    """
    Create a reference FASTA file for minimap2 alignment from existing ship sequences.

    This is used for the "check_contained_match" step in accession assignment.
    Pre-creating this file avoids writing all reference sequences to a temp file
    on every alignment check.

    Args:
        fasta_path: Path to the ship FASTA file

    Returns:
        str: Path to the created reference FASTA file (.ref.fasta)
    """
    ref_fasta_path = fasta_path + ".ref.fasta"

    try:
        # The reference file is just a copy of the main FASTA
        # We keep it separate so it can be updated independently
        if os.path.exists(fasta_path):
            import shutil

            shutil.copy2(fasta_path, ref_fasta_path)
            logger.info(f"Created reference FASTA for minimap2: {ref_fasta_path}")
            return ref_fasta_path
        else:
            logger.warning(f"Source FASTA not found: {fasta_path}")
            return None

    except Exception as e:
        logger.error(f"Failed to create reference FASTA: {e}")
        return None


def append_to_reference_fasta(fasta_path, sequences):
    """
    Append sequences to the reference FASTA file for minimap2.

    Args:
        fasta_path: Base path to FASTA file (reference is fasta_path + '.ref.fasta')
        sequences: List of (accession, sequence) tuples to append

    Returns:
        int: Number of sequences appended
    """
    ref_fasta_path = fasta_path + ".ref.fasta"

    try:
        # Deduplicate sequences by accession
        seen = {}
        for accession, sequence in sequences:
            seen[accession] = sequence
        deduped_sequences = list(seen.items())

        # Read existing accessions to avoid duplicates
        existing_accessions = set()
        if os.path.exists(ref_fasta_path):
            with open(ref_fasta_path, "r") as f:
                for line in f:
                    if line.startswith(">"):
                        accession = line[1:].split()[0]  # Get first word after >
                        existing_accessions.add(accession)

        # Append only new sequences
        appended_count = 0
        with open(ref_fasta_path, "a") as f:
            for accession, sequence in deduped_sequences:
                if accession not in existing_accessions:
                    f.write(f">{accession}\n{sequence}\n")
                    appended_count += 1

        num_skipped = len(deduped_sequences) - appended_count

        if appended_count > 0 and num_skipped > 0:
            logger.info(
                f"Updated {ref_fasta_path}: appended {appended_count} new, skipped {num_skipped} existing"
            )
        elif appended_count > 0:
            logger.info(
                f"Appended {appended_count} sequences to reference FASTA: {ref_fasta_path}"
            )
        else:
            logger.info(
                f"All {len(deduped_sequences)} sequences already exist in {ref_fasta_path} (skipped)"
            )

        return appended_count

    except Exception as e:
        logger.error(f"Failed to append to reference FASTA: {e}")
        return 0


def create_diamond_database(fasta_path, threads=2):
    # Create output path with .dmnd extension
    diamond_db = fasta_path + ".dmnd"

    cmd = [
        "diamond",
        "makedb",
        "--in",
        fasta_path,
        "--db",
        diamond_db,
        "--threads",
        str(threads),
    ]

    logger.debug(f"Creating Diamond database with command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.debug("Diamond database created successfully.")
        return diamond_db  # Return the path to the created Diamond database
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to create Diamond database: {e.stderr}")
        raise


def create_dbs():
    from src.database.sql_manager import (
        fetch_captains,
        fetch_ships,
        fetch_meta_data,
    )

    datasets = [
        ("all", False),
        ("curated", True),
    ]

    # Serialize rebuilds across processes using a file lock
    lock_file = os.path.join(
        os.path.dirname(BLAST_DB_PATHS["ship"]["all"]["nucl"]),
        ".rebuild.lock",
    )

    with FileLock(lock_file):
        for dataset_name, curated_filter in datasets:
            logger.info(f"Creating {dataset_name} ship BLAST database...")

            ship_fasta_path = BLAST_DB_PATHS["ship"][dataset_name]["nucl"]
            ship_fasta_dir = os.path.dirname(ship_fasta_path)
            os.makedirs(ship_fasta_dir, exist_ok=True)

            ship_sequences = fetch_ships(
                dereplicate=True, with_sequence=True, curated=curated_filter
            )

            if ship_sequences.empty:
                logger.error(f"No sequences found for {dataset_name} dataset")
                continue

            ship_metadata = fetch_meta_data(
                accessions=ship_sequences["ship_accession_tag"].tolist()
            )
            if ship_metadata.empty:
                logger.error(f"No metadata found for {dataset_name} dataset")
                continue

            ship_sequences_dict = {}

            for _, row in ship_sequences.iterrows():
                # Convert Series to scalar values to avoid unhashable type errors
                accession_tag = (
                    str(row["ship_accession_tag"])
                    if pd.notna(row["ship_accession_tag"])
                    else None
                )
                sequence = str(row["sequence"]) if pd.notna(row["sequence"]) else None

                if not accession_tag or not sequence:
                    continue

                metadata_rows = ship_metadata[
                    ship_metadata["ship_accession_tag"] == accession_tag
                ]
                if not metadata_rows.empty:
                    # Ensure all values are scalars, not pandas objects
                    metadata = {}
                    for col in metadata_rows.columns:
                        val = metadata_rows.iloc[0][col]
                        # Convert to Python types
                        if pd.isna(val):
                            metadata[col] = None
                        else:
                            metadata[col] = val.item() if hasattr(val, "item") else val

                    header = create_ncbi_style_header(metadata)
                    if header and sequence:
                        # Ensure header is a string to avoid unhashable type errors
                        header = str(header) if header else None
                        if header:
                            ship_sequences_dict[header] = sequence

            logger.info(f"Rebuilding {dataset_name} BLAST database...")
            write_fasta(ship_sequences_dict, ship_fasta_path)
            create_blast_database(ship_fasta_path, "nucl")
            logger.info(f"{dataset_name} BLAST database rebuilt successfully")

            # Create sourmash signatures
            logger.info(f"Creating sourmash signatures for {dataset_name} dataset...")
            try:
                create_sourmash_signatures(ship_fasta_path, seq_type="nucl")
                logger.info(f"{dataset_name} sourmash signatures created successfully")
            except Exception as e:
                logger.error(
                    f"Failed to create sourmash signatures for {dataset_name}: {e}"
                )

            # Create reference FASTA for minimap2 alignments
            logger.info(
                f"Creating reference FASTA for minimap2 ({dataset_name} dataset)..."
            )
            try:
                create_reference_fasta_for_minimap2(ship_fasta_path)
                logger.info(f"{dataset_name} reference FASTA created successfully")
            except Exception as e:
                logger.error(
                    f"Failed to create reference FASTA for {dataset_name}: {e}"
                )

    # Create captain database
    captain_fasta_path = BLAST_DB_PATHS["gene"]["tyr"]["prot"]
    captain_fasta_dir = os.path.dirname(captain_fasta_path)
    os.makedirs(captain_fasta_dir, exist_ok=True)

    captain_sequences = cache.get("all_captains")
    if captain_sequences is None:
        captain_sequences = fetch_captains(dereplicate=True, with_sequence=True)

    captain_sequences_dict = {}
    for index, row in captain_sequences.iterrows():
        captain_id_val = row.get("captain_id_col", row.get("captainID"))
        sequence_val = row["sequence"]

        if pd.notna(captain_id_val) and pd.notna(sequence_val):
            accession = str(captain_id_val)
            sequence = str(sequence_val)
            captain_sequences_dict[accession] = sequence

    # Captain DBs: skip if already exist
    if blast_db_exists(captain_fasta_path):
        logger.info("Captain BLAST DB already exists; skipping rebuild")
    else:
        write_fasta(captain_sequences_dict, captain_fasta_path)
        create_blast_database(captain_fasta_path, "prot")
        create_diamond_database(captain_fasta_path)

        # Create sourmash signatures for captain sequences
        logger.info("Creating sourmash signatures for captain sequences...")
        try:
            create_sourmash_signatures(captain_fasta_path, seq_type="prot")
            logger.info("Captain sourmash signatures created successfully")
        except Exception as e:
            logger.error(f"Failed to create captain sourmash signatures: {e}")
