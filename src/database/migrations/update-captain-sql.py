import sqlite3
import os
import sys
from pathlib import Path
from Bio import SeqIO
from contextlib import closing
from src.utils.blast_utils import run_diamond, run_hmmer

# Add project root to Python path for proper imports
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Configuration - Use absolute paths based on script location
script_dir = Path(__file__).parent  # src/database/
db_dir = script_dir / "db"
DB_PATH = str(db_dir / "starbase.sqlite")
FAA_PATH = str(db_dir / "Starships" / "captain" / "tyr" / "faa" / "faa")  # Fixed path
FNA_PATH = str(db_dir / "Starships" / "captain" / "tyr" / "fna")
PROT_BLASTDB = str(db_dir / "captain" / "tyr" / "faa" / "blastdb" / "captains.faa")
NUCL_BLASTDB = str(
    db_dir / "captain" / "tyr" / "fna" / "blastdb" / "concatenated.dedup.fna"
)


class CaptainSequenceUpdater:
    def __init__(self):
        self.db_path = DB_PATH
        self.search_paths = [FAA_PATH, FNA_PATH]

    def extract_sequence_from_file(self, fasta_file, target_header):
        """Extract sequence from a single FASTA file based on header match."""
        try:
            with open(fasta_file, "r") as f:
                records = list(SeqIO.parse(f, "fasta"))
                if not records:
                    return None, None

                # Look for exact or partial header matches
                for record in records:
                    if target_header in record.id or record.id in target_header:
                        return str(record.seq), record.id

                return None, None
        except (FileNotFoundError, IOError):
            return None, None

    def search_existing_sequences(self, directory_path, target_headers):
        """Search for existing sequences in directory using multiple possible headers."""
        if not os.path.isdir(directory_path):
            print(f"Directory does not exist: {directory_path}")
            return None, None, None

        for target_header in target_headers:
            if target_header is None:
                continue

            print(f"Searching for {target_header} in {directory_path}")

            # First try to find files that contain the target_header in the filename
            for filename in os.listdir(directory_path):
                if target_header in filename and filename.endswith(
                    (".fa", ".faa", ".fna", ".fasta")
                ):
                    file_path = os.path.join(directory_path, filename)
                    if os.path.isfile(file_path):
                        print(f"Found matching file: {filename}")
                        # Read the first sequence from this file
                        try:
                            with open(file_path, "r") as f:
                                records = list(SeqIO.parse(f, "fasta"))
                                if records:
                                    return (
                                        str(records[0].seq),
                                        records[0].id,
                                        "existing_file",
                                    )
                        except Exception as e:
                            print(f"Error reading file {filename}: {e}")
                            continue

            # If no filename match, search within files
            for filename in os.listdir(directory_path):
                file_path = os.path.join(directory_path, filename)
                if os.path.isfile(file_path) and filename.endswith(
                    (".fa", ".faa", ".fna", ".fasta")
                ):
                    sequence, header = self.extract_sequence_from_file(
                        file_path, str(target_header)
                    )
                    if sequence:
                        print(f"Found sequence in file: {filename}")
                        return sequence, header, "existing_file"

        return None, None, None

    def run_blast_search(self, target_header):
        """Run BLAST/HMMER search as fallback when sequences aren't found in existing files."""
        try:
            print(f"Running BLAST search for {target_header}")

            # Try DIAMOND first - use the protein BLAST database
            protein_filename = run_diamond(
                query_type="prot",
                input_gene="tyr",
                input_eval=0.01,
                query_fasta=PROT_BLASTDB,
                threads=2,
            )

            if protein_filename and os.path.exists(protein_filename):
                sequence, header = self.extract_sequence_from_file(
                    protein_filename, str(target_header)
                )
                if sequence:
                    print("Found sequence via DIAMOND search")
                    return sequence, header, "diamond_search"
            else:
                # Try HMMER if DIAMOND fails
                # but grab ship sequence first
                with closing(sqlite3.connect(self.db_path)) as conn:
                    cursor = conn.cursor()
                    query = """
                    SELECT sequence 
                    FROM ships 
                    LEFT JOIN joined_ships ON ships.id = joined_ships.ship_id
                    LEFT JOIN captains ON joined_ships.captain_id = captains.id
                    WHERE captains.captainID = ? AND captains.sequence IS NULL
                    """
                    cursor.execute(query, (target_header,))
                    ship_sequence = cursor.fetchone()
                    if ship_sequence:
                        sequence = ship_sequence[0]
                        return sequence, target_header, "existing_file"
                    else:
                        print(f"No ship sequence found for captain {target_header}")
                    # write a temporary file with the ship sequence
                    tmp_ship_seq = f"tmp/ship_sequence_{target_header}.fasta"
                    with open(tmp_ship_seq, "w") as f:
                        f.write(f">{target_header}\n{ship_sequence[0]}")

                    protein_filename = run_hmmer(
                        query_type="prot",
                        input_gene="tyr",
                        input_eval=0.01,
                        query_fasta=tmp_ship_seq,
                        threads=2,
                    )

            if protein_filename and os.path.exists(protein_filename):
                sequence, header = self.extract_sequence_from_file(
                    protein_filename, str(target_header)
                )
                if sequence:
                    print("Found sequence via HMMER search")
                    return sequence, header, "hmmer_search"

        except Exception as e:
            print(f"Search failed for {target_header}: {e}")

        return None, None, None

    def find_sequence_for_captain(self, captain_identifiers):
        """
        Find sequence for a captain using the priority:
        1. Search existing files in faa directory
        2. Search existing files in fna directory
        3. Run BLAST/HMMER searches
        """
        # Filter out None values and convert to strings
        valid_identifiers = [
            str(id_val) for id_val in captain_identifiers if id_val is not None
        ]

        if not valid_identifiers:
            return None, None, None

        print(f"Searching for identifiers: {valid_identifiers}")

        # Step 1 & 2: Search existing sequence files
        for search_path in self.search_paths:
            sequence, header, evidence = self.search_existing_sequences(
                search_path, valid_identifiers
            )
            if sequence:
                return sequence, header, evidence

        # Step 3: Run BLAST/HMMER searches as fallback
        for identifier in valid_identifiers:
            sequence, header, evidence = self.run_blast_search(identifier)
            if sequence:
                return sequence, header, evidence

        return None, None, None

    def resolve_captain_ship_ids(self, cursor, row_data):
        """Resolve captain and ship IDs between tables."""
        (
            accession_tag,
            accession_id,
            captainID_c,
            captainID_j,
            captain_id,
            starshipID,
            ship_id_int,
            existing_sequence,
            captain_row_id,
        ) = row_data

        # Determine ship_id if missing
        resolved_ship_id = ship_id_int
        if ship_id_int is None:
            # Try to find ship_id using captainID
            cursor.execute(
                "SELECT ship_id FROM joined_ships WHERE captainID = ?", (captainID_c,)
            )
            result = cursor.fetchone()
            if result:
                resolved_ship_id = result[0]
            else:
                cursor.execute(
                    "SELECT ship_id FROM joined_ships WHERE captain_id = ?",
                    (captain_id,),
                )
                result = cursor.fetchone()
                if result:
                    resolved_ship_id = result[0]

        # Determine best captainID to use
        captain_id_to_use = captainID_c or captainID_j or starshipID

        return resolved_ship_id, captain_id_to_use

    def update_database_tables(
        self,
        cursor,
        captain_row_id,
        accession_id,
        sequence,
        header,
        evidence,
        ship_id,
        captain_id,
    ):
        """Update both captains and joined_ships tables with new data."""
        try:
            # Update captains table
            cursor.execute(
                """
                UPDATE captains 
                SET sequence = ?, captainID = ?, ship_id = ?, evidence = ? 
                WHERE id = ?
            """,
                (sequence, captain_id, ship_id, evidence, captain_row_id),
            )

            # Update joined_ships table
            cursor.execute(
                """
                UPDATE joined_ships 
                SET captain_id = ? 
                WHERE ship_id = ?
            """,
                (captain_row_id, accession_id),
            )

            return True
        except sqlite3.Error as e:
            print(f"Database update failed for captain {captain_id}: {e}")
            return False

    def run_update(self):
        """Main function to update captain sequences."""
        with closing(sqlite3.connect(self.db_path)) as conn:
            cursor = conn.cursor()

            # Query for rows with missing sequences
            query = """
            SELECT a.accession_tag, a.id, c.captainID, j.captainID, j.captain_id, 
                   j.starshipID, c.ship_id, c.sequence, c.id
            FROM captains c
            JOIN joined_ships j ON c.id = j.captain_id
            JOIN accessions a ON c.ship_id = a.id
            WHERE c.sequence IS NULL OR c.ship_id IS NULL
            """

            cursor.execute(query)
            rows = cursor.fetchall()
            print(f"Found {len(rows)} rows with missing sequences to update")

            updated_count = 0
            failed_count = 0

            for row in rows:
                (
                    accession_tag,
                    accession_id,
                    captainID_c,
                    captainID_j,
                    captain_id,
                    starshipID,
                    ship_id_int,
                    existing_sequence,
                    captain_row_id,
                ) = row

                print(f"Processing captain row ID: {captain_row_id}")

                # Resolve IDs
                resolved_ship_id, captain_id_to_use = self.resolve_captain_ship_ids(
                    cursor, row
                )

                if resolved_ship_id is None:
                    print(f"Could not resolve ship_id for captain {captain_id_to_use}")
                    failed_count += 1
                    continue

                # Search for sequence using priority workflow
                captain_identifiers = [captainID_c, captainID_j, starshipID]
                # make the list unique
                captain_identifiers = list(set(captain_identifiers))

                sequence, header, evidence = self.find_sequence_for_captain(
                    captain_identifiers
                )

                if sequence is None:
                    print(
                        f"No sequence found for captain identifiers: {captain_identifiers}"
                    )
                    failed_count += 1
                    continue

                # Use header as captain_id if original is missing and header doesn't match accession
                final_captain_id = captain_id_to_use
                if (
                    captain_id_to_use is None
                    and header
                    and accession_tag
                    and accession_tag not in header
                ):
                    final_captain_id = header

                # Update database
                success = self.update_database_tables(
                    cursor,
                    captain_row_id,
                    accession_id,
                    sequence,
                    header,
                    evidence,
                    resolved_ship_id,
                    final_captain_id,
                )

                if success:
                    updated_count += 1
                    print(
                        f"Updated captain {final_captain_id} with sequence from {evidence}"
                    )
                else:
                    failed_count += 1

            # Commit all changes
            conn.commit()
            print(f"Update complete: {updated_count} successful, {failed_count} failed")


def main():
    updater = CaptainSequenceUpdater()
    updater.run_update()


if __name__ == "__main__":
    main()
