from datetime import datetime
from src.database.models.schema import Submission
from src.config.database import SubmissionsSession
import os
from src.utils.classification_utils import assign_accession
from src.database.sql_manager import fetch_ships


def test_create_new_submission():
    """Test creating and saving a new submission record"""
    # Create a session from the sessionmaker
    session = SubmissionsSession()

    try:
        # Load test sequence from test data file
        test_data_path = os.path.join(os.path.dirname(__file__), "test_data", "test_sequence.fasta")
        sequence_filename = os.path.basename(test_data_path)
        
        with open(test_data_path, "r") as f:
            sequence_content = f.read()

        # Create mock existing ships data for testing accession assignment
        import pandas as pd
        import hashlib
        existing_ships = pd.DataFrame({
            "accession_tag": ["SBS000001", "SBS000002"],
            "sequence": [
                "ATGCATGCATGCATGCATGC",
                "TTTTTTTTTTTTTTTTTTTT"
            ],
            "md5": [
                hashlib.md5("ATGCATGCATGCATGCATGC".encode()).hexdigest(),
                hashlib.md5("TTTTTTTTTTTTTTTTTTTT".encode()).hexdigest(),
            ],
            "rev_comp_md5": [
                hashlib.md5("GCATGCATGCATGCATGCAT".encode()).hexdigest(),
                hashlib.md5("AAAAAAAAAAAAAAAAAAAA".encode()).hexdigest(),
            ]
        })

        # Get accession and review status
        accession, needs_review = assign_accession(sequence_content, existing_ships)

        assert accession.startswith("SBS")
        assert needs_review is False  # Should be new since it's different from existing

        today = datetime.now().strftime("%Y-%m-%d")

        # Create a new submission record
        new_submission = Submission(
            seq_contents=sequence_content,
            seq_filename=sequence_filename,
            seq_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            uploader="adrian",
            evidence="manual",
            genus="Aspergillus",
            species="fischeri",
            hostchr="unknown",
            shipstart=1,
            shipend=36857,
            shipstrand="+",
            comment="Test submission",
            needs_review=False,
        )

        # Add and commit the new submission to the database
        session.add(new_submission)
        session.commit()

        # Verify the submission was created
        assert new_submission.id is not None
        assert new_submission.seq_filename == sequence_filename
        assert new_submission.uploader == "adrian"

        print(f"Submission created with ID: {new_submission.id}")

    finally:
        session.close()
