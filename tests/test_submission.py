import pytest
from datetime import datetime
from models import Submission
from src.config.database import StarbaseSession, SubmissionsSession
import os
from src.utils.classification_utils import assign_accession
from src.database.sql_manager import fetch_ships

def test_create_new_submission():
    """Test creating and saving a new submission record"""
    # Create a session from the sessionmaker
    session = SubmissionsSession()
    
    try:
        # Load test sequence from file
        test_data_path = "/home/adrian/Downloads/Hephaestus_all.fasta.split/Hephaestus_all.part_Aspfis_NFIA_048490.fasta"
        sequence_filename = os.path.basename(test_data_path)
        # test_data_path = os.path.join(os.path.dirname(__file__), "test_data", "test_sequence.fasta")
        with open(test_data_path, "r") as f:
            sequence_content = f.read()

        # Create mock existing ships data for testing accession assignment
        existing_ships = fetch_ships(curated=True)

        # Get accession and review status
        accession, needs_review = assign_accession(sequence_content, existing_ships)
        
        assert accession.startswith('SBS')
        assert needs_review is True

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
            needs_review=True,
        )

        # Add and commit the new submission to the database
        session.add(new_submission)
        session.commit()

        # Verify the submission was created
        saved_submission = (
            session.query(Submission)
            .filter_by(seq_filename=sequence_filename)
            .first()
        )
        
        assert saved_submission is not None
        assert saved_submission.seq_contents == sequence_content
        assert saved_submission.seq_filename == sequence_filename
        # Just verify the date matches today
        assert saved_submission.seq_date.split()[0] == today
        assert saved_submission.uploader == "adrian"
        assert saved_submission.evidence == "manual"
        assert saved_submission.genus == "Aspergillus"
        assert saved_submission.species == "fischeri"
        assert saved_submission.hostchr == "unknown"
        assert saved_submission.shipstart == 1
        assert saved_submission.shipend == 36857
        assert saved_submission.shipstrand == "+"
        assert saved_submission.comment == "Test submission"
        assert saved_submission.needs_review is True

    finally:
        # Cleanup - roll back the transaction and close the session
        session.rollback()
        session.close()