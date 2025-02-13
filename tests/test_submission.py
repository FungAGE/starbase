import pytest
from models import Accession, Captain
from src.config.database import StarbaseSession

def test_create_new_accession():
    """Test creating and saving a new accession record"""
    try:
        # Create a new accession record
        new_accession = Accession(
            ship_name="Enterprise",
            accession="123456",
            accession_tag="SBS123456",
            accession_new=1,
        )

        # Add and commit the new accession to the database
        StarbaseSession.add(new_accession)
        StarbaseSession.commit()

        # Verify the accession was created
        saved_accession = StarbaseSession.query(Accession).filter_by(accession="123456").first()
        assert saved_accession is not None
        assert saved_accession.ship_name == "Enterprise"
        assert saved_accession.accession_tag == "SBS123456"

    finally:
        # Cleanup
        StarbaseSession.rollback()
        StarbaseSession.close()