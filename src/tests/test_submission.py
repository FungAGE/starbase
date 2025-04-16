from models import Accession


def test_create_new_accession(db_session):
    """Test creating and saving a new accession record"""
    # Create a new accession record
    new_accession = Accession(
        ship_name="Enterprise",
        accession="123456",
        accession_tag="SBS123456",
        accession_new=1,
    )

    # Add and commit the new accession to the database
    db_session.add(new_accession)
    db_session.commit()

    # Verify the accession was created
    saved_accession = db_session.query(Accession).filter_by(accession="123456").first()
    assert saved_accession is not None
    assert saved_accession.ship_name == "Enterprise"
    assert saved_accession.accession_tag == "SBS123456"
