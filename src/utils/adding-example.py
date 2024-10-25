from models import Accession, Captain

# Create a new accession record
new_accession = Accession(
    ship_name="Enterprise",
    accession="123456",
    accession_tag="SBS123456",
    accession_new=1,
)

# Add and commit the new accession to the database
session.add(new_accession)
session.commit()
