from Bio import SeqIO
import hashlib

from sqlalchemy import (
    create_engine,
    text,
    Column,
    Integer,
    String,
    Numeric,
    VARCHAR,
    ForeignKey,
)
from sqlalchemy.orm import sessionmaker, relationship

from src.components.sql_engine import engine, Base


# set up table classes
class Accessions(Base):
    __tablename__ = "accessions"
    id = Column(Integer, primary_key=True, autoincrement=True)
    ship_name = Column(String)
    accession = Column(String)
    accession_tag = Column(String)
    accession_new = Column(Numeric)

    ships = relationship("Ships", order_by="Ships.id", back_populates="accession_obj")
    gff = relationship("Gff", order_by="Gff.ship_id", back_populates="accession_obj")
    joined_ships = relationship(
        "JoinedShips", order_by="JoinedShips.ship_id", back_populates="accession_obj"
    )


class Ships(Base):
    __tablename__ = "ships"
    id = Column(Integer, primary_key=True, autoincrement=True)
    sequence = Column(String)
    md5 = Column(String)
    accession = Column(Integer, ForeignKey("accessions.id"))

    accession_obj = relationship("Accessions", back_populates="ships")


class Captains(Base):
    __tablename__ = "captains"
    id = Column(Integer, primary_key=True, autoincrement=True)
    captainID = Column(String)
    sequence = Column(String)
    ship_id = Column(Integer)
    reviewed = Column(String)


class Genomes(Base):
    __tablename__ = "genomes"
    id = Column(Integer)
    ome = Column(VARCHAR)
    genus = Column(VARCHAR)
    species = Column(VARCHAR)
    strain = Column(VARCHAR)
    version = Column(VARCHAR)
    genomeSource = Column(VARCHAR)
    citation = Column(VARCHAR)
    biosample = Column(VARCHAR)
    acquisition_date = Column(Integer)


class Taxonomy(Base):
    __tablename__ = "taxonomy"
    id = Column(Integer)
    name = Column(VARCHAR)
    taxID = Column(VARCHAR)
    superkingdom = Column(VARCHAR)
    clade = Column(VARCHAR)
    kingdom = Column(VARCHAR)
    subkingdom = Column(VARCHAR)
    phylum = Column(VARCHAR)
    subphylum = Column(VARCHAR)
    class_ = Column(VARCHAR, name="class")
    subclass = Column(VARCHAR)
    order = Column(VARCHAR)
    suborder = Column(VARCHAR)
    family = Column(VARCHAR)
    genus = Column(VARCHAR)
    species = Column(VARCHAR)
    section = Column(VARCHAR)


class FamilyNames(Base):
    __tablename__ = "family_names"
    id = Column(Integer)
    longFamilyID = Column(VARCHAR)
    oldFamilyID = Column(VARCHAR)
    clade = Column(Integer)
    newFamilyID = Column(Integer)
    familyName = Column(VARCHAR)
    type_element_reference = Column(VARCHAR)
    notes = Column(VARCHAR)
    otherFamilyID = Column(VARCHAR)
    paper_id = Column(Integer)


class Papers(Base):
    __tablename__ = "papers"
    id = Column(Integer)
    Key = Column(String)
    ItemType = Column(String)
    PublicationYear = Column(Integer)
    Author = Column(String)
    Title = Column(String)
    PublicationTitle = Column(String)
    DOI = Column(String)
    Url = Column(String)
    AbstractNote = Column(String)
    Date = Column(String)
    starshipMentioned = Column(String)
    typePaper = Column(String)
    shortCitation = Column(VARCHAR)


class NavisHaplotype(Base):
    __tablename__ = "navis_haplotype"
    id = Column(Integer)
    starship_navis = Column(String)
    starship_haplotype = Column(String)
    ship_family_id = Column(Integer)


class Gff(Base):
    __tablename__ = "gff"
    accession = Column(VARCHAR)
    source = Column(VARCHAR)
    type = Column(VARCHAR)
    start = Column(Integer)
    end = Column(Integer)
    phase = Column(VARCHAR)
    strand = Column(VARCHAR)
    score = Column(VARCHAR)
    attributes = Column(VARCHAR)
    ship_id = Column(Integer)


class StarshipFeatures(Base):
    __tablename__ = "starship_features"
    id = Column(Integer)
    contigID = Column(VARCHAR)
    starshipID = Column(VARCHAR)
    captainID = Column(VARCHAR)
    elementBegin = Column(VARCHAR)
    elementEnd = Column(VARCHAR)
    elementLength = Column(VARCHAR)
    strand = Column(VARCHAR)
    boundaryType = Column(VARCHAR)
    emptySiteID = Column(VARCHAR)
    emptyContig = Column(VARCHAR)
    emptyBegin = Column(VARCHAR)
    emptyEnd = Column(VARCHAR)
    emptySeq = Column(VARCHAR)
    upDR = Column(VARCHAR)
    downDR = Column(VARCHAR)
    DRedit = Column(VARCHAR)
    upTIR = Column(VARCHAR)
    downTIR = Column(VARCHAR)
    TIRedit = Column(VARCHAR)
    nestedInside = Column(VARCHAR)
    containNested = Column(VARCHAR)
    ship_id = Column(Integer)


class JoinedShips(Base):
    __tablename__ = "joined_ships"
    starshipID = Column(String)
    genus = Column(String)
    species = Column(String)
    strain = Column(String)
    evidence = Column(String)
    source = Column(String)
    contigID = Column(String)
    captainID = Column(String)
    elementBegin = Column(Integer)
    elementEnd = Column(Integer)
    size = Column(Integer)
    strand = Column(String)
    boundaryType = Column(String)
    emptySiteID = Column(String)
    emptyContig = Column(String)
    emptyBegin = Column(Integer)
    emptyEnd = Column(Integer)
    emptySeq = Column(String)
    upDR = Column(String)
    downDR = Column(String)
    DRedit = Column(String)
    upTIR = Column(String)
    downTIR = Column(String)
    TIRedit = Column(String)
    nestedInside = Column(String)
    containNested = Column(String)
    dr = Column(String)
    tir = Column(String)
    starship_navis = Column(String)
    starship_haplotype = Column(String)
    target = Column(String)
    spok = Column(String)
    ars = Column(String)
    other = Column(String)
    hgt = Column(String)
    ship_family_id = Column(Integer)
    curated_status = Column(String)
    taxid = Column(Integer)
    ship_id = Column(Integer)
    genome_id = Column(String)
    ome = Column(String)
    orphan = Column(String)
    captainID_new = Column(Integer)


class Submissions(Base):
    __tablename__ = "submissions_table"
    seq_contents = Column(String)
    seq_filename = Column(String)
    seq_date = Column(String)
    anno_contents = Column(String)
    anno_filename = Column(String)
    anno_date = Column(String)
    uploader = Column(String)
    evidence = Column(String)
    genus = Column(String)
    species = Column(String)
    hostchr = Column(String)
    shipstart = Column(Integer)
    shipend = Column(Integer)
    shipstrand = Column(String)
    comment = Column(String)
    id = Column(Integer)


# TODO: add new tables
# - empty sites
# - representative ships (maybe just add a column?)
# - blast results

# relationships

Ships.accession = relationship(
    "Accessions", order_by=Accessions.id, back_populates="ships"
)

Gff.ship_id = relationship("Accessions", order_by=Accessions.id, back_populates="gff")

# Initialize the SQLite engine
Session = sessionmaker(bind=engine)


def update_md5(engine, table, id_column, seq_column, md5_column):
    # Open a connection to the database
    with engine.connect() as connection:
        # Fetch rows where the md5 column is still NULL
        select_query = text(
            f"SELECT {id_column}, {seq_column} FROM {table} WHERE {md5_column} IS NULL"
        )
        result = connection.execute(select_query).fetchall()

        # Loop through the rows and calculate the MD5 hash for each long string
        for row in result:
            id_value = row[0]
            seq = row[1]

            # Generate the MD5 hash
            md5_hash = hashlib.md5(seq.encode()).hexdigest()

            # Update the table with the generated MD5 hash
            update_query = text(
                f"UPDATE {table} SET {md5_column} = :md5_hash WHERE {id_column} = :id_value"
            )
            connection.execute(
                update_query, {"md5_hash": md5_hash, "id_value": id_value}
            )


def update_table(engine, table, id_column, seq_column, fasta_file):
    records = SeqIO.parse(fasta_file, "fasta")
    with engine.connect() as connection:
        with connection.begin():
            for record in records:
                name = str(record["id"])
                sequence = str(record["seq"])
                insert_query = text(
                    f"""
                    INSERT INTO {table} ({id_column}, {seq_column})
                    VALUES (:name, :sequence);
                    """
                )
                connection.execute(insert_query, {"name": name, "sequence": sequence})


update_table(
    engine,
    "ships",
)


# TODO: update gff
# TODO: verify attributes

# TODO: update features or joined_ships?

# TODO: update genomes

# TODO: update taxonomy

# TODO: check for duplicates
