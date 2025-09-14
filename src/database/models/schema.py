# coding: utf-8
from sqlalchemy import (
    Column,
    ForeignKey,
    Integer,
    String,
    Table,
    DateTime,
    Boolean,
    VARCHAR,
)
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata

# Association tables for many-to-many relationships
paper_family_association = Table(
    "paper_family_association",
    Base.metadata,
    Column("paper_id", Integer, ForeignKey("papers.id")),
    Column("family_id", Integer, ForeignKey("family_names.id")),
)


class Accessions(Base):
    __tablename__ = "accessions"
    id = Column(Integer, primary_key=True)
    ship_name = Column(String)
    accession_tag = Column(String, unique=True)
    version_tag = Column(Integer)

    # Relationships
    ships = relationship("Ships", back_populates="accession_obj")
    gff_entries = relationship("Gff", back_populates="accession_obj")


class Ships(Base):
    __tablename__ = "ships"
    id = Column(Integer, primary_key=True)
    sequence = Column(String)
    md5 = Column(String)
    rev_comp_md5 = Column(String)
    accession_id = Column(Integer, ForeignKey("accessions.id"))

    # Relationships
    accession_obj = relationship("Accessions", back_populates="ships")


class Captains(Base):
    __tablename__ = "captains"
    id = Column(Integer, primary_key=True)
    captainID = Column(String, unique=True)
    sequence = Column(String)
    ship_id = Column(Integer, ForeignKey("ships.id"))
    reviewed = Column(String)

    # Relationships
    ship = relationship("Ships")
    features = relationship("StarshipFeatures", back_populates="captain")


class Genome(Base):
    __tablename__ = "genomes"

    id = Column(Integer, primary_key=True)
    ome = Column(String(50))
    taxonomy_id = Column(Integer, ForeignKey("taxonomy.id"))
    version = Column(String(50))
    genomeSource = Column(String(50))
    citation = Column(String(50))
    biosample = Column(String(50))
    acquisition_date = Column(Integer)
    assembly_accession = Column(String(50))
    taxonomy = relationship("Taxonomy", back_populates="genomes")


class StarshipFeatures(Base):
    __tablename__ = "starship_features"
    id = Column(Integer, primary_key=True)
    contigID = Column(String)
    starshipID = Column(String)
    captainID = Column(String)
    elementBegin = Column(Integer)
    elementEnd = Column(Integer)
    elementLength = Column(Integer)
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
    ship_id = Column(Integer, ForeignKey("ships.id"))
    captain_id = Column(Integer, ForeignKey("captains.id"))

    # Relationships
    captain = relationship("Captains", back_populates="features")


class Papers(Base):
    __tablename__ = "papers"
    id = Column(Integer, primary_key=True)
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
    shortCitation = Column(String)

    # Relationships
    family_names = relationship(
        "FamilyNames", secondary=paper_family_association, back_populates="papers"
    )


class FamilyNames(Base):
    __tablename__ = "family_names"
    id = Column(Integer, primary_key=True)
    longFamilyID = Column(String)
    oldFamilyID = Column(String)
    clade = Column(Integer)
    newFamilyID = Column(Integer)
    familyName = Column(String)
    type_element_reference = Column(String)
    notes = Column(String)
    otherFamilyID = Column(String)
    paper_id = Column(Integer, ForeignKey("papers.id"))

    # Relationships
    papers = relationship(
        "Papers", secondary=paper_family_association, back_populates="family_names"
    )


class Taxonomy(Base):
    __tablename__ = "taxonomy"
    id = Column(Integer, primary_key=True)
    name = Column(String)
    taxID = Column(String)
    superkingdom = Column(String)
    clade = Column(String)
    kingdom = Column(String)
    subkingdom = Column(String)
    phylum = Column(String)
    subphylum = Column(String)
    class_ = Column(
        String, name="class"
    )  # Using class_ as class is a reserved keyword
    subclass = Column(String)
    order = Column(String)
    suborder = Column(String)
    family = Column(String)
    genus = Column(String)
    species = Column(String)
    section = Column(String)
    species_group = Column(String)
    subgenus = Column(String)
    strain = Column(String)

    # Relationships
    genomes = relationship("Genome", back_populates="taxonomy")


class Navis(Base):
    __tablename__ = "navis_names"
    id = Column(Integer, primary_key=True)
    navis_name = Column(String)
    previous_navis_name = Column(String)
    ship_family_id = Column(Integer, ForeignKey("family_names.id"))

    # Relationships
    family = relationship("FamilyNames")

class Haplotype(Base):
    __tablename__ = "haplotype_names"
    id = Column(Integer, primary_key=True)
    haplotype_name = Column(String)
    previous_haplotype_name = Column(String)
    navis_id = Column(Integer, ForeignKey("navis_names.id"))
    ship_family_id = Column(Integer, ForeignKey("family_names.id"))

    # Relationships
    family = relationship("FamilyNames")
    navis = relationship("Navis")


class Gff(Base):
    __tablename__ = "gff"
    id = Column(Integer, primary_key=True)
    accession_id = Column(Integer, ForeignKey("accessions.id"))
    source = Column(String)
    type = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    phase = Column(Integer)
    strand = Column(String)
    score = Column(String)
    attributes = Column(String)
    ship_id = Column(Integer, ForeignKey("ships.id"))

    # Relationships
    accession_obj = relationship("Accessions", back_populates="gff_entries")


class JoinedShips(Base):
    __tablename__ = "joined_ships"
    id = Column(Integer, primary_key=True)
    starshipID = Column(String)
    evidence = Column(String)
    source = Column(String)
    curated_status = Column(String)
    ship_family_id = Column(Integer, ForeignKey("family_names.id"))
    tax_id = Column(Integer, ForeignKey("taxonomy.id"))
    ship_id = Column(Integer, ForeignKey("ships.id"))
    genome_id = Column(Integer, ForeignKey("genomes.id"))
    captain_id = Column(Integer, ForeignKey("captains.id"))
    ship_navis_id = Column(Integer, ForeignKey("navis_names.id"))
    ship_haplotype_id = Column(Integer, ForeignKey("haplotype_names.id"))
    created_at = Column(DateTime)
    updated_at = Column(DateTime)

    # Relationships
    family = relationship("FamilyNames")
    taxonomy = relationship("Taxonomy")
    ship = relationship("Ships")
    genome = relationship("Genome")
    captain = relationship("Captains")
    navis = relationship("Navis")
    haplotype = relationship("Haplotype")

class Submission(Base):
    __tablename__ = "submissions"
    seq_contents = Column(String)
    seq_filename = Column(String)
    seq_date = Column(DateTime)
    anno_contents = Column(String)
    anno_filename = Column(String)
    anno_date = Column(DateTime)
    uploader = Column(String)
    evidence = Column(String)
    genus = Column(String)
    species = Column(String)
    hostchr = Column(String)
    shipstart = Column(Integer)
    shipend = Column(Integer)
    shipstrand = Column(String)
    comment = Column(String)
    id = Column(Integer, primary_key=True)
    accession_tag = Column(VARCHAR(255))
    needs_review = Column(Boolean)