# coding: utf-8
from sqlalchemy import (
    Column,
    ForeignKey,
    Integer,
    String,
    Table,
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

    # Relationships
    ships = relationship("Ships", back_populates="accession_obj")
    gff_entries = relationship("Gff", back_populates="accession_obj")
    joined_ships = relationship("JoinedShips", back_populates="accession_obj")
    captains = relationship("Captains", back_populates="ship")


class Ships(Base):
    __tablename__ = "ships"
    id = Column(Integer, primary_key=True)
    sequence = Column(String)
    md5 = Column(String)
    accession_id = Column(Integer, ForeignKey("accessions.id"))

    # Relationships
    accession_obj = relationship("Accessions", back_populates="ships")


class Captains(Base):
    __tablename__ = "captains"
    id = Column(Integer, primary_key=True)
    captainID = Column(String, unique=True)
    sequence = Column(String)
    ship_id = Column(Integer, ForeignKey("accessions.id"))
    reviewed = Column(String)

    # Relationships
    ship = relationship("Accessions", back_populates="captains")
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

    taxonomy = relationship("Taxonomy", back_populates="genomes")


class StarshipFeatures(Base):
    __tablename__ = "starship_features"
    id = Column(Integer, primary_key=True)
    contigID = Column(String)
    starshipID = Column(String)
    captainID = Column(String)
    elementBegin = Column(String)
    elementEnd = Column(String)
    elementLength = Column(String)
    strand = Column(String)
    boundaryType = Column(String)
    emptySiteID = Column(String)
    emptyContig = Column(String)
    emptyBegin = Column(String)
    emptyEnd = Column(String)
    emptySeq = Column(String)
    upDR = Column(String)
    downDR = Column(String)
    DRedit = Column(String)
    upTIR = Column(String)
    downTIR = Column(String)
    TIRedit = Column(String)
    nestedInside = Column(String)
    containNested = Column(String)
    ship_id = Column(Integer, ForeignKey("accessions.id"))
    captain_id = Column(Integer, ForeignKey("captains.id"))

    # Relationships
    accession = relationship("Accessions", back_populates="starship_features")
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
    shortCitation = Column(VARCHAR)

    # Relationships
    family_names = relationship(
        "FamilyNames", secondary=paper_family_association, back_populates="papers"
    )


class FamilyNames(Base):
    __tablename__ = "family_names"
    id = Column(Integer, primary_key=True)
    longFamilyID = Column(VARCHAR)
    oldFamilyID = Column(VARCHAR)
    clade = Column(Integer)
    newFamilyID = Column(Integer)
    familyName = Column(VARCHAR)
    type_element_reference = Column(VARCHAR)
    notes = Column(VARCHAR)
    otherFamilyID = Column(VARCHAR)
    paper_id = Column(Integer, ForeignKey("papers.id"))

    # Relationships
    papers = relationship(
        "Papers", secondary=paper_family_association, back_populates="family_names"
    )


class Taxonomy(Base):
    __tablename__ = "taxonomy"
    id = Column(Integer, primary_key=True)
    name = Column(VARCHAR)
    taxID = Column(VARCHAR)
    superkingdom = Column(VARCHAR)
    clade = Column(VARCHAR)
    kingdom = Column(VARCHAR)
    subkingdom = Column(VARCHAR)
    phylum = Column(VARCHAR)
    subphylum = Column(VARCHAR)
    class_ = Column(
        VARCHAR, name="class"
    )  # Using class_ as class is a reserved keyword
    subclass = Column(VARCHAR)
    order = Column(VARCHAR)
    suborder = Column(VARCHAR)
    family = Column(VARCHAR)
    genus = Column(VARCHAR)
    species = Column(VARCHAR)
    section = Column(VARCHAR)


class NavisHaplotype(Base):
    __tablename__ = "navis_haplotype"
    id = Column(Integer, primary_key=True)
    navis_name = Column(String)
    haplotype_name = Column(String)
    ship_family_id = Column(Integer, ForeignKey("family_names.id"))

    # Relationships
    family = relationship("FamilyNames")


class Gff(Base):
    __tablename__ = "gff"
    id = Column(Integer, primary_key=True)
    contigID = Column(String)
    accession = Column(String)
    source = Column(String)
    type = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    phase = Column(Integer)
    strand = Column(String)
    score = Column(String)
    attributes = Column(String)
    ship_id = Column(Integer, ForeignKey("accessions.id"))


class JoinedShips(Base):
    __tablename__ = "joined_ships"
    id = Column(Integer, primary_key=True)
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
    navis_name = Column(String)
    haplotype_name = Column(String)
    target = Column(String)
    spok = Column(String)
    ars = Column(String)
    other = Column(String)
    hgt = Column(String)
    ship_family_id = Column(Integer, ForeignKey("family_names.id"))
    curated_status = Column(String)
    taxid = Column(Integer)
    ship_id = Column(Integer, ForeignKey("accessions.id"))
    genome_id = Column(String)
    ome = Column(String)
    orphan = Column(String)
    captainID_new = Column(Integer)

    # Relationships
    accession_obj = relationship("Accessions", back_populates="joined_ships")
    family = relationship("FamilyNames")
