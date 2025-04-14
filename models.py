# coding: utf-8
from sqlalchemy import Column, Float, ForeignKey, Integer, Numeric, String, Table, Text
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata


class Accession(Base):
    __tablename__ = "accessions"

    id = Column(Integer, primary_key=True)
    ship_name = Column(Text)
    accession = Column(Text)
    accession_tag = Column(String)
    accession_new = Column(Numeric)

    captains = relationship("Captain", back_populates="ship")


class Captain(Base):
    __tablename__ = "captains"

    id = Column(Integer, primary_key=True)
    captainID = Column(Text)
    sequence = Column(Text)
    ship_id = Column(Integer, ForeignKey("accessions.id"))
    reviewed = Column(Text)

    ship = relationship("Accession", back_populates="captains")


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


class Paper(Base):
    __tablename__ = "papers"

    id = Column(Integer, primary_key=True)
    Key = Column(Text)
    ItemType = Column(Text)
    PublicationYear = Column(Integer)
    Author = Column(Text)
    Title = Column(Text)
    PublicationTitle = Column(Text)
    DOI = Column(Text)
    Url = Column(Text)
    AbstractNote = Column(Text)
    Date = Column(Text)
    starshipMentioned = Column(Text)
    typePaper = Column(Text)
    shortCitation = Column(String)

    family_names = relationship("FamilyName", back_populates="paper")


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
    _class = Column("class", String)
    subclass = Column(String)
    order = Column(String)
    suborder = Column(String)
    family = Column(String)
    genus = Column(String)
    species = Column(String)
    section = Column(String)
    species_group = Column("species group", String)
    subgenus = Column(String)
    strain = Column(String)

    genomes = relationship("Genome", back_populates="taxonomy")


class FamilyName(Base):
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

    paper = relationship("Paper", back_populates="family_names")


t_gff = Table(
    "gff",
    metadata,
    Column("accession", String),
    Column("source", String),
    Column("type", String),
    Column("start", Integer),
    Column("end", Integer),
    Column("phase", String),
    Column("strand", String),
    Column("score", String),
    Column("attributes", String),
    Column("ship_id", ForeignKey("accessions.id")),
)


t_ships = Table(
    "ships",
    metadata,
    Column("id", Integer),
    Column("sequence", Text),
    Column("md5", Text),
    Column("accession", ForeignKey("accessions.id")),
)


t_starship_features = Table(
    "starship_features",
    metadata,
    Column("id", Integer, primary_key=True),
    Column("contigID", String),
    Column("starshipID", String),
    Column("captainID", String),
    Column("elementBegin", String),
    Column("elementEnd", String),
    Column("elementLength", String),
    Column("strand", String),
    Column("boundaryType", String),
    Column("emptySiteID", String),
    Column("emptyContig", String),
    Column("emptyBegin", String),
    Column("emptyEnd", String),
    Column("emptySeq", String),
    Column("upDR", String),
    Column("downDR", String),
    Column("DRedit", String),
    Column("upTIR", String),
    Column("downTIR", String),
    Column("TIRedit", String),
    Column("nestedInside", String),
    Column("containNested", String),
    Column("ship_id", Integer, ForeignKey("accessions.id")),
    Column("captain_id", Integer, ForeignKey("captains.id")),
)


t_joined_ships = Table(
    "joined_ships",
    metadata,
    Column("starshipID", Text),
    Column("genus", Text),
    Column("species", Text),
    Column("strain", Text),
    Column("evidence", Text),
    Column("source", Text),
    Column("contigID", Text),
    Column("captainID", Text),
    Column("elementBegin", Float),
    Column("elementEnd", Float),
    Column("size", Float),
    Column("strand", Text),
    Column("boundaryType", Text),
    Column("emptySiteID", Text),
    Column("emptyContig", Text),
    Column("emptyBegin", Float),
    Column("emptyEnd", Float),
    Column("emptySeq", Text),
    Column("upDR", Text),
    Column("downDR", Text),
    Column("DRedit", Text),
    Column("upTIR", Text),
    Column("downTIR", Text),
    Column("TIRedit", Text),
    Column("nestedInside", Text),
    Column("containNested", Text),
    Column("dr", Text),
    Column("tir", Text),
    Column("starship_navis", Text),
    Column("starship_haplotype", Text),
    Column("target", Text),
    Column("spok", Text),
    Column("ars", Text),
    Column("other", Text),
    Column("hgt", Text),
    Column("ship_family_id", ForeignKey("family_names.id")),
    Column("curated_status", Text),
    Column("taxid", ForeignKey("taxonomy.id")),
    Column("ship_id", ForeignKey("accessions.id")),
    Column("genome_id", ForeignKey("genomes.id")),
    Column("ome", Text),
    Column("orphan", Text),
    Column("captainID_new", Integer),
)


t_navis_haplotype = Table(
    "navis_haplotype",
    metadata,
    Column("id", Integer, primary_key=True),
    Column("starship_navis", Text),
    Column("starship_haplotype", Text),
    Column("ship_family_id", Integer, ForeignKey("family_names.id")),
)
