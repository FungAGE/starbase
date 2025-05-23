# coding: utf-8
from sqlalchemy import (
    Column,
    ForeignKey,
    Integer,
    Numeric,
    String,
    Table,
    Text,
    CheckConstraint,
    DateTime,
    Enum,
    Index,
)
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from datetime import datetime
import enum

Base = declarative_base()
metadata = Base.metadata

# Association tables for many-to-many relationships
paper_family_association = Table(
    "paper_family_association",
    Base.metadata,
    Column("id", Integer, primary_key=True),
    Column("paper_id", Integer, ForeignKey("papers.id", ondelete="CASCADE")),
    Column("family_id", Integer, ForeignKey("family_names.id", ondelete="CASCADE")),
    Column(
        "created_at", DateTime, default=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    ),
)


class BaseModel(Base):
    """Abstract base model with common fields"""

    __abstract__ = True
    created_at = Column(
        DateTime, default=datetime.now().strftime("%Y-%m-%d %H:%M:%S"), nullable=False
    )
    updated_at = Column(
        DateTime,
        default=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        onupdate=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        nullable=False,
    )


class Accessions(BaseModel):
    __tablename__ = "accessions"
    __table_args__ = (
        Index("idx_accessions_tag", "accession_tag"),
        Index("idx_accessions_id", "id"),
    )
    id = Column(Integer, primary_key=True)
    ship_name = Column(String, nullable=False)
    accession = Column(String, unique=True, nullable=False)
    accession_tag = Column(String, nullable=True)
    accession_new = Column(Numeric, default=0)

    # Relationships
    ships = relationship("Ships", back_populates="accession_obj")
    gff_entries = relationship("Gff", back_populates="accession_obj")
    joined_ships = relationship("JoinedShips", back_populates="accession_obj")
    captains = relationship("Captains", back_populates="ship")


class Ships(Base):
    __tablename__ = "ships"
    __table_args__ = (Index("idx_ships_accession", "accession_id"),)
    id = Column(Integer, primary_key=True)
    sequence = Column(Text)
    md5 = Column(String(32))
    accession_id = Column(Integer, ForeignKey("accessions.id"))

    # Relationships
    accession_obj = relationship("Accessions", back_populates="ships")


class Captains(Base):
    __tablename__ = "captains"
    id = Column(Integer, primary_key=True)
    captainID = Column(String, unique=True)
    sequence = Column(Text)
    ship_id = Column(Integer, ForeignKey("accessions.id"))
    reviewed = Column(String)

    # Relationships
    ship = relationship("Accessions", back_populates="captains")
    features = relationship("StarshipFeatures", back_populates="captain")


class Genome(BaseModel):
    __tablename__ = "genomes"
    __table_args__ = (
        Index("idx_genome_taxonomy", "taxonomy_id"),
        Index("idx_genome_version_source", "version", "genomeSource"),
    )

    id = Column(Integer, primary_key=True)
    ome = Column(String(50))
    taxonomy_id = Column(Integer, ForeignKey("taxonomy.id"))
    version = Column(String(50))
    genomeSource = Column(String(50))
    citation = Column(String(50))
    biosample = Column(String(50))
    acquisition_date = Column(Integer)

    taxonomy = relationship("Taxonomy", back_populates="genomes")


class BoundaryType(enum.Enum):
    flank = "flank"
    insert = "insert"
    unknown = "unknown"


class StarshipFeatures(BaseModel):
    __tablename__ = "starship_features"
    id = Column(Integer, primary_key=True)
    contigID = Column(String(100), nullable=False)
    starshipID = Column(String, index=True)
    captainID = Column(String, index=True)
    elementBegin = Column(Integer)
    elementEnd = Column(Integer)
    elementLength = Column(Integer)
    strand = Column(String(1), CheckConstraint("strand IN ('+', '-')"))
    boundaryType = Column(
        Enum(BoundaryType), nullable=False, default=BoundaryType.unknown
    )
    emptySiteID = Column(String(100))
    emptyContig = Column(String(100))
    emptyBegin = Column(Integer)
    emptyEnd = Column(Integer)
    emptySeq = Column(Text)
    upDR = Column(Text)
    downDR = Column(Text)
    DRedit = Column(String(50))
    upTIR = Column(Text)
    downTIR = Column(Text)
    TIRedit = Column(String(50))
    nestedInside = Column(String(100))
    containNested = Column(String(100))
    ship_id = Column(Integer, ForeignKey("accessions.id"))
    captain_id = Column(Integer, ForeignKey("captains.id"))

    # Relationships
    accession = relationship("Accessions", back_populates="starship_features")
    captain = relationship("Captains", back_populates="features")


class Papers(BaseModel):
    __tablename__ = "papers"
    __table_args__ = (
        CheckConstraint("PublicationYear >= 1900"),
        CheckConstraint("DOI REGEXP '^10\\.\\d{4,9}/[-._;()/:\\w]+$' OR DOI IS NULL"),
    )
    id = Column(Integer, primary_key=True)
    Key = Column(String)
    ItemType = Column(String)
    PublicationYear = Column(Integer)
    Author = Column(String(255))
    Title = Column(String)
    PublicationTitle = Column(String)
    DOI = Column(String)
    Url = Column(String)
    AbstractNote = Column(Text)
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
    __table_args__ = (
        Index("idx_family_names_reference", "type_element_reference"),
        Index("idx_family_names_family", "familyName"),
    )
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


class NavisNames(Base):
    __tablename__ = "navis_names"
    id = Column(Integer, primary_key=True)
    navis_name = Column(String)
    previous_navis_name = Column(String)
    ship_family_id = Column(Integer, ForeignKey("family_names.id"))

    # Relationships
    family = relationship("FamilyNames")


class HaplotypeNames(Base):
    __tablename__ = "haplotype_names"
    id = Column(Integer, primary_key=True)
    haplotype_name = Column(String)
    previous_haplotype_name = Column(String)
    ship_family_id = Column(Integer, ForeignKey("family_names.id"))
    ship_navis_id = Column(Integer, ForeignKey("navis_names.id"))

    # Relationships
    family = relationship("FamilyNames")
    navis = relationship("NavisNames")


class Taxonomy(Base):
    __tablename__ = "taxonomy"
    __table_args__ = (
        Index("idx_taxonomy_name", "name"),
        Index("idx_taxonomy_taxid", "taxID"),
    )
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    taxID = Column(String(20))
    superkingdom = Column(String(50))
    clade = Column(String(50))
    kingdom = Column(String(50))
    subkingdom = Column(String(50))
    phylum = Column(String(50))
    subphylum = Column(String(50))
    class_ = Column(
        String(50), name="class"
    )  # Using class_ as class is a reserved keyword
    subclass = Column(String(50))
    order = Column(String(50))
    suborder = Column(String(50))
    family = Column(String(50))
    genus = Column(String(50))
    species = Column(String(50))
    section = Column(String(50))


class Gff(Base):
    __tablename__ = "gff"
    __table_args__ = (Index("idx_gff_ship_id", "ship_id"),)
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


class JoinedShips(BaseModel):
    __tablename__ = "joined_ships"
    __table_args__ = (
        Index("idx_joined_ships_curated", "curated_status"),
        Index("idx_joined_ships_ship_id", "ship_id"),
        Index("idx_joined_ships_taxid", "taxid"),
        Index("idx_joined_ships_ship_family_id", "ship_family_id"),
        Index("idx_joined_ships_genome_id", "genome_id"),
    )
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
    starship_navis = Column(String)
    starship_haplotype = Column(String)
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
