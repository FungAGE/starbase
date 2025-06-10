# coding: utf-8
from sqlalchemy import (
    Column,
    ForeignKey,
    Integer,
    String,
    Table,
    Text,
    CheckConstraint,
    DateTime,
    Enum,
    Index,
    Boolean,
)
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from datetime import datetime
import enum

Base = declarative_base()
metadata = Base.metadata


class BaseModel(Base):
    """Abstract base model with common fields"""

    __abstract__ = True
    created_at = Column(
        DateTime, default=datetime.now().strftime("%Y-%m-%d %H:%M:%S"), nullable=True
    )
    updated_at = Column(
        DateTime,
        default=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        onupdate=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        nullable=True,
    )
    deleted_at = Column(DateTime, nullable=True)
    is_deleted = Column(Boolean, default=False)


class Accessions(BaseModel):
    __tablename__ = "accessions"
    __table_args__ = (Index("idx_accession_tag", "accession_tag"),)

    id = Column(Integer, primary_key=True)
    ship_name = Column(String(255), nullable=False)
    accession_tag = Column(String(50), nullable=True, index=True)

    # One-way relationships - remove all back_populates
    ships = relationship("Ships", cascade="all, delete-orphan")
    captains = relationship("Captains", cascade="all, delete-orphan")
    starship_features = relationship("StarshipFeatures")
    gff_entries = relationship("Gff", cascade="all, delete-orphan")
    classifications = relationship("Classification")


class Ships(BaseModel):
    __tablename__ = "ships"
    __table_args__ = (Index("idx_ships_accession", "accession_id"),)
    id = Column(Integer, primary_key=True)
    sequence = Column(Text)
    md5 = Column(String(32))
    sequence_length = Column(Integer)
    evidence = Column(String(50))
    source = Column(String(50))
    curated_status = Column(String(50))
    orphan = Column(String(50))
    accession_id = Column(Integer, ForeignKey("accessions.id"), nullable=False)
    captain_id = Column(Integer, ForeignKey("captains.id"))
    tax_id = Column(Integer, ForeignKey("taxonomy.id"))
    genome_id = Column(Integer, ForeignKey("genomes.id"))
    classification_id = Column(Integer, ForeignKey("classification.id"))


class Captains(BaseModel):
    __tablename__ = "captains"
    id = Column(Integer, primary_key=True)
    captain_name = Column(String(50), unique=True)
    sequence = Column(Text)
    accession_id = Column(Integer, ForeignKey("accessions.id"))
    reviewed = Column(String)
    evidence = Column(String)


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


class BoundaryType(enum.Enum):
    flank = "flank"
    insert = "insert"
    unknown = "unknown"


class StarshipFeatures(BaseModel):
    __tablename__ = "starship_features"
    id = Column(Integer, primary_key=True)
    accession_id = Column(Integer, ForeignKey("accessions.id"))
    captain_id = Column(Integer, ForeignKey("captains.id"))
    contig = Column(String(100), nullable=False)
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
    dr = Column(String)
    tir = Column(String)
    target = Column(String)
    spok = Column(String)
    ars = Column(String)
    other = Column(String)
    hgt = Column(String)


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


class Papers(Base):
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


class FamilyNames(BaseModel):
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


class NavisNames(BaseModel):
    __tablename__ = "navis_names"
    id = Column(Integer, primary_key=True)
    navis_name = Column(String)
    previous_navis_name = Column(String)
    family_id = Column(Integer, ForeignKey("family_names.id"))


class HaplotypeNames(BaseModel):
    __tablename__ = "haplotype_names"
    id = Column(Integer, primary_key=True)
    haplotype_name = Column(String)
    previous_haplotype_name = Column(String)
    family_id = Column(Integer, ForeignKey("family_names.id"))
    navis_id = Column(Integer, ForeignKey("navis_names.id"))


class Taxonomy(BaseModel):
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
    strain = Column(String(50))
    ship = relationship("Ships", back_populates="tax")
    genome = relationship("Genome", back_populates="tax")
    classification = relationship("Classification", back_populates="tax")


class Gff(Base):
    __tablename__ = "gff"
    __table_args__ = (Index("idx_gff_accession_id", "accession_id"),)
    id = Column(Integer, primary_key=True)
    contig = Column(String(100))
    accession_id = Column(Integer, ForeignKey("accessions.id"))
    source = Column(String)
    type = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    phase = Column(Integer)
    strand = Column(String)
    score = Column(String)
    attributes = Column(String)


class Classification(BaseModel):
    __tablename__ = "classification"
    id = Column(Integer, primary_key=True)
    accession_id = Column(Integer, ForeignKey("accessions.id"))
    taxonomy_id = Column(Integer, ForeignKey("taxonomy.id"))
    genome_id = Column(Integer, ForeignKey("genomes.id"))
    family_id = Column(Integer, ForeignKey("family_names.id"))
    navis_id = Column(Integer, ForeignKey("navis_names.id"))
    haplotype_id = Column(Integer, ForeignKey("haplotype_names.id"))

    # Relationships
    accession = relationship("Accessions", back_populates="classification")
    taxonomy = relationship("Taxonomy", back_populates="classification")
    genome = relationship("Genome", back_populates="classification")
    family = relationship("FamilyNames", back_populates="classification")
    navis = relationship("NavisNames", back_populates="classification")
    haplotype = relationship("HaplotypeNames", back_populates="classification")
