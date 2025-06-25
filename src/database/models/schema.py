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


# NEW: JoinedShips model - the main table that was missing
class JoinedShips(BaseModel):
    __tablename__ = "joined_ships"
    __table_args__ = (
        Index("idx_joined_ships_starshipid", "starshipID"),
        Index("idx_joined_ships_family", "ship_family_id"),
        Index("idx_joined_ships_navis", "ship_navis_id"),
        Index("idx_joined_ships_haplotype", "ship_haplotype_id"),
        Index("idx_joined_ships_taxonomy", "tax_id"),
        Index("idx_joined_ships_genome", "genome_id"),
    )

    # Primary key - add auto-increment id
    id = Column(Integer, primary_key=True, autoincrement=True)
    starshipID = Column(String, index=True)

    # Data columns
    evidence = Column(String)
    source = Column(String)
    curated_status = Column(String)

    # Foreign key columns
    ship_family_id = Column(Integer, ForeignKey("family_names.id"))
    tax_id = Column(Integer, ForeignKey("taxonomy.id"))
    ship_id = Column(Integer, ForeignKey("ships.id"))
    genome_id = Column(Integer, ForeignKey("genomes.id"))
    captain_id = Column(Integer, ForeignKey("captains.id"))
    ship_navis_id = Column(Integer, ForeignKey("navis_names.id"))
    ship_haplotype_id = Column(Integer, ForeignKey("haplotype_names.id"))

    # Relationships
    family = relationship("FamilyNames", back_populates="joined_ships")
    taxonomy = relationship("Taxonomy", back_populates="joined_ships")
    ship = relationship("Ships", back_populates="joined_ships")
    genome = relationship("Genome", back_populates="joined_ships")
    captain = relationship("Captains", back_populates="joined_ships")
    navis = relationship("NavisNames", back_populates="joined_ships")
    haplotype = relationship("HaplotypeNames", back_populates="joined_ships")


class Accessions(BaseModel):
    __tablename__ = "accessions"
    __table_args__ = (Index("idx_accession_tag", "accession_tag"),)

    id = Column(Integer, primary_key=True)
    ship_name = Column(String(255), nullable=False)
    accession_tag = Column(String(50), nullable=True, index=True)

    # Relationships
    ships = relationship(
        "Ships", back_populates="accession", cascade="all, delete-orphan"
    )
    captains = relationship(
        "Captains", back_populates="accession", cascade="all, delete-orphan"
    )
    starship_features = relationship("StarshipFeatures", back_populates="accession")
    gff_entries = relationship(
        "Gff", back_populates="accession", cascade="all, delete-orphan"
    )
    classifications = relationship("Classification", back_populates="accession")


class Ships(BaseModel):
    __tablename__ = "ships"
    __table_args__ = (Index("idx_ships_accession", "accession_id"),)
    id = Column(Integer, primary_key=True)
    sequence = Column(Text)
    md5 = Column(String(32))
    sequence_length = Column(Integer)
    header = Column(String)
    accession = Column(String)  # Keep for backward compatibility
    accession_id = Column(Integer, ForeignKey("accessions.id"), nullable=False)

    # Relationships
    accession = relationship("Accessions", back_populates="ships")
    joined_ships = relationship("JoinedShips", back_populates="ship")


class Captains(BaseModel):
    __tablename__ = "captains"

    id = Column(Integer, primary_key=True)
    captainID = Column(String(50))  # Match database column name
    captain_name = Column(String(50), unique=True)
    sequence = Column(Text)
    ship_id = Column(Integer, ForeignKey("ships.id"))
    accession_id = Column(Integer, ForeignKey("accessions.id"))
    reviewed = Column(String)
    evidence = Column(String)

    # Relationships
    accession = relationship("Accessions", back_populates="captains")
    joined_ships = relationship("JoinedShips", back_populates="captain")


class Genome(BaseModel):  # Keep class name for compatibility but matches genomes table
    __tablename__ = "genomes"
    __table_args__ = (
        Index("idx_genomes_ome", "ome"),
        Index("idx_genomes_source", "genomeSource"),
    )

    id = Column(Integer, primary_key=True)
    ome = Column(String(50))
    genus = Column(String(50))
    species = Column(String(50))
    strain = Column(String(50))
    version = Column(String(50))
    genomeSource = Column(String(50))
    citation = Column(String(50))
    biosample = Column(String(50))
    acquisition_date = Column(Integer)
    assembly_accession = Column(Text)

    # Relationships
    joined_ships = relationship("JoinedShips", back_populates="genome")
    classifications = relationship("Classification", back_populates="genome")


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

    # Relationships
    accession = relationship("Accessions", back_populates="starship_features")


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

    # Relationships
    papers = relationship(
        "Papers", secondary=paper_family_association, back_populates="family_names"
    )
    navis_names = relationship("NavisNames", back_populates="family")
    haplotype_names = relationship("HaplotypeNames", back_populates="family")
    joined_ships = relationship("JoinedShips", back_populates="family")
    classifications = relationship("Classification", back_populates="family")


class NavisNames(BaseModel):
    __tablename__ = "navis_names"

    id = Column(Integer, primary_key=True)
    navis_name = Column(String)
    previous_navis_name = Column(String)
    ship_family_id = Column(String)  # Keep for backward compatibility
    family_id = Column(Integer, ForeignKey("family_names.id"))

    # Relationships
    family = relationship("FamilyNames", back_populates="navis_names")
    haplotype_names = relationship("HaplotypeNames", back_populates="navis")
    joined_ships = relationship("JoinedShips", back_populates="navis")
    classifications = relationship("Classification", back_populates="navis")


class HaplotypeNames(BaseModel):
    __tablename__ = "haplotype_names"

    id = Column(Integer, primary_key=True)
    haplotype_name = Column(String)
    previous_haplotype_name = Column(String)
    family_id = Column(Integer, ForeignKey("family_names.id"))
    navis_id = Column(Integer, ForeignKey("navis_names.id"))

    # Relationships
    family = relationship("FamilyNames", back_populates="haplotype_names")
    navis = relationship("NavisNames", back_populates="haplotype_names")
    joined_ships = relationship("JoinedShips", back_populates="haplotype")
    classifications = relationship("Classification", back_populates="haplotype")


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
    species_group = Column(String(50))
    subgenus = Column(String(50))
    strain = Column(String(50))

    # Relationships
    joined_ships = relationship("JoinedShips", back_populates="taxonomy")
    classifications = relationship("Classification", back_populates="taxonomy")


class Gff(Base):
    __tablename__ = "gff"
    __table_args__ = (
        Index("idx_gff_accession_id", "accession_id"),
        Index("idx_gff_ship_id", "ship_id"),
    )

    id = Column(Integer, primary_key=True)
    contig = Column(String(100))
    accession_id = Column(Integer, ForeignKey("accessions.id"))
    ship_id = Column(Integer, ForeignKey("ships.id"))
    source = Column(String)
    type = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    phase = Column(String)
    strand = Column(String)
    score = Column(String)
    attributes = Column(String)

    # Relationships
    accession = relationship("Accessions", back_populates="gff_entries")


class Classification(BaseModel):
    __tablename__ = "classification"

    id = Column(Integer, primary_key=True)
    Family = Column(String)
    Species = Column(String)
    Strain = Column(String)
    Code = Column(String)
    Navis = Column(String)
    Size = Column(Integer)
    DR = Column(String)
    aTIR = Column(String)
    Target = Column(String)
    Spok = Column(String)
    ARS = Column(String)
    Other = Column(String)
    HGT = Column(String)
    taxID = Column(String)

    # Foreign keys
    accession_id = Column(Integer, ForeignKey("accessions.id"))
    taxonomy_id = Column(Integer, ForeignKey("taxonomy.id"))
    genome_id = Column(Integer, ForeignKey("genomes.id"))
    family_id = Column(Integer, ForeignKey("family_names.id"))
    navis_id = Column(Integer, ForeignKey("navis_names.id"))
    haplotype_id = Column(Integer, ForeignKey("haplotype_names.id"))

    # Relationships
    accession = relationship("Accessions", back_populates="classifications")
    taxonomy = relationship("Taxonomy", back_populates="classifications")
    genome = relationship("Genome", back_populates="classifications")
    family = relationship("FamilyNames", back_populates="classifications")
    navis = relationship("NavisNames", back_populates="classifications")
    haplotype = relationship("HaplotypeNames", back_populates="classifications")
