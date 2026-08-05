# coding: utf-8
from sqlalchemy import (
    Column,
    ForeignKey,
    Integer,
    String,
    Text,
    DateTime,
    VARCHAR,
    Boolean,
    UniqueConstraint,
    func,
)
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata


class DatabaseVersion(Base):
    """Content version changelog — one row per version bump event."""

    __tablename__ = "database_versions"

    id = Column(Integer, primary_key=True, autoincrement=True)
    semantic_version = Column(String, nullable=False)
    description = Column(Text, nullable=True)
    created_at = Column(DateTime, server_default=func.now())
    created_by = Column(String, default="system")


class Accessions(Base):
    __tablename__ = "accessions"
    id = Column(Integer, primary_key=True)
    accession_tag = Column(String, unique=True, nullable=False)
    version_tag = Column(Integer, nullable=False)

    # Relationships
    ships = relationship("Ships", back_populates="accession_obj")
    gff_entries = relationship("Gff", back_populates="accession_obj")
    joined_ship = relationship("JoinedShips", back_populates="accession", uselist=False)


class ShipAccessions(Base):
    __tablename__ = "ship_accessions"
    id = Column(Integer, primary_key=True)
    ship_accession_tag = Column(String, unique=True, nullable=False)
    ship_version_tag = Column(Integer, nullable=False, default=1)
    ship_id = Column(Integer, ForeignKey("ships.id"), unique=True, nullable=False)
    ship_accession_display = Column(String, nullable=True)

    # Relationships
    ship = relationship("Ships", back_populates="ship_accession")


class Ships(Base):
    __tablename__ = "ships"
    id = Column(Integer, primary_key=True)
    sequence = Column(String)
    md5 = Column(String)
    sequence_length = Column(Integer)
    header = Column(String)
    rev_comp_md5 = Column(String)
    accession_id = Column(Integer, ForeignKey("accessions.id"))
    type_ship = Column(String)

    # Relationships
    accession_obj = relationship("Accessions", back_populates="ships")
    ship_accession = relationship(
        "ShipAccessions", back_populates="ship", uselist=False
    )


class Captains(Base):
    __tablename__ = "captains"
    id = Column(Integer, primary_key=True)
    captainID = Column(String, unique=True)
    sequence = Column(String)
    md5 = Column(String)
    sequence_length = Column(Integer)
    ship_id = Column(Integer, ForeignKey("ships.id"))
    reviewed = Column(String)
    evidence = Column(String)

    # Relationships
    ship = relationship("Ships")
    features = relationship("StarshipFeatures", back_populates="captain")


class Genome(Base):
    __tablename__ = "genomes"

    id = Column(Integer, primary_key=True)
    ome = Column(String(50))
    version = Column(String(50))
    genomeSource = Column(String(50))
    citation = Column(String(50))
    biosample = Column(String(50))
    acquisition_date = Column(Integer)
    assembly_accession = Column(String(50))
    taxonomy_id = Column(Integer, ForeignKey("taxonomy.id"))

    # Relationships
    taxonomy = relationship("Taxonomy", back_populates="genomes")


class StarshipFeatures(Base):
    __tablename__ = "starship_features"
    id = Column(Integer, primary_key=True)
    accession_id = Column(Integer, ForeignKey("accessions.id"))
    ship_id = Column(Integer, ForeignKey("ships.id"))
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
    dr = Column(String)
    tir = Column(String)
    target = Column(String)
    spok = Column(String)
    ars = Column(String)
    other = Column(String)
    hgt = Column(String)
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
    class_ = Column(String, name="class")  # Using class_ as class is a reserved keyword
    subclass = Column(VARCHAR)
    order = Column(VARCHAR)
    suborder = Column(VARCHAR)
    family = Column(VARCHAR)
    genus = Column(VARCHAR)
    species = Column(VARCHAR)
    section = Column(VARCHAR)
    species_group = Column(VARCHAR)
    subgenus = Column(VARCHAR)
    strain = Column(VARCHAR)

    # Relationships
    genomes = relationship("Genome", back_populates="taxonomy")


class Navis(Base):
    __tablename__ = "navis_names"
    id = Column(Integer, primary_key=True)
    navis_name = Column(String)
    previous_navis_name = Column(String)
    ship_family_id = Column(Integer, ForeignKey("family_names.id"))
    activity = Column(
        Integer, default=1, nullable=True
    )  # 0 = inactive, don't show in modals

    # Relationships
    family = relationship("FamilyNames")


class Haplotype(Base):
    __tablename__ = "haplotype_names"
    id = Column(Integer, primary_key=True)
    haplotype_name = Column(String)
    previous_haplotype_name = Column(String)
    navis_id = Column(Integer, ForeignKey("navis_names.id"))
    ship_family_id = Column(Integer, ForeignKey("family_names.id"))
    activity = Column(
        Integer, default=1, nullable=True
    )  # 0 = inactive, don't show in modals

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
    starshipID = Column(
        String, nullable=False
    )  # Every ship needs an ID (duplicates allowed)
    evidence = Column(String)
    source = Column(String)
    curated_status = Column(String)

    # Direct link to accession (when sequence data is available)
    accession_id = Column(Integer, ForeignKey("accessions.id"), nullable=True)
    # Direct link to ship-level (SSB) accession; mirrors accession_id (SSA, group-level).
    # Previously only reachable via ship_id -> ships.id -> ship_accessions.ship_id.
    ship_accession_id = Column(Integer, ForeignKey("ship_accessions.id"), nullable=True)

    # Links to classification and annotation data
    ship_id = Column(Integer, ForeignKey("ships.id"))
    ship_family_id = Column(Integer, ForeignKey("family_names.id"))
    tax_id = Column(Integer, ForeignKey("taxonomy.id"))
    genome_id = Column(Integer, ForeignKey("genomes.id"))
    captain_id = Column(Integer, ForeignKey("captains.id"))
    ship_navis_id = Column(Integer, ForeignKey("navis_names.id"))
    ship_haplotype_id = Column(Integer, ForeignKey("haplotype_names.id"))
    created_at = Column(DateTime)
    updated_at = Column(DateTime)
    is_deleted = Column(Boolean, nullable=False, default=False, server_default="0")

    # Relationships
    accession = relationship("Accessions", back_populates="joined_ship")
    ship_accession = relationship("ShipAccessions")
    ship = relationship("Ships")
    family = relationship("FamilyNames")
    taxonomy = relationship("Taxonomy")
    genome = relationship("Genome")
    captain = relationship("Captains")
    navis = relationship("Navis")
    haplotype = relationship("Haplotype")
    quality_tags = relationship(
        "ShipQualityTags", back_populates="joined_ship", cascade="all, delete-orphan"
    )


class ShipQualityTags(Base):
    __tablename__ = "ship_quality_tags"
    id = Column(Integer, primary_key=True)
    joined_ship_id = Column(Integer, ForeignKey("joined_ships.id"), nullable=False)
    tag_type = Column(String(50), nullable=False)
    tag_value = Column(
        String(100), nullable=True
    )  # Optional: for tags that need values
    created_at = Column(DateTime, nullable=False, server_default="CURRENT_TIMESTAMP")
    created_by = Column(String(50), nullable=True, default="auto")

    # Relationships
    joined_ship = relationship("JoinedShips", back_populates="quality_tags")

    # Unique constraint to prevent duplicate tags per ship
    __table_args__ = ({"sqlite_autoincrement": True},)


class Submission(Base):
    """Submission records in the submissions database (queue for processing)."""

    __tablename__ = "submissions"
    id = Column(Integer, primary_key=True, autoincrement=True)
    seq_contents = Column(Text, nullable=False)
    seq_filename = Column(String(255), nullable=False)
    seq_date = Column(String(50), nullable=True)
    anno_contents = Column(Text, nullable=True)
    anno_filename = Column(String(255), nullable=True)
    anno_date = Column(String(50), nullable=True)
    uploader = Column(String(255), nullable=True)
    evidence = Column(String(100), nullable=True)
    genus = Column(String(255), nullable=True)
    species = Column(String(255), nullable=True)
    strain = Column(String(255), nullable=True)
    hostchr = Column(String(255), nullable=True)
    assembly_accession = Column(String(50), nullable=True)
    shipstart = Column(Integer, nullable=True)
    shipend = Column(Integer, nullable=True)
    shipstrand = Column(String(10), nullable=True)
    comment = Column(Text, nullable=True)
    ship_accession_tag = Column(String(50), nullable=True)
    accession_tag = Column(String(50), nullable=True)
    needs_review = Column(Boolean, default=False, nullable=True)
    submission_group_id = Column(String(36), nullable=True)
    processing_status = Column(String(20), default="pending", nullable=True)
    # Classification from BLAST prefill (optional)
    classification_source = Column(String(50), nullable=True)
    classification_family = Column(String(100), nullable=True)
    classification_navis = Column(String(100), nullable=True)
    classification_haplotype = Column(String(100), nullable=True)
    closest_match = Column(String(50), nullable=True)
    classification_confidence = Column(String(20), nullable=True)


class Annotation(Base):
    """Curated gene/protein annotation, deduplicated by sequence.

    Ported from MAS4starships' starship.Annotation. No FK to a ship directly --
    reachable only via GeneFeature (one Annotation can be shared by genes on
    multiple ships, deduped by identical protein sequence).
    """

    __tablename__ = "annotations"
    id = Column(Integer, primary_key=True)
    sequence = Column(Text, nullable=False, unique=True)
    annotation = Column(String(255), nullable=True, default="No Annotation")
    public_notes = Column(Text, nullable=True, default="")
    private_notes = Column(Text, nullable=True, default="")
    # 0=GREEN, 1=YELLOW, 2=RED, 3=REVIEW_NAME, 4=N_A, 5=ORANGE, 7=UNANNOTATED
    flag = Column(Integer, nullable=False, default=7)
    assigned_to = Column(String(255), nullable=True)
    created_at = Column(DateTime, nullable=False, server_default=func.now())
    updated_at = Column(DateTime, nullable=True)

    gene_features = relationship("GeneFeature", back_populates="annotation")
    history = relationship(
        "AnnotationHistory", back_populates="annotation", cascade="all, delete-orphan"
    )


class GeneFeature(Base):
    """One gene/CDS/repeat/tRNA call within a ship.

    Not the same thing as StarshipFeatures (the starship element's own
    begin/end within its host genome) -- this is per-gene, within a ship.
    Ported from MAS4starships' starship.Feature.
    """

    __tablename__ = "gene_features"
    id = Column(Integer, primary_key=True)
    joined_ship_id = Column(Integer, ForeignKey("joined_ships.id"), nullable=False)
    annotation_id = Column(Integer, ForeignKey("annotations.id"), nullable=True)
    start = Column(Integer, nullable=False, default=0)
    stop = Column(Integer, nullable=False, default=0)
    type = Column(String(50), nullable=True)  # gene, CDS, Repeat Region, tRNA
    strand = Column(String(1), nullable=True)

    joined_ship = relationship("JoinedShips")
    annotation = relationship("Annotation", back_populates="gene_features")


class AnnotationHistory(Base):
    """Flat audit trail for Annotation edits.

    Not a django-simple-history clone -- one row per mutation with before/after
    values, written by the curation API on every save.
    """

    __tablename__ = "annotation_history"
    id = Column(Integer, primary_key=True)
    annotation_id = Column(Integer, ForeignKey("annotations.id"), nullable=False)
    changed_by = Column(String(255), nullable=True)
    changed_at = Column(DateTime, nullable=False, server_default=func.now())
    old_flag = Column(Integer, nullable=True)
    new_flag = Column(Integer, nullable=True)
    old_annotation = Column(Text, nullable=True)
    new_annotation = Column(Text, nullable=True)
    old_public_notes = Column(Text, nullable=True)
    new_public_notes = Column(Text, nullable=True)
    old_private_notes = Column(Text, nullable=True)
    new_private_notes = Column(Text, nullable=True)

    annotation = relationship("Annotation", back_populates="history")


# The 4 homology-search result tables share one shape: one row per
# (annotation, database) pair. Result generation is a deferred slice -- these
# tables exist so the viewer has somewhere to read from. `result` holds raw
# hit data (not a file path -- the backend owns the DB directly, no separate
# file store to manage). status: 0=complete, 1=running, 2=error.


class BlastpResult(Base):
    __tablename__ = "blastp_results"
    id = Column(Integer, primary_key=True)
    annotation_id = Column(Integer, ForeignKey("annotations.id"), nullable=False)
    database = Column(String(50), nullable=False)  # swissprot, protein, nr
    result = Column(Text, nullable=True)
    run_date = Column(DateTime, nullable=True)
    status = Column(Integer, nullable=False, default=0)
    __table_args__ = (UniqueConstraint("annotation_id", "database"),)


class RpsblastResult(Base):
    __tablename__ = "rpsblast_results"
    id = Column(Integer, primary_key=True)
    annotation_id = Column(Integer, ForeignKey("annotations.id"), nullable=False)
    database = Column(String(50), nullable=False)  # cdd
    result = Column(Text, nullable=True)
    run_date = Column(DateTime, nullable=True)
    status = Column(Integer, nullable=False, default=0)
    __table_args__ = (UniqueConstraint("annotation_id", "database"),)


class HhsearchResult(Base):
    __tablename__ = "hhsearch_results"
    id = Column(Integer, primary_key=True)
    annotation_id = Column(Integer, ForeignKey("annotations.id"), nullable=False)
    database = Column(String(50), nullable=False)  # pdb
    result = Column(Text, nullable=True)
    run_date = Column(DateTime, nullable=True)
    status = Column(Integer, nullable=False, default=0)
    __table_args__ = (UniqueConstraint("annotation_id", "database"),)


class InterproResult(Base):
    __tablename__ = "interpro_results"
    id = Column(Integer, primary_key=True)
    annotation_id = Column(Integer, ForeignKey("annotations.id"), nullable=False)
    database = Column(String(50), nullable=False)  # interpro
    result = Column(Text, nullable=True)
    run_date = Column(DateTime, nullable=True)
    status = Column(Integer, nullable=False, default=0)
    __table_args__ = (UniqueConstraint("annotation_id", "database"),)
