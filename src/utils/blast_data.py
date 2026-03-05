"""
Consolidated data models for the BLAST/classification pipeline.

This file contains all dataclasses, state management, and adapters needed for
the pipeline. Everything is consolidated here for simplicity and maintainability.

Key Components
==============
1. **Consolidated Data**:
   - SequenceAnalysis: Single comprehensive data structure
   - ClassificationResult: Clean classification data with enums
   - BlastResult: Simplified BLAST data
   - WorkflowConfig: Unified configuration

2. **Legacy Models** (Backward compatibility):
   - ClassificationData, BlastData, WorkflowState: Original dataclasses
   - FetchShipParams, FetchCaptainParams: Original config classes

3. **State Management**:
   - ConsolidatedPipelineState: New simplified state manager
   - PipelineState: Legacy state manager
   - Migration adapters for gradual transition

4. **Serialization**:
   All dataclasses expose to_dict()/from_dict() for JSON compatibility with Dash stores.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, List, Dict, Any
from src.config.logging import get_logger
from src.config.cache import cache, DEFAULT_CACHE_TIMEOUT

logger = get_logger(__name__)

# Cache key prefix for pipeline state (shared across workers via filesystem cache)
_BLAST_PIPELINE_CACHE_KEY_PREFIX = "blast_pipeline:"

# =============================================================================
# CONSOLIDATED DATA MODELS
# =============================================================================


class SequenceType(Enum):
    """Sequence type enumeration"""

    NUCLEOTIDE = "nucl"
    PROTEIN = "prot"


class WorkflowStatus(Enum):
    """Workflow status enumeration"""

    INITIALIZED = "initialized"
    RUNNING = "running"
    COMPLETE = "complete"
    FAILED = "failed"


class MatchStage(Enum):
    """Classification match stage enumeration"""

    EXACT = "exact"
    CONTAINED = "contained"
    SIMILAR = "similar"
    FAMILY = "family"
    NAVIS = "navis"
    HAPLOTYPE = "haplotype"


class ConfidenceLevel(Enum):
    """Classification confidence level"""

    HIGH = "High"
    MEDIUM = "Medium"
    LOW = "Low"


@dataclass
class ClassificationResult:
    """
    Consolidated classification result.

    This replaces ClassificationData and contains all classification information
    in a single, clean structure.
    """

    # Match information
    stage: Optional[MatchStage] = None
    closest_match: Optional[str] = None
    confidence: Optional[ConfidenceLevel] = None
    match_details: Optional[str] = None

    # Taxonomic classification
    family: Optional[str] = None
    navis: Optional[str] = None
    haplotype: Optional[str] = None

    # Technical details
    sequence_type: SequenceType = SequenceType.NUCLEOTIDE

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return {
            "stage": self.stage.value if self.stage else None,
            "closest_match": self.closest_match,
            "confidence": self.confidence.value if self.confidence else None,
            "match_details": self.match_details,
            "family": self.family,
            "navis": self.navis,
            "haplotype": self.haplotype,
            "sequence_type": self.sequence_type.value,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ClassificationResult":
        """Create from dictionary"""
        if not data:
            return None

        return cls(
            stage=MatchStage(data["stage"]) if data.get("stage") else None,
            closest_match=data.get("closest_match"),
            confidence=ConfidenceLevel(data["confidence"])
            if data.get("confidence")
            else None,
            match_details=data.get("match_details"),
            family=data.get("family"),
            navis=data.get("navis"),
            haplotype=data.get("haplotype"),
            sequence_type=SequenceType(data.get("sequence_type", "nucl")),
        )

    def is_empty(self) -> bool:
        """Check if classification is empty"""
        return all(
            v is None
            for v in [
                self.stage,
                self.closest_match,
                self.family,
                self.navis,
                self.haplotype,
            ]
        )


@dataclass
class BlastResult:
    """
    Consolidated BLAST result information.

    This replaces BlastData and contains only essential BLAST information.
    """

    # Input data
    sequence: Optional[str] = None
    sequence_type: SequenceType = SequenceType.NUCLEOTIDE
    fasta_file: Optional[str] = None

    # BLAST outputs
    blast_content: Optional[str] = None
    blast_file: Optional[str] = None
    blast_hits: Optional[List[Dict[str, Any]]] = None  # Parsed BLAST hits

    # Processing status
    processed: bool = False
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return {
            "sequence": self.sequence,
            "sequence_type": self.sequence_type.value,
            "fasta_file": self.fasta_file,
            "blast_content": self.blast_content,
            "blast_file": self.blast_file,
            "blast_hits": self.blast_hits,
            "processed": self.processed,
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "BlastResult":
        """Create from dictionary"""
        if not data:
            return None

        return cls(
            sequence=data.get("sequence"),
            sequence_type=SequenceType(data.get("sequence_type", "nucl")),
            fasta_file=data.get("fasta_file"),
            blast_content=data.get("blast_content"),
            blast_file=data.get("blast_file"),
            blast_hits=data.get("blast_hits"),
            processed=data.get("processed", False),
            error=data.get("error"),
        )


@dataclass
class WorkflowConfig:
    """
    Simplified workflow configuration.

    This replaces FetchShipParams and FetchCaptainParams with a single config object.
    """

    # Ship database settings
    ship_curated_only: bool = False
    ship_include_sequence: bool = True
    ship_dereplicate: bool = True

    # Captain database settings
    captain_curated_only: bool = True
    captain_include_sequence: bool = True

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            "ship_curated_only": self.ship_curated_only,
            "ship_include_sequence": self.ship_include_sequence,
            "ship_dereplicate": self.ship_dereplicate,
            "captain_curated_only": self.captain_curated_only,
            "captain_include_sequence": self.captain_include_sequence,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "WorkflowConfig":
        """Create from dictionary"""
        if not data:
            return cls()  # Return default config

        return cls(
            ship_curated_only=data.get("ship_curated_only", False),
            ship_include_sequence=data.get("ship_include_sequence", True),
            ship_dereplicate=data.get("ship_dereplicate", True),
            captain_curated_only=data.get("captain_curated_only", True),
            captain_include_sequence=data.get("captain_include_sequence", True),
        )


@dataclass
class SequenceAnalysis:
    """
    Consolidated sequence analysis state.

    This replaces SequenceState, WorkflowState, BlastData, and ClassificationData
    with a single, comprehensive data structure.
    """

    # Identity
    sequence_id: str

    # Input data
    sequence: Optional[str] = None
    sequence_header: Optional[str] = None
    sequence_type: SequenceType = SequenceType.NUCLEOTIDE

    # BLAST analysis
    blast_result: Optional[BlastResult] = None

    # Classification analysis
    classification: Optional[ClassificationResult] = None

    # Workflow state
    status: WorkflowStatus = WorkflowStatus.INITIALIZED
    config: WorkflowConfig = field(default_factory=WorkflowConfig)

    # Progress tracking
    current_stage: Optional[str] = None
    stage_progress: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    # Results and errors
    error: Optional[str] = None
    completed_at: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return {
            "sequence_id": self.sequence_id,
            "sequence": self.sequence,
            "sequence_header": self.sequence_header,
            "sequence_type": self.sequence_type.value,
            "blast_result": self.blast_result.to_dict() if self.blast_result else None,
            "classification": self.classification.to_dict()
            if self.classification
            else None,
            "status": self.status.value,
            "config": self.config.to_dict(),
            "current_stage": self.current_stage,
            "stage_progress": self.stage_progress,
            "error": self.error,
            "completed_at": self.completed_at,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SequenceAnalysis":
        """Create from dictionary"""
        if not data:
            return None

        return cls(
            sequence_id=data["sequence_id"],
            sequence=data.get("sequence"),
            sequence_header=data.get("sequence_header"),
            sequence_type=SequenceType(data.get("sequence_type", "nucl")),
            blast_result=BlastResult.from_dict(data.get("blast_result")),
            classification=ClassificationResult.from_dict(data.get("classification")),
            status=WorkflowStatus(data.get("status", "initialized")),
            config=WorkflowConfig.from_dict(data.get("config", {})),
            current_stage=data.get("current_stage"),
            stage_progress=data.get("stage_progress", {}),
            error=data.get("error"),
            completed_at=data.get("completed_at"),
        )

    # Convenience methods
    def is_complete(self) -> bool:
        """Check if analysis is complete"""
        return self.status == WorkflowStatus.COMPLETE

    def has_error(self) -> bool:
        """Check if analysis has an error"""
        return self.error is not None or self.status == WorkflowStatus.FAILED

    def has_classification(self) -> bool:
        """Check if analysis has classification results"""
        return self.classification is not None and not self.classification.is_empty()

    def set_error(self, error: str):
        """Set error and update status"""
        self.error = error
        self.status = WorkflowStatus.FAILED

    def set_complete(self):
        """Mark analysis as complete"""
        self.status = WorkflowStatus.COMPLETE
        from datetime import datetime

        self.completed_at = datetime.now().isoformat()


# Legacy compatibility functions for gradual migration
def convert_from_legacy_blast_data(blast_data_dict: Dict[str, Any]) -> BlastResult:
    """Convert legacy BlastData dict to new BlastResult"""
    if not blast_data_dict:
        return None

    return BlastResult(
        sequence=blast_data_dict.get("sequence"),
        sequence_type=SequenceType(blast_data_dict.get("seq_type", "nucl")),
        fasta_file=blast_data_dict.get("fasta_file"),
        blast_content=blast_data_dict.get("blast_content"),
        blast_file=blast_data_dict.get("blast_file"),
        blast_hits=blast_data_dict.get("blast_df"),
        processed=blast_data_dict.get("processed", False),
        error=blast_data_dict.get("error"),
    )


def convert_from_legacy_classification_data(
    classification_dict: Dict[str, Any],
) -> ClassificationResult:
    """Convert legacy ClassificationData dict to new ClassificationResult"""
    if not classification_dict:
        return None

    # Map legacy field names
    stage_mapping = {
        "exact": MatchStage.EXACT,
        "contained": MatchStage.CONTAINED,
        "similar": MatchStage.SIMILAR,
        "family": MatchStage.FAMILY,
        "navis": MatchStage.NAVIS,
        "haplotype": MatchStage.HAPLOTYPE,
    }

    confidence_mapping = {
        "High": ConfidenceLevel.HIGH,
        "Medium": ConfidenceLevel.MEDIUM,
        "Low": ConfidenceLevel.LOW,
    }

    legacy_source = classification_dict.get("source")
    legacy_confidence = classification_dict.get("confidence")

    return ClassificationResult(
        stage=stage_mapping.get(legacy_source) if legacy_source else None,
        closest_match=classification_dict.get("closest_match"),
        confidence=confidence_mapping.get(legacy_confidence)
        if legacy_confidence
        else None,
        match_details=classification_dict.get("match_details"),
        family=classification_dict.get("family"),
        navis=classification_dict.get("navis"),
        haplotype=classification_dict.get("haplotype"),
        sequence_type=SequenceType(classification_dict.get("seq_type", "nucl")),
    )


# =============================================================================
# LEGACY MODELS (Backward compatibility)
# =============================================================================


@dataclass
class FetchShipParams:
    curated: bool = False
    with_sequence: bool = True
    dereplicate: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "curated": self.curated,
            "with_sequence": self.with_sequence,
            "dereplicate": self.dereplicate,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "FetchShipParams":
        return cls(**data)


@dataclass
class FetchCaptainParams:
    curated: bool = False
    with_sequence: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "curated": self.curated,
            "with_sequence": self.with_sequence,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "FetchCaptainParams":
        return cls(**data)


@dataclass
class ClassificationData:
    seq_type: Optional[str] = "nucl"
    fasta_file: Optional[str] = None
    source: Optional[str] = None
    family: Optional[str] = None
    navis: Optional[str] = None
    haplotype: Optional[str] = None
    closest_match: Optional[str] = None
    match_details: Optional[str] = None
    confidence: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "seq_type": self.seq_type,
            "fasta_file": self.fasta_file,
            "source": self.source,
            "family": self.family,
            "navis": self.navis,
            "haplotype": self.haplotype,
            "closest_match": self.closest_match,
            "match_details": self.match_details,
            "confidence": self.confidence,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ClassificationData":
        if data is None:
            return None
        return cls(**data)

    def is_empty(self) -> bool:
        """Check if the classification data is empty (no meaningful data)"""
        return all(
            getattr(self, field) is None
            for field in [
                "seq_type",
                "fasta_file",
                "source",
                "family",
                "navis",
                "haplotype",
                "closest_match",
                "match_details",
                "confidence",
            ]
        )


@dataclass
class BlastData:
    seq_type: Optional[str] = "nucl"  # Default to nucleotide sequences
    processed_sequences: Optional[List[int]] = None
    sequence_results: Optional[Dict[str, Any]] = None
    total_sequences: int = 0
    blast_df: Optional[List[Dict[str, Any]]] = None
    blast_content: Optional[str] = None
    blast_file: Optional[str] = None
    sequence: Optional[str] = None
    fasta_file: Optional[str] = None
    error: Optional[str] = None
    processed: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Convert the data class to a dictionary for backward compatibility"""
        result = {
            "seq_type": self.seq_type,
            "processed_sequences": self.processed_sequences,
            "sequence_results": self.sequence_results,
            "total_sequences": self.total_sequences,
            "blast_df": self.blast_df,
            "blast_content": self.blast_content,
            "blast_file": self.blast_file,
            "sequence": self.sequence,
            "fasta_file": self.fasta_file,
            "error": self.error,
            "processed": self.processed,
        }
        return result

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "BlastData":
        """Create a BlastData instance from a dictionary"""
        if data is None:
            return None

        # Create the class instance with extracted data
        return cls(
            seq_type=data.get("seq_type", "nucl"),
            processed_sequences=data.get("processed_sequences", []),
            sequence_results=data.get("sequence_results", []),
            total_sequences=data.get("total_sequences", 0),
            blast_df=data.get("blast_df"),
            blast_content=data.get("blast_content"),
            blast_file=data.get("blast_file"),
            sequence=data.get("sequence"),
            fasta_file=data.get("fasta_file"),
            error=data.get("error"),
            processed=data.get("processed", False),
        )


@dataclass
class WorkflowState:
    complete: bool = False
    error: Optional[str] = None
    found_match: bool = False
    match_stage: Optional[str] = None
    match_result: Optional[str] = None
    classification_data: Optional[Dict[str, Any]] = None
    stages: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    task_id: str = ""
    status: str = "initialized"
    workflow_started: bool = True
    current_stage: Optional[str] = None
    current_stage_idx: int = 0
    fetch_ship_params: FetchShipParams = field(default_factory=FetchShipParams)
    fetch_captain_params: FetchCaptainParams = field(default_factory=FetchCaptainParams)

    def to_dict(self) -> Dict[str, Any]:
        """Convert WorkflowState to dictionary for JSON serialization."""
        return {
            "complete": self.complete,
            "error": self.error,
            "found_match": self.found_match,
            "match_stage": self.match_stage,
            "match_result": self.match_result,
            "classification_data": self.classification_data,
            "stages": self.stages,
            "task_id": self.task_id,
            "status": self.status,
            "workflow_started": self.workflow_started,
            "current_stage": self.current_stage,
            "current_stage_idx": self.current_stage_idx,
            "fetch_ship_params": self.fetch_ship_params.to_dict(),
            "fetch_captain_params": self.fetch_captain_params.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "WorkflowState":
        """Create a WorkflowState instance from a dictionary"""
        if data is None:
            return None

        # Handle nested dictionaries
        fetch_ship_params = FetchShipParams.from_dict(data.get("fetch_ship_params", {}))
        fetch_captain_params = FetchCaptainParams.from_dict(
            data.get("fetch_captain_params", {})
        )

        return cls(
            complete=data.get("complete", False),
            error=data.get("error"),
            found_match=data.get("found_match", False),
            match_stage=data.get("match_stage"),
            match_result=data.get("match_result"),
            classification_data=data.get("classification_data"),
            stages=data.get("stages", {}),
            task_id=data.get("task_id", ""),
            status=data.get("status", "initialized"),
            workflow_started=data.get("workflow_started", True),
            current_stage=data.get("current_stage"),
            current_stage_idx=data.get("current_stage_idx", 0),
            fetch_ship_params=fetch_ship_params,
            fetch_captain_params=fetch_captain_params,
        )

    def set_classification(self, classification_data: "ClassificationData") -> None:
        """Set classification results on the workflow state"""
        self.found_match = True
        self.match_stage = classification_data.source
        self.match_result = classification_data.closest_match
        self.classification_data = (
            classification_data.to_dict() if classification_data else None
        )


# =============================================================================
# STATE MANAGEMENT
# =============================================================================


@dataclass
class SequenceState:
    """State for a single sequence in the pipeline"""

    sequence_id: str
    blast_data: Optional[BlastData] = None
    workflow_state: Optional[WorkflowState] = None
    classification_data: Optional[ClassificationData] = None
    processed: bool = False
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return {
            "sequence_id": self.sequence_id,
            "blast_data": self.blast_data.to_dict() if self.blast_data else None,
            "workflow_state": self.workflow_state.to_dict()
            if self.workflow_state
            else None,
            "classification_data": self.classification_data.to_dict()
            if self.classification_data
            else None,
            "processed": self.processed,
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SequenceState":
        """Create from dictionary"""
        return cls(
            sequence_id=data["sequence_id"],
            blast_data=BlastData.from_dict(data["blast_data"])
            if data.get("blast_data")
            else None,
            workflow_state=WorkflowState.from_dict(data["workflow_state"])
            if data.get("workflow_state")
            else None,
            classification_data=ClassificationData.from_dict(
                data["classification_data"]
            )
            if data.get("classification_data")
            else None,
            processed=data.get("processed", False),
            error=data.get("error"),
        )


class PipelineState:
    """
    Centralized state manager for the classification pipeline.

    This eliminates data duplication by providing a single source of truth
    for all sequence data, workflow states, and classification results.

    When running with multiple workers (e.g. production), state is
    backed by the filesystem cache so all workers share the same pipeline
    state for a given submission (keyed by session blast_submission_id).
    """

    def __init__(self):
        self._sequences: Dict[str, SequenceState] = {}
        self._active_sequence_id: Optional[str] = None
        # Optional cache backing: when set, state is persisted after each mutation
        self._cache_key: Optional[str] = None
        self._cache_timeout: int = DEFAULT_CACHE_TIMEOUT

    def _maybe_persist(self) -> None:
        """Persist state to filesystem cache when cache backing is configured."""
        if not self._cache_key:
            return
        try:
            cache.set(
                self._cache_key,
                self.to_dict(),
                timeout=self._cache_timeout,
            )
        except Exception as e:
            logger.warning("Failed to persist pipeline state to cache: %s", e)

    def add_sequence(
        self, sequence_id: str, clear_existing: bool = False
    ) -> SequenceState:
        """Add a new sequence to track"""
        if sequence_id in self._sequences:
            if clear_existing:
                logger.info(
                    f"Clearing existing sequence {sequence_id} and creating new one"
                )
                # Clear the existing sequence state
                del self._sequences[sequence_id]
            else:
                logger.warning(
                    f"Sequence {sequence_id} already exists, returning existing"
                )
                return self._sequences[sequence_id]

        sequence_state = SequenceState(sequence_id=sequence_id)
        self._sequences[sequence_id] = sequence_state

        # Set as active sequence
        old_active_id = self._active_sequence_id
        self._active_sequence_id = sequence_id
        if old_active_id != sequence_id:
            logger.debug(
                f"Changed active sequence ID from {old_active_id} to {sequence_id}"
            )

        logger.debug(f"Added sequence {sequence_id} to pipeline state")
        self._maybe_persist()
        return sequence_state

    def add_sequence_without_activation(
        self, sequence_id: str, clear_existing: bool = False
    ) -> SequenceState:
        """Add a new sequence to track without changing the active sequence ID"""
        if sequence_id in self._sequences:
            if clear_existing:
                logger.info(
                    f"Clearing existing sequence {sequence_id} and creating new one"
                )
                # Clear the existing sequence state
                del self._sequences[sequence_id]
            else:
                logger.warning(
                    f"Sequence {sequence_id} already exists, returning existing"
                )
                return self._sequences[sequence_id]

        sequence_state = SequenceState(sequence_id=sequence_id)
        self._sequences[sequence_id] = sequence_state

        # Don't change the active sequence ID
        logger.debug(
            f"Added sequence {sequence_id} to pipeline state (without activation)"
        )
        self._maybe_persist()
        return sequence_state

    def get_sequence(self, sequence_id: str) -> Optional[SequenceState]:
        """Get sequence state by ID"""
        if sequence_id:
            return self._sequences.get(sequence_id)
        return None

    def resolve_sequence_id(self, sequence_id: str) -> Optional[str]:
        """
        Resolve sequence ID by trying different formats.
        This helps when there are mismatches between tab IDs and base IDs.
        """
        if not sequence_id:
            return None

        # Direct match
        if sequence_id in self._sequences:
            return sequence_id

        # Try without tab suffix
        if "_tab_" in sequence_id:
            base_sequence_id = sequence_id.split("_tab_")[0]
            if base_sequence_id in self._sequences:
                return base_sequence_id

        # Try with tab suffix
        for stored_id in self._sequences.keys():
            if stored_id.startswith(sequence_id) or sequence_id.startswith(stored_id):
                return stored_id

        return None

    def get_active_sequence(self) -> Optional[SequenceState]:
        """Get the currently active sequence"""
        if self._active_sequence_id:
            return self._sequences.get(self._active_sequence_id)
        return None

    def set_active_sequence(self, sequence_id: str):
        """Set the active sequence"""
        if sequence_id in self._sequences:
            old_active_id = self._active_sequence_id
            self._active_sequence_id = sequence_id
            if old_active_id != sequence_id:
                logger.debug(
                    f"Explicitly set active sequence ID from {old_active_id} to {sequence_id}"
                )
        else:
            logger.error(f"Cannot set active sequence to {sequence_id} - not found")
        self._maybe_persist()

    def update_blast_data(self, sequence_id: str, blast_data: BlastData):
        """Update BLAST data for a sequence"""
        if sequence_id not in self._sequences:
            self.add_sequence_without_activation(sequence_id)

        self._sequences[sequence_id].blast_data = blast_data
        logger.debug(f"Updated BLAST data for sequence {sequence_id}")
        self._maybe_persist()

    def update_workflow_state(self, sequence_id: str, workflow_state: WorkflowState):
        """Update workflow state for a sequence"""
        if sequence_id not in self._sequences:
            self.add_sequence_without_activation(sequence_id)

        self._sequences[sequence_id].workflow_state = workflow_state
        logger.debug(f"Updated workflow state for sequence {sequence_id}")
        self._maybe_persist()

    def update_classification_data(
        self, sequence_id: str, classification_data: ClassificationData
    ):
        """Update classification data for a sequence - SINGLE SOURCE OF TRUTH"""
        if sequence_id not in self._sequences:
            self.add_sequence_without_activation(sequence_id)

        sequence_state = self._sequences[sequence_id]
        sequence_state.classification_data = classification_data

        # Automatically sync to workflow state if it exists
        if sequence_state.workflow_state:
            sequence_state.workflow_state.classification_data = (
                classification_data.to_dict()
            )

        logger.debug(f"Updated classification data for sequence {sequence_id}")
        logger.debug(
            f"Classification: family={classification_data.family}, navis={classification_data.navis}, haplotype={classification_data.haplotype}"
        )
        self._maybe_persist()

    def get_classification_data(self, sequence_id: str) -> Optional[ClassificationData]:
        """Get classification data for a sequence - SINGLE SOURCE OF TRUTH"""
        sequence_state = self.get_sequence(sequence_id)
        if sequence_state:
            return sequence_state.classification_data
        return None

    def mark_sequence_processed(self, sequence_id: str):
        """Mark a sequence as processed"""
        if sequence_id in self._sequences:
            self._sequences[sequence_id].processed = True
            logger.debug(f"Marked sequence {sequence_id} as processed")
            self._maybe_persist()

    def set_sequence_error(self, sequence_id: str, error: str):
        """Set an error for a sequence"""
        if sequence_id in self._sequences:
            self._sequences[sequence_id].error = error
            logger.error(f"Set error for sequence {sequence_id}: {error}")
            self._maybe_persist()

    def clear_all_sequences(self):
        """Clear all sequences from the pipeline state"""
        logger.info("Clearing all sequences from pipeline state")
        self._sequences.clear()
        self._active_sequence_id = None
        self._maybe_persist()

    def start_new_submission(self, sequence_id: str) -> SequenceState:
        """Start a new submission by clearing old state and adding the sequence"""
        logger.info(f"Starting new submission with sequence ID: {sequence_id}")
        # Clear all previous sequences to ensure clean state
        self.clear_all_sequences()
        # Add the new sequence
        sequence_state = self.add_sequence(sequence_id)
        logger.info(
            f"New submission initialized - active sequence: {self._active_sequence_id}"
        )
        self._maybe_persist()
        return sequence_state

    def get_processed_sequences(self) -> List[str]:
        """Get list of processed sequence IDs"""
        return [seq_id for seq_id, state in self._sequences.items() if state.processed]

    def to_dict(self) -> Dict[str, Any]:
        """Convert entire state to dictionary for serialization"""
        return {
            "sequences": {
                seq_id: state.to_dict() for seq_id, state in self._sequences.items()
            },
            "active_sequence_id": self._active_sequence_id,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "PipelineState":
        """Create from dictionary"""
        state = cls()
        if "sequences" in data:
            for seq_id, seq_data in data["sequences"].items():
                state._sequences[seq_id] = SequenceState.from_dict(seq_data)
        state._active_sequence_id = data.get("active_sequence_id")
        return state

    # Legacy compatibility methods for gradual migration
    def to_blast_data_dict(self, sequence_id: str = None) -> Optional[Dict[str, Any]]:
        """Convert to legacy BlastData format for backward compatibility"""
        if sequence_id is None:
            sequence_id = self._active_sequence_id

        sequence_state = self.get_sequence(sequence_id)
        if sequence_state and sequence_state.blast_data:
            blast_dict = sequence_state.blast_data.to_dict()

            # Inject classification data into sequence results if available
            if sequence_state.classification_data and "sequence_results" in blast_dict:
                for seq_key in blast_dict["sequence_results"]:
                    blast_dict["sequence_results"][seq_key]["classification"] = (
                        sequence_state.classification_data.to_dict()
                    )

            return blast_dict
        else:
            logger.debug("to_blast_data_dict: no sequence_state or blast_data")
        return None

    def to_workflow_state_dict(
        self, sequence_id: str = None
    ) -> Optional[Dict[str, Any]]:
        """Convert to legacy WorkflowState format for backward compatibility"""
        if sequence_id is None:
            sequence_id = self._active_sequence_id

        sequence_state = self.get_sequence(sequence_id)
        if sequence_state and sequence_state.workflow_state:
            return sequence_state.workflow_state.to_dict()
        return None


# Global instance - single source of truth (used when no session submission_id)
_legacy_pipeline_state = PipelineState()


def get_legacy_pipeline_state() -> PipelineState:
    """
    Get the pipeline state for the current request.

    When the request has a session with blast_submission_id, state is loaded
    from/saved to the filesystem cache so all workers share the same state
    (fixes multi-worker deployments where text-input submission and polling
    hit different workers).
    """
    try:
        from flask import session
    except ImportError:
        return _legacy_pipeline_state

    submission_id = session.get("blast_submission_id") if session else None
    if not submission_id:
        return _legacy_pipeline_state

    cache_key = f"{_BLAST_PIPELINE_CACHE_KEY_PREFIX}{submission_id}"
    try:
        data = cache.get(cache_key)
        if data:
            state = PipelineState.from_dict(data)
            state._cache_key = cache_key
            state._cache_timeout = DEFAULT_CACHE_TIMEOUT
            return state
        # New submission: return empty state that will be persisted when we mutate
        state = PipelineState()
        state._cache_key = cache_key
        state._cache_timeout = DEFAULT_CACHE_TIMEOUT
        return state
    except Exception as e:
        logger.warning("Failed to load pipeline state from cache: %s", e)
        return _legacy_pipeline_state


def reset_legacy_pipeline_state():
    """Reset the global legacy pipeline state (useful for testing)"""
    global _legacy_pipeline_state
    _legacy_pipeline_state = PipelineState()


class ConsolidatedPipelineState:
    """
    Simplified pipeline state manager using consolidated data models.

    This eliminates data duplication and provides a much cleaner API.
    """

    def __init__(self):
        self._analyses: Dict[str, SequenceAnalysis] = {}
        self._active_sequence_id: Optional[str] = None

    # Core sequence management
    def add_sequence(
        self, sequence_id: str, sequence: str = None, header: str = None
    ) -> SequenceAnalysis:
        """Add a new sequence for analysis"""
        if sequence_id in self._analyses:
            logger.warning(f"Sequence {sequence_id} already exists, returning existing")
            return self._analyses[sequence_id]

        analysis = SequenceAnalysis(
            sequence_id=sequence_id, sequence=sequence, sequence_header=header
        )
        self._analyses[sequence_id] = analysis

        # Set as active if it's the first sequence
        if self._active_sequence_id is None:
            self._active_sequence_id = sequence_id

        logger.debug(f"Added sequence {sequence_id} to pipeline state")
        return analysis

    def get_analysis(self, sequence_id: str) -> Optional[SequenceAnalysis]:
        """Get sequence analysis by ID"""
        return self._analyses.get(sequence_id)

    def get_active_analysis(self) -> Optional[SequenceAnalysis]:
        """Get the currently active sequence analysis"""
        if self._active_sequence_id:
            return self._analyses.get(self._active_sequence_id)
        return None

    def set_active_sequence(self, sequence_id: str):
        """Set the active sequence"""
        if sequence_id in self._analyses:
            self._active_sequence_id = sequence_id
            logger.debug(f"Set active sequence to {sequence_id}")
        else:
            logger.error(f"Cannot set active sequence to {sequence_id} - not found")

    # Status and progress management
    def update_status(self, sequence_id: str, status: WorkflowStatus):
        """Update workflow status for a sequence"""
        analysis = self.get_analysis(sequence_id)
        if analysis:
            analysis.status = status
            logger.debug(f"Updated status for {sequence_id}: {status.value}")

    def set_error(self, sequence_id: str, error: str):
        """Set error for a sequence"""
        analysis = self.get_analysis(sequence_id)
        if analysis:
            analysis.set_error(error)
            logger.error(f"Set error for {sequence_id}: {error}")

    def mark_complete(self, sequence_id: str):
        """Mark sequence analysis as complete"""
        analysis = self.get_analysis(sequence_id)
        if analysis:
            analysis.set_complete()
            logger.debug(f"Marked {sequence_id} as complete")

    # Query methods
    def get_completed_sequences(self) -> List[str]:
        """Get list of completed sequence IDs"""
        return [
            seq_id
            for seq_id, analysis in self._analyses.items()
            if analysis.is_complete()
        ]

    def get_failed_sequences(self) -> List[str]:
        """Get list of failed sequence IDs"""
        return [
            seq_id
            for seq_id, analysis in self._analyses.items()
            if analysis.has_error()
        ]

    def get_sequences_with_classification(self) -> List[str]:
        """Get list of sequence IDs that have classification results"""
        return [
            seq_id
            for seq_id, analysis in self._analyses.items()
            if analysis.has_classification()
        ]

    # Serialization
    def to_dict(self) -> Dict[str, any]:
        """Convert entire state to dictionary"""
        return {
            "analyses": {
                seq_id: analysis.to_dict()
                for seq_id, analysis in self._analyses.items()
            },
            "active_sequence_id": self._active_sequence_id,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, any]) -> "ConsolidatedPipelineState":
        """Create from dictionary"""
        state = cls()
        if "analyses" in data:
            for seq_id, analysis_data in data["analyses"].items():
                state._analyses[seq_id] = SequenceAnalysis.from_dict(analysis_data)
        state._active_sequence_id = data.get("active_sequence_id")
        return state

    # Legacy compatibility methods for backward compatibility with Dash stores
    def to_legacy_blast_data_dict(
        self, sequence_id: str = None
    ) -> Optional[Dict[str, any]]:
        """Convert to legacy BlastData format for backward compatibility"""
        if sequence_id is None:
            sequence_id = self._active_sequence_id

        analysis = self.get_analysis(sequence_id)
        if not analysis or not analysis.blast_result:
            return None

        blast_result = analysis.blast_result

        # Create legacy format
        legacy_dict = {
            "seq_type": blast_result.sequence_type.value,
            "sequence": blast_result.sequence,
            "fasta_file": blast_result.fasta_file,
            "blast_content": blast_result.blast_content,
            "blast_file": blast_result.blast_file,
            "blast_df": blast_result.blast_hits,
            "processed": blast_result.processed,
            "error": blast_result.error,
            # Legacy multi-sequence fields
            "processed_sequences": [0] if blast_result.processed else [],
            "sequence_results": {
                "0": {
                    "sequence": blast_result.sequence,
                    "blast_content": blast_result.blast_content,
                    "blast_file": blast_result.blast_file,
                    "blast_df": blast_result.blast_hits,
                    "processed": blast_result.processed,
                    "error": blast_result.error,
                }
            },
            "total_sequences": 1,
        }

        # Add classification data if available
        if analysis.classification:
            classification_dict = {
                "seq_type": analysis.classification.sequence_type.value,
                "source": analysis.classification.stage.value
                if analysis.classification.stage
                else None,
                "family": analysis.classification.family,
                "navis": analysis.classification.navis,
                "haplotype": analysis.classification.haplotype,
                "closest_match": analysis.classification.closest_match,
                "match_details": analysis.classification.match_details,
                "confidence": analysis.classification.confidence.value
                if analysis.classification.confidence
                else None,
            }
            legacy_dict["sequence_results"]["0"]["classification"] = classification_dict

        return legacy_dict

    def to_legacy_workflow_state_dict(
        self, sequence_id: str = None
    ) -> Optional[Dict[str, any]]:
        """Convert to legacy WorkflowState format for backward compatibility"""
        if sequence_id is None:
            sequence_id = self._active_sequence_id

        analysis = self.get_analysis(sequence_id)
        if not analysis:
            return None

        # Create legacy format
        legacy_dict = {
            "complete": analysis.is_complete(),
            "error": analysis.error,
            "found_match": analysis.has_classification(),
            "match_stage": analysis.classification.stage.value
            if analysis.classification and analysis.classification.stage
            else None,
            "match_result": analysis.classification.closest_match
            if analysis.classification
            else None,
            "classification_data": analysis.classification.to_dict()
            if analysis.classification
            else None,
            "stages": analysis.stage_progress,
            "task_id": analysis.sequence_id,
            "status": analysis.status.value,
            "workflow_started": analysis.status != WorkflowStatus.INITIALIZED,
            "current_stage": analysis.current_stage,
            "current_stage_idx": 0,  # Simplified for legacy compatibility
            "fetch_ship_params": {
                "curated": analysis.config.ship_curated_only,
                "with_sequence": analysis.config.ship_include_sequence,
                "dereplicate": analysis.config.ship_dereplicate,
            },
            "fetch_captain_params": {
                "curated": analysis.config.captain_curated_only,
                "with_sequence": analysis.config.captain_include_sequence,
            },
        }

        return legacy_dict

    def to_legacy_classification_data_dict(
        self, sequence_id: str = None
    ) -> Optional[Dict[str, any]]:
        """Convert to legacy ClassificationData format for backward compatibility"""
        if sequence_id is None:
            sequence_id = self._active_sequence_id

        analysis = self.get_analysis(sequence_id)
        if not analysis or not analysis.classification:
            return None

        classification = analysis.classification
        return {
            "seq_type": classification.sequence_type.value,
            "source": classification.stage.value if classification.stage else None,
            "family": classification.family,
            "navis": classification.navis,
            "haplotype": classification.haplotype,
            "closest_match": classification.closest_match,
            "match_details": classification.match_details,
            "confidence": classification.confidence.value
            if classification.confidence
            else None,
            "fasta_file": analysis.blast_result.fasta_file
            if analysis.blast_result
            else None,
        }


# Global instance - single source of truth
_consolidated_pipeline_state = ConsolidatedPipelineState()


def get_consolidated_pipeline_state() -> ConsolidatedPipelineState:
    """Get the global consolidated pipeline state instance"""
    return _consolidated_pipeline_state


def reset_consolidated_pipeline_state():
    """Reset the global consolidated pipeline state (useful for testing)"""
    global _consolidated_pipeline_state
    _consolidated_pipeline_state = ConsolidatedPipelineState()


class DashStateAdapter:
    """
    Adapter that connects Dash stores to the centralized pipeline state.

    This eliminates the need for multiple stores and ensures all data
    comes from a single source of truth.
    """

    def __init__(self):
        self.pipeline_state = get_legacy_pipeline_state()

    def get_blast_data_store(self, sequence_id: str = None) -> Optional[Dict[str, Any]]:
        """Get data for blast-data-store"""
        result = self.pipeline_state.to_blast_data_dict(sequence_id)
        logger.debug(
            f"get_blast_data_store({sequence_id}) returned: {result is not None}"
        )
        return result

    def update_blast_data_store(
        self, sequence_id: str, blast_data_dict: Dict[str, Any]
    ):
        """Update from blast-data-store"""
        from src.utils.blast_data import BlastData

        if blast_data_dict:
            blast_data = BlastData.from_dict(blast_data_dict)
            self.pipeline_state.update_blast_data(sequence_id, blast_data)

    def get_workflow_state_store(
        self, sequence_id: str = None
    ) -> Optional[Dict[str, Any]]:
        """Get data for workflow-state-store"""
        return self.pipeline_state.to_workflow_state_dict(sequence_id)

    def update_workflow_state_store(
        self, sequence_id: str, workflow_state_dict: Dict[str, Any]
    ):
        """Update from workflow-state-store"""
        from src.utils.blast_data import WorkflowState

        if workflow_state_dict:
            workflow_state = WorkflowState.from_dict(workflow_state_dict)
            self.pipeline_state.update_workflow_state(sequence_id, workflow_state)

    def get_classification_data_store(
        self, sequence_id: str = None
    ) -> Optional[Dict[str, Any]]:
        """Get data for classification-data-store"""
        classification_data = self.pipeline_state.get_classification_data(
            sequence_id or self.pipeline_state._active_sequence_id
        )
        return classification_data.to_dict() if classification_data else None

    def update_classification_data_store(
        self, sequence_id: str, classification_data_dict: Dict[str, Any]
    ):
        """Update from classification-data-store"""
        from src.utils.blast_data import ClassificationData

        if classification_data_dict:
            classification_data = ClassificationData.from_dict(classification_data_dict)
            self.pipeline_state.update_classification_data(
                sequence_id, classification_data
            )

    def sync_all_stores(self, sequence_id: str) -> Dict[str, Any]:
        """
        Sync all stores from the centralized state.
        Returns a dict with all store data.
        """
        return {
            "blast_data": self.get_blast_data_store(sequence_id),
            "workflow_state": self.get_workflow_state_store(sequence_id),
            "classification_data": self.get_classification_data_store(sequence_id),
        }

    def get_sequence_classification_for_ui(
        self, sequence_id: str = None
    ) -> Optional[Dict[str, Any]]:
        """
        Get classification data specifically formatted for UI display.
        This replaces the complex logic in blast.py callbacks.
        """
        if sequence_id is None:
            sequence_id = self.pipeline_state._active_sequence_id

        if not sequence_id:
            return None

        sequence_state = self.pipeline_state.get_sequence(sequence_id)
        if not sequence_state:
            return None

        # Always prefer the centralized classification data
        if sequence_state.classification_data:
            classification_dict = sequence_state.classification_data.to_dict()
            logger.debug(
                f"Returning classification data for UI (sequence {sequence_id}): {classification_dict}"
            )
            return classification_dict

        # Fallback to workflow state classification data if available
        if (
            sequence_state.workflow_state
            and sequence_state.workflow_state.classification_data
        ):
            logger.debug("Falling back to workflow state classification data")
            return sequence_state.workflow_state.classification_data

        logger.debug(f"No classification data found for sequence {sequence_id}")
        return None


# Global adapter instance
_dash_adapter = DashStateAdapter()


def get_dash_adapter() -> DashStateAdapter:
    """Get the global Dash state adapter"""
    return _dash_adapter


class MigrationAdapter:
    """
    Adapter that provides backward compatibility during migration.

    This allows existing code to continue working while gradually migrating
    to the consolidated data models.
    """

    def __init__(self):
        self.consolidated_state = get_consolidated_pipeline_state()

    # New consolidated API methods
    def add_sequence_analysis(
        self, sequence_id: str, sequence: str = None, header: str = None
    ) -> SequenceAnalysis:
        """Add a new sequence using the consolidated model"""
        return self.consolidated_state.add_sequence(sequence_id, sequence, header)

    def get_sequence_analysis(self, sequence_id: str) -> Optional[SequenceAnalysis]:
        """Get sequence analysis using the consolidated model"""
        return self.consolidated_state.get_analysis(sequence_id)

    def update_blast_result(self, sequence_id: str, blast_result: BlastResult):
        """Update BLAST result using consolidated model"""
        analysis = self.consolidated_state.get_analysis(sequence_id)
        if analysis:
            analysis.blast_result = blast_result
            logger.debug(f"Updated BLAST result for {sequence_id}")

    def update_classification_result(
        self, sequence_id: str, classification: ClassificationResult
    ):
        """Update classification result using consolidated model"""
        analysis = self.consolidated_state.get_analysis(sequence_id)
        if analysis:
            analysis.classification = classification
            logger.debug(f"Updated classification for {sequence_id}")

    # Legacy compatibility methods (for existing code)
    def get_blast_data_store(self, sequence_id: str = None) -> Optional[Dict[str, Any]]:
        """Get data for blast-data-store (legacy compatibility)"""
        return self.consolidated_state.to_legacy_blast_data_dict(sequence_id)

    def update_blast_data_store(
        self, sequence_id: str, blast_data_dict: Dict[str, Any]
    ):
        """Update from blast-data-store (legacy compatibility)"""
        if blast_data_dict:
            # Convert legacy dict to new model
            blast_result = convert_from_legacy_blast_data(blast_data_dict)
            if blast_result:
                self.update_blast_result(sequence_id, blast_result)

    def get_workflow_state_store(
        self, sequence_id: str = None
    ) -> Optional[Dict[str, Any]]:
        """Get data for workflow-state-store (legacy compatibility)"""
        return self.consolidated_state.to_legacy_workflow_state_dict(sequence_id)

    def update_workflow_state_store(
        self, sequence_id: str, workflow_state_dict: Dict[str, Any]
    ):
        """Update from workflow-state-store (legacy compatibility)"""
        if workflow_state_dict:
            analysis = self.consolidated_state.get_analysis(sequence_id)
            if analysis:
                # Update status
                status_map = {
                    "initialized": WorkflowStatus.INITIALIZED,
                    "running": WorkflowStatus.RUNNING,
                    "complete": WorkflowStatus.COMPLETE,
                    "failed": WorkflowStatus.FAILED,
                }
                if "status" in workflow_state_dict:
                    analysis.status = status_map.get(
                        workflow_state_dict["status"], WorkflowStatus.INITIALIZED
                    )

                # Update error
                if "error" in workflow_state_dict:
                    analysis.error = workflow_state_dict["error"]

                # Update current stage
                if "current_stage" in workflow_state_dict:
                    analysis.current_stage = workflow_state_dict["current_stage"]

                # Update stage progress
                if "stages" in workflow_state_dict:
                    analysis.stage_progress = workflow_state_dict["stages"]

    def get_classification_data_store(
        self, sequence_id: str = None
    ) -> Optional[Dict[str, Any]]:
        """Get data for classification-data-store (legacy compatibility)"""
        return self.consolidated_state.to_legacy_classification_data_dict(sequence_id)

    def update_classification_data_store(
        self, sequence_id: str, classification_data_dict: Dict[str, Any]
    ):
        """Update from classification-data-store (legacy compatibility)"""
        if classification_data_dict:
            # Convert legacy dict to new model
            classification = convert_from_legacy_classification_data(
                classification_data_dict
            )
            if classification:
                self.update_classification_result(sequence_id, classification)

    def sync_all_stores(self, sequence_id: str) -> Dict[str, Any]:
        """
        Sync all stores from the consolidated state (legacy compatibility).
        Returns a dict with all store data.
        """
        return {
            "blast_data": self.get_blast_data_store(sequence_id),
            "workflow_state": self.get_workflow_state_store(sequence_id),
            "classification_data": self.get_classification_data_store(sequence_id),
        }

    def get_sequence_classification_for_ui(
        self, sequence_id: str = None
    ) -> Optional[Dict[str, Any]]:
        """
        Get classification data specifically formatted for UI display (legacy compatibility).
        """
        if sequence_id is None:
            sequence_id = self.consolidated_state._active_sequence_id

        if not sequence_id:
            return None

        analysis = self.consolidated_state.get_analysis(sequence_id)
        if not analysis or not analysis.classification:
            return None

        classification = analysis.classification
        classification_dict = {
            "seq_type": classification.sequence_type.value,
            "source": classification.stage.value if classification.stage else None,
            "family": classification.family,
            "navis": classification.navis,
            "haplotype": classification.haplotype,
            "closest_match": classification.closest_match,
            "match_details": classification.match_details,
            "confidence": classification.confidence.value
            if classification.confidence
            else None,
            "fasta_file": analysis.blast_result.fasta_file
            if analysis.blast_result
            else None,
        }

        logger.debug(f"Returning classification data for UI: {classification_dict}")
        return classification_dict

    # Convenience methods for common operations
    def create_sequence_from_legacy_data(
        self,
        sequence_id: str,
        blast_data_dict: Dict[str, Any] = None,
        classification_data_dict: Dict[str, Any] = None,
    ) -> SequenceAnalysis:
        """Create a new sequence analysis from legacy data structures"""

        # Extract basic sequence info
        sequence = None
        header = None
        if blast_data_dict:
            sequence = blast_data_dict.get("sequence")
            # Try to extract header from fasta_file if available
            if blast_data_dict.get("fasta_file"):
                try:
                    with open(blast_data_dict["fasta_file"], "r") as f:
                        first_line = f.readline().strip()
                        if first_line.startswith(">"):
                            header = first_line[1:]  # Remove '>'
                except Exception as e:
                    logger.error(f"Error extracting header from fasta file: {e}")
                    raise e

        # Create the analysis
        analysis = self.add_sequence_analysis(sequence_id, sequence, header)

        # Convert and add BLAST data
        if blast_data_dict:
            blast_result = convert_from_legacy_blast_data(blast_data_dict)
            if blast_result:
                self.update_blast_result(sequence_id, blast_result)

        # Convert and add classification data
        if classification_data_dict:
            classification = convert_from_legacy_classification_data(
                classification_data_dict
            )
            if classification:
                self.update_classification_result(sequence_id, classification)

        return analysis


# Global adapter instance
_migration_adapter = MigrationAdapter()


def get_migration_adapter() -> MigrationAdapter:
    """Get the global migration adapter instance"""
    return _migration_adapter


# =============================================================================
# MIGRATION UTILITIES
# =============================================================================


def convert_sequence_analysis_to_legacy_blast_data(
    analysis: SequenceAnalysis,
) -> Dict[str, Any]:
    """
    Convert a SequenceAnalysis object to legacy BlastData format.

    This utility helps during the migration period by allowing new code that uses
    SequenceAnalysis to be compatible with existing code that expects BlastData format.
    """
    if not analysis:
        return None

    result = {
        "seq_type": analysis.sequence_type.value,
        "sequence": analysis.sequence,
        "fasta_file": analysis.blast_result.fasta_file
        if analysis.blast_result
        else None,
        "blast_content": analysis.blast_result.blast_content
        if analysis.blast_result
        else None,
        "blast_file": analysis.blast_result.blast_file
        if analysis.blast_result
        else None,
        "blast_df": analysis.blast_result.blast_hits if analysis.blast_result else None,
        "processed": analysis.is_complete(),
        "error": analysis.error,
        "processed_sequences": [0] if analysis.is_complete() else [],
        "total_sequences": 1,
        "sequence_results": {
            "0": {
                "sequence": analysis.sequence,
                "blast_content": analysis.blast_result.blast_content
                if analysis.blast_result
                else None,
                "blast_file": analysis.blast_result.blast_file
                if analysis.blast_result
                else None,
                "blast_df": analysis.blast_result.blast_hits
                if analysis.blast_result
                else None,
                "processed": analysis.is_complete(),
                "error": analysis.error,
            }
        },
    }

    # Add classification data if available
    if analysis.classification:
        classification_dict = {
            "seq_type": analysis.classification.sequence_type.value,
            "source": analysis.classification.stage.value
            if analysis.classification.stage
            else None,
            "family": analysis.classification.family,
            "navis": analysis.classification.navis,
            "haplotype": analysis.classification.haplotype,
            "closest_match": analysis.classification.closest_match,
            "match_details": analysis.classification.match_details,
            "confidence": analysis.classification.confidence.value
            if analysis.classification.confidence
            else None,
            "fasta_file": analysis.blast_result.fasta_file
            if analysis.blast_result
            else None,
        }
        result["sequence_results"]["0"]["classification"] = classification_dict

    return result


def convert_sequence_analysis_to_legacy_workflow_state(
    analysis: SequenceAnalysis,
) -> Dict[str, Any]:
    """
    Convert a SequenceAnalysis object to legacy WorkflowState format.
    """
    if not analysis:
        return None

    return {
        "complete": analysis.is_complete(),
        "error": analysis.error,
        "found_match": analysis.has_classification(),
        "match_stage": analysis.classification.stage.value
        if analysis.classification and analysis.classification.stage
        else None,
        "match_result": analysis.classification.closest_match
        if analysis.classification
        else None,
        "classification_data": analysis.classification.to_dict()
        if analysis.classification
        else None,
        "stages": analysis.stage_progress,
        "task_id": analysis.sequence_id,
        "status": analysis.status.value,
        "workflow_started": analysis.status != WorkflowStatus.INITIALIZED,
        "current_stage": analysis.current_stage,
        "current_stage_idx": 0,
        "fetch_ship_params": {
            "curated": analysis.config.ship_curated_only,
            "with_sequence": analysis.config.ship_include_sequence,
            "dereplicate": analysis.config.ship_dereplicate,
        },
        "fetch_captain_params": {
            "curated": analysis.config.captain_curated_only,
            "with_sequence": analysis.config.captain_include_sequence,
        },
    }


def enable_unified_processing():
    """Enable using the new unified data models globally"""
    global _use_unified_processing
    _use_unified_processing = True
    logger.info("Enabled unified processing with consolidated data models")


def disable_unified_processing():
    """Disable using the new unified data models globally"""
    global _use_unified_processing
    _use_unified_processing = False
    logger.info("Disabled unified processing, using legacy data models")


def is_unified_processing_enabled():
    """Check if unified processing is enabled"""
    return _use_unified_processing


# Global flag for migration - ENABLED by default for complete migration
_use_unified_processing = True


def safe_convert_sequence_analysis_to_legacy(
    analysis: SequenceAnalysis, tab_idx: int = 0
) -> Optional[Dict[str, Any]]:
    """
    Safely convert SequenceAnalysis to legacy format with proper error handling.

    Args:
        analysis: The SequenceAnalysis object to convert
        tab_idx: The tab index for multi-sequence results (default: 0)

    Returns:
        Dictionary in legacy BlastData format, or None if conversion fails
    """
    if not analysis:
        logger.warning("Cannot convert None analysis to legacy format")
        return None

    try:
        # Convert to legacy format
        legacy_dict = convert_sequence_analysis_to_legacy_blast_data(analysis)

        if not legacy_dict:
            logger.warning("Conversion to legacy format returned None")
            return None

        # Update the sequence_results key to match the tab index if needed
        if tab_idx != 0 and "sequence_results" in legacy_dict:
            # Move from "0" to str(tab_idx)
            seq_data = legacy_dict["sequence_results"]["0"]
            legacy_dict["sequence_results"] = {str(tab_idx): seq_data}
            legacy_dict["processed_sequences"] = (
                [tab_idx] if analysis.is_complete() else []
            )

        logger.debug(
            f"Successfully converted SequenceAnalysis to legacy format for tab {tab_idx}"
        )
        return legacy_dict

    except Exception as e:
        logger.error(f"Error converting SequenceAnalysis to legacy format: {e}")
        return None
