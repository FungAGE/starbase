from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any


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
    source: Optional[str] = None
    family: Optional[str] = None
    navis: Optional[str] = None
    haplotype: Optional[str] = None
    closest_match: Optional[str] = None
    match_details: Optional[str] = None
    confidence: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
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
            for field in ['source', 'family', 'navis', 'haplotype', 'closest_match', 'match_details', 'confidence']
        )


@dataclass
class BlastData:
    seq_type: str
    fetch_ship_params: FetchShipParams = field(default_factory=FetchShipParams)
    fetch_captain_params: FetchCaptainParams = field(default_factory=FetchCaptainParams)
    blast_df: Optional[List[Dict[str, Any]]] = None
    blast_content: Optional[str] = None
    blast_file: Optional[str] = None
    sequence: Optional[str] = None
    fasta_file: Optional[str] = None
    error: Optional[str] = None
    classification: Optional[ClassificationData] = None
    processed: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Convert the data class to a dictionary for backward compatibility"""
        result = {
            "seq_type": self.seq_type,
            "fetch_ship_params": self.fetch_ship_params.to_dict(),
            "fetch_captain_params": self.fetch_captain_params.to_dict(),
            "blast_df": self.blast_df,
            "blast_content": self.blast_content,
            "blast_file": self.blast_file,
            "sequence": self.sequence,
            "fasta_file": self.fasta_file,
            "error": self.error,
            "processed": self.processed,
        }

        if self.classification:
            result["classification"] = self.classification.to_dict()

        return result

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "BlastData":
        """Create a BlastData instance from a dictionary"""
        # Handle nested dictionaries
        fetch_ship_params = FetchShipParams.from_dict(data.get("fetch_ship_params", {}))
        fetch_captain_params = FetchCaptainParams.from_dict(data.get("fetch_captain_params", {}))

        classification_data = ClassificationData.from_dict(data.get("classification"))

        # Create the class instance with extracted data
        return cls(
            seq_type=data.get("seq_type", "nucl"),
            fetch_ship_params=fetch_ship_params,
            fetch_captain_params=fetch_captain_params,
            blast_df=data.get("blast_df"),
            blast_content=data.get("blast_content"),
            blast_file=data.get("blast_file"),
            sequence=data.get("sequence"),
            fasta_file=data.get("fasta_file"),
            error=data.get("error"),
            classification=classification_data,
            processed=data.get("processed", False),
        )


@dataclass
class WorkflowState:
    complete: bool = False
    error: Optional[str] = None
    found_match: bool = False
    match_stage: Optional[str] = None
    classification_data: Optional[ClassificationData] = None
    stages: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    task_id: str = ""
    status: str = "initialized"
    workflow_started: bool = True
    current_stage: Optional[str] = None
    current_stage_idx: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Convert WorkflowState to dictionary for JSON serialization."""
        return {
            "complete": self.complete,
            "error": self.error,
            "found_match": self.found_match,
            "match_stage": self.match_stage,
            "classification_data": self.classification_data.to_dict() if self.classification_data else None,
            "stages": self.stages,
            "task_id": self.task_id,
            "status": self.status,
            "workflow_started": self.workflow_started,
            "current_stage": self.current_stage,
            "current_stage_idx": self.current_stage_idx,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "WorkflowState":
        """Create a WorkflowState instance from a dictionary"""
        classification_data = ClassificationData.from_dict(data.get("classification_data"))
        
        return cls(
            complete=data.get("complete", False),
            error=data.get("error"),
            found_match=data.get("found_match", False),
            match_stage=data.get("match_stage"),
            classification_data=classification_data,
            stages=data.get("stages", {}),
            task_id=data.get("task_id", ""),
            status=data.get("status", "initialized"),
            workflow_started=data.get("workflow_started", True),
            current_stage=data.get("current_stage"),
            current_stage_idx=data.get("current_stage_idx", 0),
        )

    def set_classification(self, classification: ClassificationData) -> None:
        """Set the classification data and update related fields"""
        self.classification_data = classification
        if classification and classification.source:
            self.match_stage = classification.source
        if classification and not classification.is_empty():
            self.found_match = True

    def get_classification_dict(self) -> Optional[Dict[str, Any]]:
        """Get classification data as dictionary for backward compatibility"""
        if self.classification_data:
            return self.classification_data.to_dict()
        return None