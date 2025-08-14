"""
Data models for the classification pipeline.

Overview
========
- the reason I've created dataclasses, is to make the flow of data through the blast/classification pipelines easier to handle, especially with dash callbacks.
- dataclasses are JSON-friendly and have to_dict/from_dict helpers so they can be used in dash stores and moved across process boundaries.

Key dataclasses
---------
- WorkflowState: holds the current state and outputs for the multi-stage classification workflow. 
  - Tracks progress and status for each stage
  - whether a match was found
  - the stage producing a match
  - attaches the blast and classiciation outputs
- This is the main dataclass that is used to generate UI components.
  
- BlastData: holds the inputs and outputs from blast
  - sequence/FASTA path
  - BLAST results (either content, file, or a  normalized list/dict form)
  - error messages flag
  - optional
  
- ClassificationData: A compact result payload describing what was inferred for
  the query (source of the decision, family/navis/haplotype when applicable,
  closest_match, match_details, and a coarse confidence level). This object can
  be shown directly in the UI.

- Dataclasses that hold parameters for methods used during classification
    - FetchShipParams, FetchCaptainParams
        - contain booleans
        - embedded inside BlastData


Serialization
-------------
All dataclasses expose to_dict()/from_dict() to ensure stable, explicit JSON
shapes. Prefer using these helpers instead of vars()/asdict() to keep control of
nested structures and backwards compatibility.
"""

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
            for field in ['seq_type', 'fasta_file', 'source', 'family', 'navis', 'haplotype', 'closest_match', 'match_details', 'confidence']
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
        fetch_captain_params = FetchCaptainParams.from_dict(data.get("fetch_captain_params", {}))

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
        self.classification_data = classification_data.to_dict() if classification_data else None