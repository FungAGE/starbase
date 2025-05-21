from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any


@dataclass
class FetchShipParams:
    curated: bool = False
    with_sequence: bool = True
    dereplicate: bool = True


@dataclass
class FetchCaptainParams:
    curated: bool = False
    with_sequence: bool = True


@dataclass
class ClassificationData:
    source: Optional[str] = None
    family: Optional[str] = None
    navis: Optional[str] = None
    haplotype: Optional[str] = None
    closest_match: Optional[str] = None
    match_details: Optional[str] = None
    confidence: Optional[str] = None


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
            "fetch_ship_params": {
                "curated": self.fetch_ship_params.curated,
                "with_sequence": self.fetch_ship_params.with_sequence,
                "dereplicate": self.fetch_ship_params.dereplicate,
            },
            "fetch_captain_params": {
                "curated": self.fetch_captain_params.curated,
                "with_sequence": self.fetch_captain_params.with_sequence,
            },
            "blast_df": self.blast_df,
            "blast_content": self.blast_content,
            "blast_file": self.blast_file,
            "sequence": self.sequence,
            "fasta_file": self.fasta_file,
            "error": self.error,
            "processed": self.processed,
        }

        if self.classification:
            result["classification"] = {
                "source": self.classification.source,
                "family": self.classification.family,
                "navis": self.classification.navis,
                "haplotype": self.classification.haplotype,
                "closest_match": self.classification.closest_match,
                "match_details": self.classification.match_details,
                "confidence": self.classification.confidence,
            }

        return result

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "BlastData":
        """Create a BlastData instance from a dictionary"""
        # Handle nested dictionaries
        fetch_ship_params = FetchShipParams(**data.get("fetch_ship_params", {}))
        fetch_captain_params = FetchCaptainParams(
            **data.get("fetch_captain_params", {})
        )

        classification_data = None
        if "classification" in data and data["classification"]:
            classification_data = ClassificationData(**data["classification"])

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
