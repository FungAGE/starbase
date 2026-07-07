"""Tests for submission workflow refactor."""

import pandas as pd
import pytest

from src.utils.web_submission_adapter import (
    WebValidationError,
    encode_single_fasta,
    infer_metadata_strand,
    merge_accession_reference_pools,
    parse_location_metadata_tsv,
    parse_submission_fasta,
    submission_row_to_validated_data,
)
from src.components.submission_location_grid import rows_to_locations


def test_parse_submission_fasta_plain_text():
    fasta = ">ship1\nATGCATGC\n"
    seq_id, sequence = parse_submission_fasta(fasta)
    assert seq_id == "ship1"
    assert sequence == "ATGCATGC"


def test_parse_submission_fasta_base64_upload():
    encoded = encode_single_fasta("ship2", "TTTTAAAA")
    seq_id, sequence = parse_submission_fasta(encoded)
    assert seq_id == "ship2"
    assert sequence == "TTTTAAAA"


def test_merge_accession_reference_pools_deduplicates_by_accession():
    main_ships = pd.DataFrame(
        {
            "accession_tag": ["SSA000001", "SSA000002"],
            "accession_display": ["SSA000001", "SSA000002"],
            "sequence": ["ATGC", "GGCCTT"],
        }
    )
    staging_ships = pd.DataFrame(
        {
            "accession_tag": ["SSA000002", "SSA000003"],
            "accession_display": ["SSA000002", "SSA000003"],
            "sequence": ["GGCCTT", "AAAA"],
        }
    )
    combined = merge_accession_reference_pools(main_ships, staging_ships)
    assert len(combined) == 3
    assert set(combined["accession_tag"]) == {"SSA000001", "SSA000002", "SSA000003"}


def test_submission_row_to_validated_data_maps_fields():
    row = {
        "seq_contents": ">starshipA\nATGCATGCATGC\n",
        "seq_filename": "upload.fa",
        "uploader": "curator@example.com",
        "evidence": "manual",
        "genus": "Alternaria",
        "species": "alternata",
        "strain": "strain1",
        "hostchr": "scaffold_1",
        "shipstart": 100,
        "shipend": 5000,
        "shipstrand": "-",
        "assembly_accession": "GCA_000001305.1",
        "comment": "test",
        "classification_family": "FamA",
        "classification_navis": "NavA",
        "classification_haplotype": "HapA",
        "classification_source": "exact",
        "closest_match": "SSA000001",
        "classification_confidence": "high",
    }
    validated = submission_row_to_validated_data(row)
    assert validated["seq_id"] == "starshipA"
    assert validated["sequence"] == "ATGCATGCATGC"
    assert validated["strand_radio"] == 2
    assert validated["classification"]["closest_match"] == "SSA000001"


def test_infer_metadata_strand_from_coordinates():
    assert infer_metadata_strand(None, 5000, 1000) == "-"
    assert infer_metadata_strand(None, 1000, 5000) == "+"
    assert infer_metadata_strand("+", 5000, 1000) == "+"


def test_rows_to_locations_infers_strand():
    locations = rows_to_locations(
        [{"ship_header": "s1", "shipstart": "5000", "shipend": "1000", "strand": ""}]
    )
    assert locations[0]["strand_radio"] == 2


def test_parse_location_metadata_tsv_partial_columns():
    tsv = (
        "ship_header\tgenus\tspecies\tassembly_accession\thostchr\tshipstart\tshipend\n"
        "seq_a\tAlternaria\talternata\tGCA_000001305.1\tchr1\t5000\t1000\n"
    )
    rows, count = parse_location_metadata_tsv(tsv, ["seq_a", "seq_b"])
    assert count == 1
    assert rows[0]["genus"] == "Alternaria"
    assert rows[0]["assembly_accession"] == "GCA_000001305.1"
    assert rows[0]["strand"] == "-"
    assert rows[1]["ship_header"] == "seq_b"
    assert rows[1]["genus"] == ""


def test_parse_location_metadata_tsv_rejects_unknown_columns():
    tsv = "sequence,genus,species\nship1,Fusarium,oxysporum\n"
    with pytest.raises(WebValidationError, match="Unknown metadata column"):
        parse_location_metadata_tsv(tsv, ["ship1"])


@pytest.mark.integration
def test_insert_submission_sets_group_id():
    """Grouped uploads share submission_group_id and start as pending."""
    import datetime
    from src.database.sql_engine import get_submissions_session
    from src.database.models.schema import Submission
    from src.pages.submit import insert_submission

    group_id = "test-group-integration-001"
    with get_submissions_session() as session:
        try:
            row_id = insert_submission(
                seq_contents=">ship\nATGC\n",
                seq_filename="test.fa",
                seq_date=datetime.datetime.now().timestamp(),
                anno_contents=None,
                anno_filename=None,
                anno_date=None,
                uploader="test@example.com",
                evidence="manual",
                genus="Test",
                species="species",
                hostchr="chr1",
                shipstart=1,
                shipend=100,
                shipstrand="+",
                comment="workflow test",
                accession=None,
                needs_review=True,
                submission_group_id=group_id,
                processing_status="pending",
            )
            row = session.get(Submission, row_id)
            assert row is not None
            assert row.submission_group_id == group_id
            assert row.processing_status == "pending"
            assert row.accession_tag is None
        finally:
            session.query(Submission).filter(
                Submission.submission_group_id == group_id
            ).delete()
            session.commit()
