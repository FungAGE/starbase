"""Tests for genome coordinate validation helpers."""

import pytest
import requests

from src.database.cleanup.utils.validate_genome_coordinates import (
    ContigAlignment,
    _contig_index_from_reports,
    _coords_already_match,
    _is_ncbi_assembly,
    _lookup_contig,
    _parse_minimap2_paf,
    _parse_submission_fasta,
    _process_ncbi_row,
    _sequence_matches_genomic,
    _validate_against_assembly,
    _validate_local_row,
    fetch_ncbi_contig_index,
    fetch_ncbi_sequence_region,
)

# Stable public reference for live NCBI checks (E. coli K-12 MG1655).
_LIVE_ASSEMBLY = "GCA_000005845.2"
_LIVE_CONTIG = "NC_000913.3"
_LIVE_BEGIN = 1000
_LIVE_END = 1099


def test_is_ncbi_assembly():
    assert _is_ncbi_assembly("GCA_000001305.1")
    assert _is_ncbi_assembly("GCF_019924465.1")
    assert not _is_ncbi_assembly("Zymps1")
    assert not _is_ncbi_assembly("")


def test_local_missing_genome_for_coords():
    row = {
        "joined_ship_id": 1,
        "starshipID": "SHIP1",
        "ship_id": 10,
        "genome_id": None,
        "assembly_accession": "",
        "genome_source": "",
        "contig_id": "NC_000001.11",
        "element_begin": 100,
        "element_end": 500,
        "element_length": 401,
    }
    issues = _validate_local_row(row, source="main")
    types = {i["issue_type"] for i in issues}
    assert "coords_without_genome" in types


def test_local_sequence_span_mismatch():
    row = {
        "joined_ship_id": 2,
        "starshipID": "SHIP2",
        "ship_id": 11,
        "genome_id": 5,
        "assembly_accession": "GCA_000001305.1",
        "genome_source": "ncbi",
        "contig_id": "NC_000001.11",
        "element_begin": 100,
        "element_end": 109,
        "element_length": 10,
        "ship_sequence": "A" * 50,
    }
    issues = _validate_local_row(row, source="main")
    assert any(i["issue_type"] == "sequence_span_mismatch" for i in issues)


def test_ncbi_contig_and_range_validation():
    index = _contig_index_from_reports(
        [
            {
                "accession": "NC_000001.11",
                "refseq_accession": "NC_000001.11",
                "length": 1000,
                "name": "1",
            }
        ]
    )
    info = _lookup_contig("NC_000001.11", index)
    assert info is not None
    assert info.length == 1000
    assert info.seq_accession == "NC_000001.11"

    row = {
        "joined_ship_id": 3,
        "starshipID": "SHIP3",
        "ship_id": 12,
        "assembly_accession": "GCA_000001305.1",
        "genome_source": "ncbi",
        "contig_id": "NC_000001.11",
        "element_begin": 900,
        "element_end": 1100,
    }
    issues = _validate_against_assembly(row, index, source="main")
    assert any(i["issue_type"] == "coordinates_out_of_range" for i in issues)

    row_ok = dict(row)
    row_ok["element_end"] = 950
    assert not _validate_against_assembly(row_ok, index, source="main")


def test_sequence_match_exact():
    seq = "ATGCATGCATGCATGC"
    matched, kind = _sequence_matches_genomic(seq, seq)
    assert matched
    assert kind == "exact"


def test_sequence_match_reverse_complement():
    from src.utils.seq_utils import revcomp

    seq = "ATGCATGCATGCATGC"
    matched, kind = _sequence_matches_genomic(str(revcomp(seq)), seq)
    assert matched
    assert kind == "exact"


def test_sequence_mismatch():
    matched, kind = _sequence_matches_genomic("ATGCATGC", "TTTTTTTT")
    assert not matched
    assert kind == ""


def test_parse_minimap2_paf_forward():
    fields = [
        "query",
        "16",
        "0",
        "16",
        "+",
        "contig",
        "100",
        "10",
        "26",
        "16",
        "16",
    ]
    alignment = _parse_minimap2_paf(fields, query_len=16)
    assert alignment is not None
    assert alignment.ref_begin == 11
    assert alignment.ref_end == 26
    assert alignment.strand == "+"
    assert alignment.coverage == 1.0
    assert alignment.identity == 1.0
    assert alignment.is_perfect


def test_contig_alignment_minus_strand_coords():
    alignment = ContigAlignment(
        ref_begin=100,
        ref_end=199,
        strand="-",
        coverage=1.0,
        identity=1.0,
    )
    assert alignment.element_begin == 199
    assert alignment.element_end == 100
    assert alignment.element_length == 100


def test_coords_already_match():
    alignment = ContigAlignment(
        ref_begin=1000,
        ref_end=1099,
        strand="+",
        coverage=1.0,
        identity=1.0,
    )
    assert _coords_already_match(1000, 1099, "+", alignment)
    assert not _coords_already_match(1000, 1099, "-", alignment)


def test_parse_submission_fasta():
    fasta = ">seq1\nATGC\nTTAA\n>seq2\nGGGG\n"
    assert _parse_submission_fasta(fasta) == "ATGCTTAA"


@pytest.mark.integration
def test_ncbi_full_sequence_validation_against_live_genbank():
    """Fetch a real interval from GenBank and verify full sequence match/mismatch."""
    session = requests.Session()
    ncbi_cache: dict = {}
    seq_cache: dict = {}

    index = fetch_ncbi_contig_index(_LIVE_ASSEMBLY, ncbi_cache, session)
    assert index is not None
    info = _lookup_contig(_LIVE_CONTIG, index)
    assert info is not None

    genomic = fetch_ncbi_sequence_region(
        info.seq_accession, _LIVE_BEGIN, _LIVE_END, seq_cache, session
    )
    span = _LIVE_END - _LIVE_BEGIN + 1
    assert genomic and len(genomic) == span

    base_row = {
        "joined_ship_id": 0,
        "starshipID": "LIVE_TEST",
        "ship_id": 0,
        "assembly_accession": _LIVE_ASSEMBLY,
        "genome_source": "ncbi",
        "contig_id": _LIVE_CONTIG,
        "element_begin": _LIVE_BEGIN,
        "element_end": _LIVE_END,
        "element_length": span,
        "strand": "+",
    }

    def run(ship_sequence: str):
        stats: dict = {}
        issues = _process_ncbi_row(
            {**base_row, "ship_sequence": ship_sequence},
            source="main",
            validate_ncbi=True,
            validate_sequences=True,
            ncbi_cache=ncbi_cache,
            seq_cache=seq_cache,
            session=session,
            stats=stats,
        )
        return issues, stats

    good_issues, good_stats = run(genomic)
    assert good_issues == []
    assert good_stats.get("sequence_matches") == 1

    bad_issues, bad_stats = run("A" * span)
    assert any(i["issue_type"] == "sequence_mismatch" for i in bad_issues)
    assert bad_stats.get("sequence_mismatches") == 1
