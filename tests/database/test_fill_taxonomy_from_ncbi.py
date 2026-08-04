"""Tests for NCBI taxonomy backfill helpers."""

from unittest.mock import MagicMock, patch

from src.database.cleanup.utils.fill_taxonomy_from_ncbi import (
    ASSEMBLY_ACCESSION_RE,
    BIOSAMPLE_ACCESSION_RE,
    NcbiRateLimiter,
    _parse_organism,
    fetch_taxonomy_from_assembly,
)


def test_assembly_accession_pattern():
    assert ASSEMBLY_ACCESSION_RE.match("GCF_028828025.1")
    assert ASSEMBLY_ACCESSION_RE.match("GCA_000001305.1")
    assert not ASSEMBLY_ACCESSION_RE.match("Zymps1")


def test_biosample_accession_pattern():
    assert BIOSAMPLE_ACCESSION_RE.match("SAMN12345678")
    assert BIOSAMPLE_ACCESSION_RE.match("SAMEA1234567")
    assert not BIOSAMPLE_ACCESSION_RE.match("Zymps1")


def test_parse_organism():
    assert _parse_organism("Penicillium rubens") == ("Penicillium", "rubens")
    assert _parse_organism("Escherichia") == ("Escherichia", None)


def test_fetch_taxonomy_from_assembly_uses_cache():
    http = MagicMock()
    limiter = NcbiRateLimiter(min_interval=0)
    cache = {}

    search_resp = MagicMock()
    search_resp.json.return_value = {"esearchresult": {"idlist": ["123"]}}
    search_resp.raise_for_status = MagicMock()

    summary_resp = MagicMock()
    summary_resp.json.return_value = {
        "result": {
            "uids": ["123"],
            "123": {"taxid": "5075", "organism": "Penicillium rubens"},
        }
    }
    summary_resp.raise_for_status = MagicMock()

    http.get.side_effect = [search_resp, summary_resp]

    info = fetch_taxonomy_from_assembly(
        "GCF_028828025.1", http, cache, limiter, api_key="test-key"
    )
    assert info["taxid"] == "5075"
    assert info["genus"] == "Penicillium"
    assert info["species"] == "rubens"

    info2 = fetch_taxonomy_from_assembly(
        "GCF_028828025.1", http, cache, limiter, api_key="test-key"
    )
    assert info2 == info
    assert http.get.call_count == 2


def test_fill_taxonomy_from_ncbi_requires_api_key(monkeypatch):
    from src.database.cleanup.utils import fill_taxonomy_from_ncbi as mod

    monkeypatch.setattr(mod, "NCBI_API_KEY", None)
    try:
        mod.fill_taxonomy_from_ncbi(dry_run=True)
    except ValueError as exc:
        assert "NCBI_API_KEY" in str(exc)
    else:
        raise AssertionError("expected ValueError when NCBI_API_KEY is missing")


def test_fill_taxonomy_from_ncbi_dry_run(monkeypatch):
    from src.database.cleanup.utils import fill_taxonomy_from_ncbi as mod

    genome = MagicMock()
    genome.id = 1
    genome.assembly_accession = "GCF_028828025.1"
    genome.biosample = None
    genome.taxonomy_id = None
    genome.ome = "test"

    js = MagicMock()
    js.id = 10
    js.tax_id = None
    js.genome_id = 1
    js.starshipID = "Test_1"

    session = MagicMock()
    session.query.return_value.options.return_value.filter.return_value.all.return_value = [
        genome
    ]
    session.query.return_value.options.return_value.filter.return_value.filter.return_value.all.return_value = [
        js
    ]

    tax_row = MagicMock()
    tax_row.id = 99

    with (
        patch.object(mod, "StarbaseSession", return_value=session),
        patch.object(
            mod,
            "_resolve_genome_taxonomy_info",
            return_value={
                "taxid": "5075",
                "organism": "Penicillium rubens",
                "genus": "Penicillium",
                "species": "rubens",
                "lookup": "assembly",
                "lookup_key": "GCF_028828025.1",
            },
        ),
        patch.object(mod, "_find_taxonomy", return_value=tax_row),
    ):
        summary, rows, counts = mod.fill_taxonomy_from_ncbi(
            dry_run=True, overwrite=False, api_key="test-key"
        )

    assert "DRY RUN" in summary
    assert counts["genomes_updated"] == 1
    assert counts["joined_ships_updated"] == 1
    assert rows[0]["taxonomy_id"] == 99
    session.commit.assert_not_called()
