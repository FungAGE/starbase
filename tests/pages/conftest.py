import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock


@pytest.fixture
def sample_meta_dataframe():
    """Sample metadata DataFrame for testing"""
    return pd.DataFrame({
        "accession_tag": ["SSA000001", "SSA000002", "SSA000003", "SSA000004"],
        "familyName": ["Family1", "Family1", "Family2", "Family2"],
        "curated_status": ["curated", "uncurated", "curated", "uncurated"],
        "elementLength": [1000, 2000, 1500, 2500],
        "upDR": ["ATGC", "ATGC", "TTAA", "TTAA"],
        "downDR": ["GCAT", "GCAT", "AATT", "AATT"],
        "name": ["Species1", "Species2", "Species3", "Species4"],
        "order": ["Order1", "Order1", "Order2", "Order2"],
        "family": ["Family1", "Family1", "Family2", "Family2"],
        "genus": ["Genus1", "Genus1", "Genus2", "Genus2"],
        "subkingdom": ["Subkingdom1", "Subkingdom1", "Subkingdom2", "Subkingdom2"],
        "phylum": ["Phylum1", "Phylum1", "Phylum2", "Phylum2"],
        "subphylum": ["Subphylum1", "Subphylum1", "Subphylum2", "Subphylum2"],
        "class": ["Class1", "Class1", "Class2", "Class2"],
        "subclass": ["Subclass1", "Subclass1", "Subclass2", "Subclass2"],
        "suborder": ["Suborder1", "Suborder1", "Suborder2", "Suborder2"],
        "contigID": ["contig_1", "contig_2", "contig_3", "contig_4"],
        "ship_id": [1, 2, 3, 4],
    })


@pytest.fixture
def sample_papers_dataframe():
    """Sample papers DataFrame for testing"""
    return pd.DataFrame({
        "familyName": ["Family1", "Family2"],
        "type_element_reference": ["Ref1", "Ref2"],
        "Url": ["http://example.com/1", "http://example.com/2"],
        "shortCitation": ["Author1 et al.", "Author2 et al."],
    })


@pytest.fixture
def sample_meta_dict(sample_meta_dataframe):
    """Sample metadata as dictionary list for testing"""
    return sample_meta_dataframe.to_dict("records")


@pytest.fixture
def sample_papers_dict(sample_papers_dataframe):
    """Sample papers data as dictionary list for testing"""
    return sample_papers_dataframe.to_dict("records")


@pytest.fixture
def sample_download_rows():
    """Sample rows for download testing"""
    return [
        {"accession_tag": "SSA000001", "familyName": "Family1"},
        {"accession_tag": "SSA000002", "familyName": "Family2"},
    ]


@pytest.fixture
def sample_ships_dataframe():
    """Sample ships DataFrame for download testing"""
    return pd.DataFrame({
        "accession_tag": ["SSA000001", "SSA000002"],
        "familyName": ["Family1", "Family2"],
        "sequence": ["ATGCATGCATGC", "TTGGTTGGTTGG"],
        "name": ["Species1", "Species2"],
        "order": ["Order1", "Order2"],
        "family": ["Family1", "Family2"],
        "genus": ["Genus1", "Genus2"],
        "species": ["species1", "species2"],
        "strain": ["strain1", "strain2"],
        "taxID": [12345, 67890],
        "assembly_accession": ["GCA_000001", "GCA_000002"],
        "genomeSource": ["Source1", "Source2"],
        "contigID": ["Contig1", "Contig2"],
        "elementBegin": [100, 200],
        "elementEnd": [500, 600],
        "size": [400, 400],
    })


@pytest.fixture
def mock_cache():
    """Mock cache object for testing"""
    cache = MagicMock()
    cache.get.return_value = None
    cache.set.return_value = None
    return cache


@pytest.fixture
def mock_fetch_meta_data():
    """Mock fetch_meta_data function"""
    with pytest.MonkeyPatch().context() as m:
        mock_func = MagicMock()
        m.setattr("src.pages.wiki.fetch_meta_data", mock_func)
        yield mock_func


@pytest.fixture
def mock_fetch_paper_data():
    """Mock fetch_paper_data function"""
    with pytest.MonkeyPatch().context() as m:
        mock_func = MagicMock()
        m.setattr("src.pages.wiki.fetch_paper_data", mock_func)
        yield mock_func


@pytest.fixture
def mock_fetch_ships():
    """Mock fetch_ships function"""
    with pytest.MonkeyPatch().context() as m:
        mock_func = MagicMock()
        m.setattr("src.pages.wiki.fetch_ships", mock_func)
        yield mock_func


@pytest.fixture
def mock_make_dl_table():
    """Mock make_dl_table function"""
    with pytest.MonkeyPatch().context() as m:
        mock_func = MagicMock()
        m.setattr("src.pages.wiki.make_dl_table", mock_func)
        yield mock_func


@pytest.fixture
def mock_create_sunburst_plot():
    """Mock create_sunburst_plot function"""
    with pytest.MonkeyPatch().context() as m:
        mock_func = MagicMock()
        m.setattr("src.pages.wiki.create_sunburst_plot", mock_func)
        yield mock_func


@pytest.fixture
def mock_clean_contig_ids():
    """Mock clean_contigIDs function"""
    with pytest.MonkeyPatch().context() as m:
        mock_func = MagicMock(side_effect=lambda x: x.replace("contig_", "clean_"))
        m.setattr("src.pages.wiki.clean_contigIDs", mock_func)
        yield mock_func


@pytest.fixture
def mock_create_ncbi_style_header():
    """Mock create_ncbi_style_header function"""
    with pytest.MonkeyPatch().context() as m:
        mock_func = MagicMock(side_effect=lambda row: f">{row['accession_tag']} [family={row['familyName']}]")
        m.setattr("src.pages.wiki.create_ncbi_style_header", mock_func)
        yield mock_func


@pytest.fixture
def empty_dataframe():
    """Empty DataFrame for testing edge cases"""
    return pd.DataFrame()


@pytest.fixture
def dataframe_with_nans():
    """DataFrame with NaN values for testing"""
    return pd.DataFrame({
        "accession_tag": ["SSA000001", None, "SSA000003"],
        "familyName": ["Family1", "Family2", None],
        "elementLength": [1000, np.nan, 1500],
        "upDR": ["ATGC", None, "TTAA"],
        "downDR": ["GCAT", "GCAT", None],
    })


@pytest.fixture
def large_dataframe():
    """Large DataFrame for performance testing"""
    n_rows = 1000
    return pd.DataFrame({
        "accession_tag": [f"SSA{i:06d}" for i in range(n_rows)],
        "familyName": [f"Family{i % 10}" for i in range(n_rows)],
        "curated_status": ["curated" if i % 2 == 0 else "uncurated" for i in range(n_rows)],
        "elementLength": [1000 + i for i in range(n_rows)],
        "name": [f"Species{i}" for i in range(n_rows)],
        "order": [f"Order{i % 5}" for i in range(n_rows)],
        "family": [f"Family{i % 5}" for i in range(n_rows)],
        "genus": [f"Genus{i % 20}" for i in range(n_rows)],
    })
