import pytest
import pandas as pd

from src.utils.classification_utils import generate_md5_hash
from src.utils.seq_utils import clean_sequence, revcomp


def _mock_ships_df():
    seq851 = "ATGCATGCATGC"
    seq904 = "ATGCATGCATGCATGC"
    seq596 = "ATGCATGCATTT"
    return pd.DataFrame(
        {
            "accession_tag": [
                "SSA002851",
                "SSA002904",
                "SSA002596",
                "SSA000001",
                "SSA000002",
            ],
            "accession_display": [
                "SSA002851.1",
                "SSA002904.1",
                "SSA002596.1",
                "SSA000001.1",
                "SSA000002.1",
            ],
            "sequence": [
                seq851,
                seq904,
                seq596,
                "ATGCATGCATGCATGCATGC",
                "TTTTTTTTTTTTTTTTTTTT",
            ],
            "familyName": ["Prometheus"] * 5,
            "md5": [
                generate_md5_hash(clean_sequence(s))
                for s in [
                    seq851,
                    seq904,
                    seq596,
                    "ATGCATGCATGCATGCATGC",
                    "TTTTTTTTTTTTTTTTTTTT",
                ]
            ],
            "rev_comp_md5": [
                generate_md5_hash(clean_sequence(revcomp(s)))
                for s in [
                    seq851,
                    seq904,
                    seq596,
                    "ATGCATGCATGCATGCATGC",
                    "TTTTTTTTTTTTTTTTTTTT",
                ]
            ],
        }
    )


@pytest.fixture
def test_ships_df():
    return _mock_ships_df()


@pytest.fixture
def test_sequence():
    return "ATGCATGCATGC"


@pytest.fixture
def test_sequence_revcomp(test_sequence):
    return revcomp(test_sequence)


@pytest.fixture
def test_contained_sequence():
    return "ATGCATGC"


@pytest.fixture
def test_similar_sequence():
    return "ATGCATGCATTT"


@pytest.fixture
def test_haplotype_ships_df():
    seq_a = "ATGCATGCATGCATGCATGC"
    seq_b = "ATGCATGCATGCATGCATGT"
    return pd.DataFrame(
        {
            "accession_tag": ["SSA002851", "SSA002904", "SSA002596"],
            "accession_display": ["SSA002851.1", "SSA002904.1", "SSA002596.1"],
            "sequence": [seq_a, seq_a, seq_b],
            "haplotype_name": ["2", "var22", "Ph1h1"],
            "captainID": ["captain_130", "captain_1285", "captain_1247"],
            "md5": [
                generate_md5_hash(clean_sequence(s)) for s in [seq_a, seq_a, seq_b]
            ],
            "rev_comp_md5": [
                generate_md5_hash(clean_sequence(revcomp(s)))
                for s in [seq_a, seq_a, seq_b]
            ],
        }
    )


@pytest.fixture
def test_haplotype_sequence(test_haplotype_ships_df):
    return test_haplotype_ships_df.iloc[0]["sequence"]


@pytest.fixture
def test_similarities():
    return [
        ("query_sequence", "SSA002851.1", 0.99),
        ("query_sequence", "SSA002904.1", 0.85),
        ("query_sequence", "SSA002596.1", 0.85),
        ("SSA002851.1", "SSA002904.1", 0.98),
    ]
