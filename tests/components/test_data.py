from unittest.mock import patch
import pandas as pd
from src.database.sql_manager import fetch_meta_data, fetch_ships, fetch_captains
from src.components.data import create_ship_accession_modal_data, create_accession_modal_data

def test_successful_modal_data_creation(test_single_accession_meta_data):
    """Test successful creation of ship accession modal data"""
    # Mock fetch_meta_data to return test data
    with patch("src.components.data.fetch_meta_data") as mock_fetch:
        mock_fetch.return_value = test_single_accession_meta_data

        with patch("src.components.data.get_quality_tags") as mock_quality:
            mock_quality.return_value = [
                {"tag_type": "missing_direct_repeats", "tag_value": None}
            ]

            result = create_ship_accession_modal_data(test_single_accession_meta_data["ship_accession_tag"])

            assert isinstance(result, dict)

            assert result["title"] == test_single_accession_meta_data["ship_accession_tag"]
            assert result["familyName"] == "Prometheus"
            assert result["curated_status"] == "curated"


def test_successful_modal_data_creation(test_multiple_accession_meta_data):
    """Test successful creation of accession modal data"""
    # Mock fetch_meta_data to return test data
    with patch("src.components.data.fetch_meta_data") as mock_fetch:
        mock_fetch.return_value = test_multiple_accession_meta_data

        with patch("src.components.data.get_quality_tags") as mock_quality:
            mock_quality.return_value = []

            result = create_accession_modal_data("SSA002851.1")

            assert isinstance(result, dict)
            assert result["title"] == "Starship Accession: SSA002851.1"
            assert result["familyName"] == "Prometheus"
            assert result["genomes_present"] == "1"


def test_single_accession_data(single_accession_meta_data):
    """Legacy fixture providing raw DataFrame - kept for backward compatibility"""
    
    meta_df = fetch_meta_data(accessions=single_accession_meta_data["ship_accession_tag"])
    assert not meta_df.empty
    assert meta_df.equals(pd.DataFrame(single_accession_meta_data))


def test_multiple_accession_data(multiple_accession_meta_data):
    """Legacy fixture providing raw DataFrame - kept for backward compatibility"""
    meta_df = fetch_meta_data(accessions=multiple_accession_meta_data["ship_accession_tag"])
    assert not meta_df.empty
    assert meta_df.equals(pd.DataFrame(multiple_accession_meta_data))


def test_single_accession_multiple_ships(single_accession_multiple_ships_meta_data):
    meta_df = fetch_meta_data(accessions=single_accession_multiple_ships_meta_data["ship_accession_tag"])
    assert not meta_df.empty
    assert meta_df.equals(pd.DataFrame(single_accession_multiple_ships_meta_data))


def test_ship_accession_modal_data(single_accession_meta_data):
    """Fixture providing processed ship accession modal data dictionary"""
    from src.components.data import create_ship_accession_modal_data
    return create_ship_accession_modal_data(single_accession_meta_data["ship_accession_tag"])


def test_accession_modal_data(multiple_accession_meta_data):
    """Fixture providing processed accession modal data dictionary"""
    from src.components.data import create_accession_modal_data
    return create_accession_modal_data(multiple_accession_meta_data["ship_accession_tag"])


def test_ships_df(test_ships_df):
    try:
        ships_df = fetch_ships(
            accessions=["SSA002851", "SSA002904", "SSA002596"],
            with_sequence=True,
        )
        if not ships_df.empty:
            return ships_df
    except Exception:
        pass


()
def test_captains_df():
    try:
        # Try a few different accession tags that are more likely to work
        for accession in ["SSA002851", "SSA002904", "SSA002596"]:
            captains_df = fetch_captains(
                accessions=[accession], with_sequence=True
            )
            if not captains_df.empty:
                return captains_df
    except Exception:
        pass



def test_sequence():
    # looking for a real sequence from the database
    try:
        ship_df = fetch_ships(
            accessions=["SSA002851"], with_sequence=True
        )
        if not ship_df.empty:
            sequence = ship_df.iloc[0]["sequence"]
            return sequence
    except Exception:
        pass
    return "ATGCATGCATGC"  # Mock sequence



def test_sequence_revcomp(test_sequence):
    try:
        # return the reverse complement of the real sequence
        complement = str.maketrans("ATGC", "TACG")
        revcomp = test_sequence.translate(complement)[::-1]
        return revcomp
    except Exception:
        pass
    return "GCATGCATGCAT"  # Mock sequence



def test_contained_sequence(test_sequence):
    try:
        # return a contained subsequence of the real sequence
        return test_sequence[
            50:-50
        ]  # Take from position 50 to 50 positions from the end
    except Exception:
        pass
    return "ATGCATGC"  # Mock sequence



def test_similar_sequence(test_sequence):
    # introduce a small mutation to create a similar sequence
    try:
        return test_sequence[:10] + "A" + test_sequence[11:]  # Change one base
    except Exception:
        pass
    return "ATGCATGCAA"  # Mock sequence


# separate test fixtures for haplotype matching

def test_haplotype_ships_df():
    try:
        ships_df = fetch_ships(
            accessions=["SSA002851", "SSA002904", "SSA002596"],
            with_sequence=True,
        )
        if not ships_df.empty:
            return ships_df
    except Exception:
        pass

    return pd.DataFrame(
        {
            "accession_tag": ["SSA002851", "SSA002904", "SSA002596"],
            "accession_display": ["SSA002851.1", "SSA002904.1", "SSA002596.1"],
            "sequence": [
                "ATGCATGCATGCATGCATGC",  # Base sequence
                "ATGCATGCATGCATGCATGC",  # Identical sequence (should have 1.0 similarity)
                "ATGCATGCATGCATGCATGT",  # One base different (should have high similarity)
            ],
            "haplotype_name": ["2", "var22", "Ph1h1"],  # Add haplotype information
            "captainID": [
                "captain_130",
                "captain_1285",
                "captain_1247",
            ],  # Add captain IDs
            "md5": [
                generate_md5_hash(clean_sequence("ATGCATGCATGCATGCATGC")),
                generate_md5_hash(clean_sequence("ATGCATGCATGCATGCATGC")),
                generate_md5_hash(clean_sequence("ATGCATGCATGCATGCATGT")),
            ],
            "rev_comp_md5": [
                generate_md5_hash(clean_sequence("GCATGCATGCATGCATGCAT")),
                generate_md5_hash(clean_sequence("GCATGCATGCATGCATGCAT")),
                generate_md5_hash(clean_sequence("ACATGCATGCATGCATGCAT")),
            ],
        }
    )



def test_haplotype_sequence(test_haplotype_ships_df):
    # return the sequence with the haplotype
    try:
        sequence = test_haplotype_ships_df.iloc[0]["sequence"]  # SSA002851
        if sequence is not None:
            return sequence
    except Exception:
        pass
    return "ATGCATGCATGC"  # Mock sequence



def test_similarities(test_haplotype_sequence, test_haplotype_ships_df):
    """Generate similarities data using check_similar_match."""
    try:
        import tempfile
        from src.utils.classification_utils import calculate_similarities
        from src.utils.seq_utils import write_multi_fasta
        import pandas as pd

        # Create a combined DataFrame that includes both the query sequence and existing ships
        # Add the query sequence as a new row
        query_row = pd.DataFrame(
            {
                "accession_display": ["query_sequence"],
                "sequence": [test_haplotype_sequence],
            }
        )

        # Combine with existing ships
        if "accession_display" in test_haplotype_ships_df.columns:
            combined_df = pd.concat(
                [query_row, test_haplotype_ships_df], ignore_index=True
            )
        else:
            # If no accession_display, create it from accession_tag
            ships_with_display = test_haplotype_ships_df.copy()
            ships_with_display["accession_display"] = ships_with_display[
                "accession_tag"
            ]
            combined_df = pd.concat([query_row, ships_with_display], ignore_index=True)

        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
        write_multi_fasta(
            combined_df,
            tmp_fasta,
            sequence_col="sequence",
            id_col="accession_display",
        )

        # calculate_similarities returns a nested dictionary, but cluster_sequences expects a list of tuples
        similarities_dict = calculate_similarities(
            fasta_file=tmp_fasta,
            seq_type="nucl",
        )

        # Convert nested dictionary to list of tuples format
        similarities_list = []
        for seq_id1, inner_dict in similarities_dict.items():
            for seq_id2, similarity in inner_dict.items():
                if seq_id1 != seq_id2:  # Skip self-comparisons
                    similarities_list.append((seq_id1, seq_id2, similarity))

        return similarities_list
    except Exception:
        pass

    # Fallback to mock similarities data, containing 3 sequences
    return [
        ("query_sequence", "SSA002851.1", 0.99),
        ("query_sequence", "SSA002904.1", 0.85),
        ("query_sequence", "SSA002596.1", 0.85),
        ("SSA002851.1", "SSA002904.1", 0.98),
    ]
