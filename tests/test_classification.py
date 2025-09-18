from src.utils.seq_utils import write_temp_fasta
from src.utils.classification_utils import (
    assign_accession,
    check_exact_match,
    check_contained_match,
    check_similar_match,
    generate_new_accession,
)

def test_exact_match(test_ships_df, test_sequence, test_sequence_revcomp):
    """Test matching of exact sequences using md5sums"""

    test_sequence_fasta = write_temp_fasta(header="sequence", sequence=test_sequence)
    test_sequence_revcomp_fasta = write_temp_fasta(header="revcomp", sequence=test_sequence_revcomp)

    result = check_exact_match(fasta=test_sequence_fasta, existing_ships=test_ships_df)
    result_revcomp = check_exact_match(fasta=test_sequence_revcomp_fasta, existing_ships=test_ships_df)
    assert result == "SBS000002.1"  # Updated to match actual database data with version
    assert result_revcomp == "SBS000002.1"  # Updated to match actual database data with version

def test_contained_match(test_ships_df, test_contained_sequence):
    """Test contained sequence matching."""
    # Test contained match - should return SBS000002 as it's longer
    result = check_contained_match(test_contained_sequence, test_ships_df)
    assert result == "SBS000002.1"  # Now matches the longer sequence

def test_similar_match(
    test_ships_df, test_similar_sequence
):
    """Test similar sequence matching using k-mer comparison."""
    # Test similar match
    result, similarities = check_similar_match(test_similar_sequence, test_ships_df, threshold=0.8)
    assert result == "SBS000002.1"


def test_generate_new_accession(test_ships_df):
    """Test new accession number generation."""
    result = generate_new_accession(test_ships_df)
    # check that the result is an accession that doesn't yet exist in ships_df
    # create a list to check against
    existing_accessions = [ship.accession_tag for ship in test_ships_df.itertuples()]
    assert result not in existing_accessions

def test_assign_accession_exact_match(test_ships_df, test_sequence):
    """Test full workflow with exact match."""
    accession, needs_review = assign_accession(
        test_sequence, existing_ships=test_ships_df
    )
    assert accession == "SBS000002.1"
    assert needs_review is False


def test_assign_accession_contained_match(test_ships_df, test_contained_sequence):
    """Test full workflow with contained match."""
    accession, needs_review = assign_accession(
        test_contained_sequence, existing_ships=test_ships_df
    )
    assert accession == "SBS000002.1"
    assert needs_review is True

def test_assign_accession_similar_match(
    test_ships_df, test_similar_sequence
):
    """Test full workflow with similar match."""
    accession, needs_review = assign_accession(
        test_similar_sequence, existing_ships=test_ships_df, threshold=0.8
    )
    assert accession == "SBS000002.1"
    assert needs_review is True


def test_assign_accession_new(test_ships_df):
    """Test full workflow with new sequence."""
    accession, needs_review = assign_accession("TTTTTTTT", existing_ships=test_ships_df)
    assert accession == "SBS000003.1"
    assert needs_review is True

# TODO: add test for adding a new accession version