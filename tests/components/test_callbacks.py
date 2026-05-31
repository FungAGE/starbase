import pytest

pytestmark = pytest.mark.skip(
    reason="create_modal was removed; accession modals are rendered client-side"
)


def test_create_modal_ship_accession(ship_accession_modal_data):
    assert ship_accession_modal_data is not None


def test_create_modal_accession(accession_modal_data):
    assert accession_modal_data is not None


def test_create_modal_error_handling():
    pass
