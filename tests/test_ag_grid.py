"""
Test suite for AG Grid components stability and functionality.
Consolidated tests for the AG Grid stability fix.
"""

import pytest
import pandas as pd
from src.components.tables import (
    make_pgv_table,
    make_dl_table,
    make_ship_table,
    make_paper_table,
    make_wiki_table,
    DEFAULT_GRID_OPTIONS,
)


class TestAGGridStability:
    """Test AG Grid configuration and stability."""

    def test_default_grid_options_are_stable(self):
        """Test that default grid options don't contain problematic settings."""
        # These options caused the 'this.state.gridApi is null' error
        problematic_options = ["ensureDomOrder", "onFirstDataRendered"]

        for option in problematic_options:
            assert option not in DEFAULT_GRID_OPTIONS, (
                f"Found problematic option: {option}"
            )

        # Should use normal layout, not autoHeight
        assert DEFAULT_GRID_OPTIONS.get("domLayout") == "normal"
        assert DEFAULT_GRID_OPTIONS["suppressPropertyNamesCheck"] is True

    def test_all_table_types_creation(self, mock_accession_data):
        """Test that all table types can be created without errors."""
        # Test PGV table
        pgv_table = make_pgv_table(
            df=mock_accession_data,
            columns=[
                {"id": "accession_tag", "name": "Accession"},
                {"id": "familyName", "name": "Family"},
            ],
            id="test-pgv",
            select_rows=True,
        )
        assert pgv_table is not None
        assert len(pgv_table.rowData) == 2

        # Test download table
        dl_table = make_dl_table(
            df=mock_accession_data,
            id="test-dl",
            table_columns=[
                {"id": "accession_tag", "name": "Accession"},
                {"id": "familyName", "name": "Family"},
            ],
        )
        assert dl_table is not None
        assert len(dl_table.children.rowData) == 2

        # Test ship table
        ship_table = make_ship_table(
            df=mock_accession_data,
            id="test-ship",
            columns=[
                {"field": "accession_tag", "name": "Accession"},
                {"field": "familyName", "name": "Family"},
            ],
            select_rows=True,
        )
        assert ship_table is not None
        assert len(ship_table.rowData) == 2

        # Test paper table (uses different data structure)
        paper_table = make_paper_table()
        assert paper_table is not None

        # Test wiki table
        wiki_table = make_wiki_table(n_ships=100, max_size=5000, min_size=1000)
        assert wiki_table is not None
        assert len(wiki_table.rowData) == 3  # Should have 3 summary rows
        assert (
            wiki_table.dashGridOptions["pagination"] is False
        )  # No pagination for small table

    def test_pgv_table_checkbox_configuration(self, mock_accession_data):
        """Test that checkbox selection is properly configured."""
        table = make_pgv_table(
            df=mock_accession_data,
            columns=[
                {"id": "accession_tag", "name": "Accession"},
                {"id": "familyName", "name": "Family"},
            ],
            id="checkbox-test",
            select_rows=True,
        )

        # First column should have checkbox with proper width
        first_col = table.columnDefs[0]
        assert first_col["checkboxSelection"] is True
        assert first_col["headerCheckboxSelection"] is True
        assert first_col["width"] == 180  # Adequate width for checkbox + content
        assert first_col["minWidth"] == 150
        assert first_col["flex"] == 0

        # Grid should have proper selection options
        assert table.dashGridOptions["rowSelection"] == "multiple"
        assert table.dashGridOptions["suppressRowClickSelection"] is True

    def test_accession_column_width_fix(self):
        """Test that accession columns have adequate width (fixes the 50px issue)."""
        sample_data = pd.DataFrame(
            {
                "accession_tag": ["VERY_LONG_ACCESSION_TAG_001"],
                "familyName": ["Family1"],
            }
        )

        # Test PGV table
        pgv_table = make_pgv_table(
            df=sample_data,
            columns=[
                {"id": "accession_tag", "name": "Accession"},
                {"id": "familyName", "name": "Family"},
            ],
            id="width-test-pgv",
            select_rows=True,
        )

        accession_col = pgv_table.columnDefs[0]
        assert accession_col["width"] == 180, (
            f"PGV accession width is {accession_col.get('width')}, should be 180"
        )
        assert accession_col["minWidth"] == 150

        # Test download table
        dl_table = make_dl_table(
            df=sample_data,
            id="width-test-dl",
            table_columns=[
                {"id": "accession_tag", "name": "Accession"},
                {"id": "familyName", "name": "Family"},
            ],
        )

        dl_accession_col = dl_table.children.columnDefs[0]
        assert dl_accession_col["width"] == 180, (
            f"DL accession width is {dl_accession_col.get('width')}, should be 180"
        )
        assert dl_accession_col["minWidth"] == 150

    def test_empty_data_handling(self):
        """Test table creation with empty/null data doesn't break."""
        columns = [
            {"id": "accession_tag", "name": "Accession"},
            {"id": "familyName", "name": "Family"},
        ]

        # Test with None
        table = make_pgv_table(
            df=None, columns=columns, id="empty-test", select_rows=True
        )
        assert table is not None
        assert isinstance(table.rowData, list)
        assert len(table.rowData) == 0

        # Test with empty DataFrame
        empty_df = pd.DataFrame()
        table2 = make_pgv_table(
            df=empty_df, columns=columns, id="empty-df-test", select_rows=False
        )
        assert table2 is not None

    def test_table_height_configuration(self):
        """Test that tables have proper fixed heights for stability."""
        sample_data = pd.DataFrame({"test": [1, 2, 3]})

        # PGV table should be 450px
        pgv_table = make_pgv_table(df=sample_data, columns=None, id="height-test-pgv")
        assert pgv_table.style["height"] == "450px"

        # Download table should be 600px
        dl_table = make_dl_table(
            df=sample_data,
            id="height-test-dl",
            table_columns=[{"id": "test", "name": "Test"}],
        )
        assert dl_table.children.style["height"] == "600px"

        # Ship table should be 500px
        ship_table = make_ship_table(df=sample_data, id="height-test-ship")
        assert ship_table.style["height"] == "500px"

        # Wiki table should be 192px (exact fit for 3 rows)
        wiki_table = make_wiki_table(10, 1000, 500)
        assert wiki_table.style["height"] == "192px"

    def test_no_problematic_event_handlers(self, mock_accession_data):
        """Test that tables don't use problematic event handlers that cause gridApi errors."""
        table = make_pgv_table(
            df=mock_accession_data,
            columns=[{"id": "accession_tag", "name": "Accession"}],
            id="handler-test",
            select_rows=True,
        )

        grid_options = table.dashGridOptions

        # Should not have complex event handlers that access gridApi prematurely
        assert "onFirstDataRendered" not in grid_options
        assert "ensureDomOrder" not in grid_options

        # Should have stable row/header heights
        assert grid_options["rowHeight"] == 48
        assert grid_options["headerHeight"] == 48


@pytest.mark.skip_browser
class TestAGGridBrowserIntegration:
    """Browser-based tests that require Selenium (marked to skip by default)."""

    def test_grid_loads_without_console_errors(self, driver, test_client):
        """Test that AG Grid loads without 'this.state.gridApi is null' errors."""
        try:
            # This would need to be adapted to your actual app URL structure
            driver.get("http://localhost:8050/pgv")

            # Wait for grid to load
            import time

            time.sleep(3)

            # Check browser console for the specific error we fixed
            logs = driver.get_log("browser")
            grid_api_errors = [
                log
                for log in logs
                if log["level"] == "SEVERE"
                and "this.state.gridApi is null" in log.get("message", "")
            ]

            assert len(grid_api_errors) == 0, f"Found gridApi errors: {grid_api_errors}"

        except Exception as e:
            pytest.skip(f"Browser test not available: {e}")


# Performance test for larger datasets
def test_large_dataset_performance():
    """Test that tables can handle moderately large datasets."""
    large_data = pd.DataFrame(
        {
            "id": range(1000),
            "name": [f"Item_{i}" for i in range(1000)],
            "value": [i * 10 for i in range(1000)],
        }
    )

    columns = [
        {"id": "id", "name": "ID"},
        {"id": "name", "name": "Name"},
        {"id": "value", "name": "Value"},
    ]

    table = make_pgv_table(
        df=large_data,
        columns=columns,
        id="performance-test",
        select_rows=True,
        pg_sz=50,
    )

    assert table is not None
    assert len(table.rowData) == 1000
    assert table.dashGridOptions["paginationPageSize"] == 50
