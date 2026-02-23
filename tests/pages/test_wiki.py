import pytest
import numpy as np
import pandas as pd
from unittest.mock import patch

# Mock dash.register_page before importing the wiki module
with patch("dash.register_page"):
    from src.pages.wiki import (
        create_accordion_item,
        load_initial_data,
        generate_download_helper,
        update_table_stats,
        update_download_selected_button,
        populate_search_components,
    )


class TestWikiCoreFunctions:
    """Test core Wiki page functions that don't require Dash context"""

    def test_create_accordion_item_basic(self):
        """Test basic accordion item creation"""
        meta_df = pd.DataFrame(
            {
                "familyName": ["Family1", "Family1"],
                "accession_tag": ["SSA000001", "SSA000002"],
                "elementLength": [1000, 2000],
                "upDR": ["ATGC", "ATGC"],
                "downDR": ["GCAT", "GCAT"],
            }
        )
        papers_df = pd.DataFrame(
            {
                "familyName": ["Family1"],
                "type_element_reference": ["Ref1"],
                "Url": ["http://example.com"],
            }
        )

        result = create_accordion_item(meta_df, papers_df, "Family1")

        assert result is not None
        assert hasattr(result, "title")
        assert result.title == "Family1"

    def test_create_accordion_item_nan_category(self):
        """Test that nan category returns None"""
        meta_df = pd.DataFrame({"familyName": ["Family1"]})
        papers_df = pd.DataFrame({"familyName": ["Family1"]})

        result = create_accordion_item(meta_df, papers_df, "nan")
        assert result is None

    @patch("src.pages.wiki.cache")
    @patch("src.pages.wiki.fetch_meta_data")
    def test_load_initial_data_success(self, mock_fetch, mock_cache):
        """Test successful data loading"""
        mock_df = pd.DataFrame({"test": [1, 2, 3]})
        mock_cache.get.return_value = mock_df

        result = load_initial_data()

        assert result == mock_df.to_dict("records")
        mock_fetch.assert_not_called()

    @patch("src.pages.wiki.cache")
    @patch("src.pages.wiki.fetch_meta_data")
    def test_load_initial_data_cache_miss(self, mock_fetch, mock_cache):
        """Test data loading when cache is empty"""
        mock_cache.get.return_value = None
        mock_df = pd.DataFrame({"test": [1, 2, 3]})
        mock_fetch.return_value = mock_df

        result = load_initial_data()

        assert result == mock_df.to_dict("records")
        mock_fetch.assert_called_once()
        mock_cache.set.assert_called_once()

    @patch("src.pages.wiki.fetch_ships")
    def test_generate_download_helper_success(self, mock_fetch_ships):
        """Test successful download generation"""
        rows = [
            {"accession_tag": "SSA000001", "familyName": "Family1"},
            {"accession_tag": "SSA000002", "familyName": "Family2"},
        ]

        mock_df = pd.DataFrame(
            {
                "accession_tag": ["SSA000001", "SSA000002"],
                "familyName": ["Family1", "Family2"],
                "sequence": ["ATGCATGC", "TTGGTTGG"],
            }
        )
        mock_fetch_ships.return_value = mock_df

        result, count = generate_download_helper(rows, True, False)

        assert result is not None
        assert "content" in result
        assert "filename" in result
        assert "type" in result
        assert count == 2
        assert "ATGCATGC" in result["content"]
        assert "TTGGTTGG" in result["content"]

    @patch("src.pages.wiki.fetch_ships")
    def test_generate_download_helper_no_rows(self, mock_fetch_ships):
        """Test download helper with no rows"""
        result, count = generate_download_helper([], True, False)

        assert result is None
        assert count == 0
        mock_fetch_ships.assert_not_called()

    def test_update_table_stats_with_data(self):
        """Test table stats with valid data"""
        meta_data = [
            {"accession_tag": "SSA000001", "curated_status": "curated"},
            {"accession_tag": "SSA000002", "curated_status": "uncurated"},
            {"accession_tag": "SSA000003", "curated_status": "curated"},
        ]

        result = update_table_stats(meta_data, None, False, False)
        assert result == "Showing 3 records"

    def test_update_table_stats_with_curated_filter(self):
        """Test table stats with curated filter"""
        meta_data = [
            {"accession_tag": "SSA000001", "curated_status": "curated"},
            {"accession_tag": "SSA000002", "curated_status": "uncurated"},
            {"accession_tag": "SSA000003", "curated_status": "curated"},
        ]

        result = update_table_stats(meta_data, None, True, False)
        assert result == "Showing 2 records"

    def test_update_table_stats_no_data(self):
        """Test table stats with no data"""
        result = update_table_stats(None, None, False, False)
        assert result == "No data available"

    def test_update_download_selected_button_with_selection(self):
        """Test button state with selected rows"""
        selected_rows = [
            {"accession_tag": "SSA000001"},
            {"accession_tag": "SSA000002"},
        ]

        result = update_download_selected_button(selected_rows)
        assert result is False  # Button should be enabled

    def test_update_download_selected_button_no_selection(self):
        """Test button state with no selected rows"""
        selected_rows = []

        result = update_download_selected_button(selected_rows)
        assert result is True  # Button should be disabled

    def test_populate_search_components_with_data(self):
        """Test populating search components with data"""
        meta_data = [
            {
                "name": "Species1",
                "family": "Family1",
                "genus": "Genus1",
                "familyName": "StarshipFamily1",
            },
            {
                "name": "Species2",
                "family": "Family2",
                "genus": "Genus2",
                "familyName": "StarshipFamily2",
            },
        ]

        taxa_data, family_data = populate_search_components(meta_data)

        assert len(taxa_data) > 0
        assert len(family_data) > 0
        assert any(item["value"] == "Species1" for item in taxa_data)
        assert any(item["value"] == "StarshipFamily1" for item in family_data)

    def test_populate_search_components_no_data(self):
        """Test populating search components with no data"""
        taxa_data, family_data = populate_search_components(None)

        assert taxa_data == []
        assert family_data == []

    def test_populate_search_components_empty_data(self):
        """Test populating search components with empty data"""
        taxa_data, family_data = populate_search_components([])

        assert taxa_data == []
        assert family_data == []


class TestWikiDataProcessing:
    """Test data processing functions"""

    def test_accordion_creation_with_complete_data(self):
        """Test accordion creation with complete data"""
        meta_data = pd.DataFrame(
            {
                "accession_tag": ["SSA000001", "SSA000002", "SSA000003"],
                "familyName": ["Family1", "Family1", "Family2"],
                "curated_status": ["curated", "uncurated", "curated"],
                "elementLength": [1000, 2000, 1500],
                "upDR": ["ATGC", "ATGC", "TTAA"],
                "downDR": ["GCAT", "GCAT", "AATT"],
                "name": ["Species1", "Species2", "Species3"],
                "order": ["Order1", "Order1", "Order2"],
                "family": ["Family1", "Family1", "Family2"],
            }
        )

        papers_data = pd.DataFrame(
            {
                "familyName": ["Family1", "Family2"],
                "type_element_reference": ["Ref1", "Ref2"],
                "Url": ["http://example.com/1", "http://example.com/2"],
            }
        )

        # Test creating accordion items for each family
        family1_item = create_accordion_item(meta_data, papers_data, "Family1")
        family2_item = create_accordion_item(meta_data, papers_data, "Family2")

        assert family1_item is not None
        assert family2_item is not None
        assert family1_item.title == "Family1"
        assert family2_item.title == "Family2"

    def test_table_stats_integration(self):
        """Test table stats with complete data"""
        meta_data = [
            {"accession_tag": "SSA000001", "curated_status": "curated"},
            {"accession_tag": "SSA000002", "curated_status": "uncurated"},
            {"accession_tag": "SSA000003", "curated_status": "curated"},
        ]

        # Test without filters
        result = update_table_stats(meta_data, None, False, False)
        assert result == "Showing 3 records"

        # Test with curated filter
        result = update_table_stats(meta_data, None, True, False)
        assert result == "Showing 2 records"

        # Test with dereplicated filter
        result = update_table_stats(meta_data, None, False, True)
        assert result == "Showing 3 records"

    def test_download_helper_with_duplicate_accessions(self):
        """Test download helper with duplicate accession tags"""
        rows = [
            {"accession_tag": "SSA000001", "familyName": "Family1"},
            {"accession_tag": "SSA000001", "familyName": "Family1"},  # Duplicate
        ]

        with patch("src.pages.wiki.fetch_ships") as mock_fetch:
            mock_df = pd.DataFrame(
                {
                    "accession_tag": ["SSA000001", "SSA000001"],
                    "familyName": ["Family1", "Family1"],
                    "sequence": ["ATGCATGC", "ATGCATGC"],
                }
            )
            mock_fetch.return_value = mock_df

            result, count = generate_download_helper(rows, True, False)

            assert result is not None
            # The function deduplicates by accession_tag and sequence, so count should be 1
            assert count == 1


class TestWikiEdgeCases:
    """Test edge cases and error handling"""

    def test_accordion_item_with_missing_columns(self):
        """Test accordion item creation with missing columns - should handle gracefully"""
        df = pd.DataFrame(
            {
                "familyName": ["Family1"],
                "accession_tag": ["SSA000001"],
            }
        )
        papers = pd.DataFrame(
            {
                "familyName": ["Family1"],
                "type_element_reference": ["Ref1"],  # Add required column
                "Url": ["http://example.com"],
            }
        )

        # Should handle missing elementLength column gracefully
        with pytest.raises(KeyError):
            create_accordion_item(df, papers, "Family1")

    def test_accordion_item_with_nan_values(self):
        """Test accordion item creation with NaN values"""
        df = pd.DataFrame(
            {
                "familyName": ["Family1", "Family1"],
                "accession_tag": ["SSA000001", None],
                "elementLength": [1000, np.nan],
                "upDR": ["ATGC", None],
                "downDR": ["GCAT", None],
            }
        )
        papers = pd.DataFrame(
            {
                "familyName": ["Family1"],
                "type_element_reference": ["Ref1"],
                "Url": ["http://example.com"],
            }
        )

        result = create_accordion_item(df, papers, "Family1")
        assert result is not None

    def test_generate_download_helper_exception(self):
        """Test download helper when exception occurs"""
        rows = [{"accession_tag": "SSA000001", "familyName": "Family1"}]

        with patch("src.pages.wiki.fetch_ships") as mock_fetch:
            mock_fetch.side_effect = Exception("Database error")

            result, count = generate_download_helper(rows, True, False)

            assert result is None
            assert count == 0

    def test_update_table_stats_exception(self):
        """Test table stats when exception occurs"""
        # Pass invalid data to trigger exception
        result = update_table_stats("invalid_data", None, False, False)
        assert result == "Error loading data"

    def test_load_initial_data_exception(self):
        """Test loading data when exception occurs"""
        with patch("src.pages.wiki.cache") as mock_cache:
            mock_cache.get.side_effect = Exception("Cache error")

            result = load_initial_data()

            assert result is None


class TestWikiSearchComponents:
    """Test search component functionality"""

    def test_search_components_with_various_taxonomy_levels(self):
        """Test search components with various taxonomy levels"""
        meta_data = [
            {
                "name": "Species1",
                "subkingdom": "Subkingdom1",
                "phylum": "Phylum1",
                "class": "Class1",
                "order": "Order1",
                "family": "Family1",
                "genus": "Genus1",
                "familyName": "StarshipFamily1",
            },
            {
                "name": "Species2",
                "subkingdom": "Subkingdom2",
                "phylum": "Phylum2",
                "class": "Class2",
                "order": "Order2",
                "family": "Family2",
                "genus": "Genus2",
                "familyName": "StarshipFamily2",
            },
        ]

        taxa_data, family_data = populate_search_components(meta_data)

        # Check that all taxonomy levels are included
        all_values = [item["value"] for item in taxa_data]
        assert "Species1" in all_values
        assert "Subkingdom1" in all_values
        assert "Phylum1" in all_values
        assert "Class1" in all_values
        assert "Order1" in all_values
        assert "Family1" in all_values
        assert "Genus1" in all_values

        # Check family data
        family_values = [item["value"] for item in family_data]
        assert "StarshipFamily1" in family_values
        assert "StarshipFamily2" in family_values

    def test_search_components_with_missing_columns(self):
        """Test search components with missing columns"""
        meta_data = [
            {
                "name": "Species1",
                "familyName": "StarshipFamily1",
            },  # Missing some columns
            {"family": "Family2", "genus": "Genus2"},  # Missing name and familyName
        ]

        taxa_data, family_data = populate_search_components(meta_data)

        # Should still work with missing columns
        assert len(taxa_data) > 0
        assert len(family_data) > 0
        assert any(item["value"] == "Species1" for item in taxa_data)
        assert any(item["value"] == "StarshipFamily1" for item in family_data)

    def test_search_components_data_format(self):
        """Test that search components return correct format"""
        meta_data = [
            {"name": "Species1", "familyName": "StarshipFamily1"},
        ]

        taxa_data, family_data = populate_search_components(meta_data)

        # Check format
        assert isinstance(taxa_data, list)
        assert isinstance(family_data, list)

        if taxa_data:
            assert "value" in taxa_data[0]
            assert "label" in taxa_data[0]

        if family_data:
            assert "value" in family_data[0]
            assert "label" in family_data[0]
