import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
import dash
from dash import html

# Mock dash.register_page before importing the wiki module
with patch('dash.register_page'):
    from src.pages.wiki import (
        load_meta_data,
        load_paper_data,
        create_accordion,
        create_search_results,
        populate_search_components,
        handle_taxa_and_family_search,
        update_search_sunburst,
        generate_download_all,
        generate_download_selected,
    )


class TestLoadMetaDataCallback:
    """Test the load_meta_data callback function"""
    
    @patch('src.pages.wiki.cache')
    @patch('src.pages.wiki.fetch_meta_data')
    def test_load_meta_data_cache_hit(self, mock_fetch, mock_cache):
        """Test loading meta data when cache has data"""
        mock_df = pd.DataFrame({"test": [1, 2, 3]})
        mock_cache.get.return_value = mock_df
        
        result = load_meta_data("http://example.com")
        
        assert result == mock_df.to_dict("records")
        mock_fetch.assert_not_called()
    
    @patch('src.pages.wiki.cache')
    @patch('src.pages.wiki.fetch_meta_data')
    def test_load_meta_data_cache_miss(self, mock_fetch, mock_cache):
        """Test loading meta data when cache is empty"""
        mock_cache.get.return_value = None
        mock_df = pd.DataFrame({"test": [1, 2, 3]})
        mock_fetch.return_value = mock_df
        
        result = load_meta_data("http://example.com")
        
        assert result == mock_df.to_dict("records")
        mock_fetch.assert_called_once()
        mock_cache.set.assert_called_once()
    
    def test_load_meta_data_no_url(self):
        """Test loading meta data with no URL"""
        with pytest.raises(dash.exceptions.PreventUpdate):
            load_meta_data(None)
    
    @patch('src.pages.wiki.cache')
    @patch('src.pages.wiki.fetch_meta_data')
    def test_load_meta_data_fetch_failure(self, mock_fetch, mock_cache):
        """Test loading meta data when fetch fails"""
        mock_cache.get.return_value = None
        mock_fetch.return_value = None
        
        result = load_meta_data("http://example.com")
        
        assert result == []
    
    @patch('src.pages.wiki.cache')
    @patch('src.pages.wiki.fetch_meta_data')
    def test_load_meta_data_with_contig_cleaning(self, mock_fetch, mock_cache):
        """Test loading meta data with contigID cleaning"""
        mock_cache.get.return_value = None
        mock_df = pd.DataFrame({
            "test": [1, 2, 3],
            "contigID": ["contig_1", "contig_2", "contig_3"]
        })
        mock_fetch.return_value = mock_df
        
        with patch('src.pages.wiki.clean_contigIDs') as mock_clean:
            mock_clean.side_effect = lambda x: x.replace("contig_", "clean_")
            
            result = load_meta_data("http://example.com")
            
            assert result == mock_df.to_dict("records")
            mock_clean.assert_called()


class TestLoadPaperDataCallback:
    """Test the load_paper_data callback function"""
    
    @patch('src.pages.wiki.cache')
    @patch('src.pages.wiki.fetch_paper_data')
    def test_load_paper_data_cache_hit(self, mock_fetch, mock_cache):
        """Test loading paper data when cache has data"""
        mock_df = pd.DataFrame({"test": [1, 2, 3]})
        mock_cache.get.return_value = mock_df
        
        result = load_paper_data("http://example.com")
        
        assert result == mock_df.to_dict("records")
        mock_fetch.assert_not_called()
    
    @patch('src.pages.wiki.cache')
    @patch('src.pages.wiki.fetch_paper_data')
    def test_load_paper_data_cache_miss(self, mock_fetch, mock_cache):
        """Test loading paper data when cache is empty"""
        mock_cache.get.return_value = None
        mock_df = pd.DataFrame({"test": [1, 2, 3]})
        mock_fetch.return_value = mock_df
        
        result = load_paper_data("http://example.com")
        
        assert result == mock_df.to_dict("records")
        mock_fetch.assert_called_once()
    
    def test_load_paper_data_no_url(self):
        """Test loading paper data with no URL"""
        with pytest.raises(dash.exceptions.PreventUpdate):
            load_paper_data(None)
    
    @patch('src.pages.wiki.fetch_paper_data')
    def test_load_paper_data_sqlalchemy_error(self, mock_fetch):
        """Test loading paper data with SQLAlchemy error"""
        from sqlalchemy.exc import SQLAlchemyError
        mock_fetch.side_effect = SQLAlchemyError("Database error")
        
        result = load_paper_data("http://example.com")
        
        assert result == []


class TestCreateAccordionCallback:
    """Test the create_accordion callback function"""
    
    def test_create_accordion_with_data(self):
        """Test accordion creation with valid data"""
        meta_data = [
            {"familyName": "Family1", "accession_tag": "SSA000001"},
            {"familyName": "Family2", "accession_tag": "SSA000002"},
        ]
        paper_data = [
            {"familyName": "Family1", "type_element_reference": "Ref1"},
            {"familyName": "Family2", "type_element_reference": "Ref2"},
        ]
        
        with patch('src.pages.wiki.create_accordion_item') as mock_create:
            mock_create.return_value = MagicMock()
            
            result = create_accordion(meta_data, paper_data)
            
            assert result is not None
            mock_create.assert_called()
    
    def test_create_accordion_none_inputs(self):
        """Test accordion creation with None inputs"""
        with pytest.raises(dash.exceptions.PreventUpdate):
            create_accordion(None, None)
    
    def test_create_accordion_one_none_input(self):
        """Test accordion creation with one None input"""
        meta_data = [{"familyName": "Family1"}]
        
        with pytest.raises(dash.exceptions.PreventUpdate):
            create_accordion(meta_data, None)


class TestCreateSearchResultsCallback:
    """Test the create_search_results callback function"""
    
    def test_create_search_results_with_filtered_data(self):
        """Test search results with filtered data"""
        filtered_data = [
            {"accession_tag": "SSA000001", "familyName": "Family1", "curated_status": "curated"},
            {"accession_tag": "SSA000002", "familyName": "Family1", "curated_status": "uncurated"},
        ]
        cached_data = [
            {"accession_tag": "SSA000003", "familyName": "Family2", "curated_status": "curated"},
        ]
        
        with patch('src.pages.wiki.make_dl_table') as mock_table:
            mock_table.return_value = html.Div("Mock Table")
            
            result = create_search_results(filtered_data, cached_data, False, False)
            
            assert result is not None
            mock_table.assert_called_once()
    
    def test_create_search_results_no_data(self):
        """Test search results with no data"""
        result = create_search_results(None, None, False, False)
        
        assert "Start a search to see results" in str(result)
    
    def test_create_search_results_empty_data(self):
        """Test search results with empty data"""
        result = create_search_results([], None, False, False)
        
        assert "No results match your search criteria" in str(result)
    
    def test_create_search_results_with_curated_filter(self):
        """Test search results with curated filter"""
        data = [
            {"accession_tag": "SSA000001", "familyName": "Family1", "curated_status": "curated"},
            {"accession_tag": "SSA000002", "familyName": "Family1", "curated_status": "uncurated"},
        ]
        
        with patch('src.pages.wiki.make_dl_table') as mock_table:
            mock_table.return_value = html.Div("Mock Table")
            
            result = create_search_results(data, None, True, False)
            
            assert result is not None
            # Should filter to only curated items
            call_args = mock_table.call_args[0]
            table_data = call_args[0]
            assert len(table_data) == 1  # Only curated item should remain


class TestPopulateSearchComponentsCallback:
    """Test the populate_search_components callback function"""
    
    def test_populate_search_components_with_data(self):
        """Test populating search components with data"""
        meta_data = [
            {"name": "Species1", "family": "Family1", "genus": "Genus1", "familyName": "StarshipFamily1"},
            {"name": "Species2", "family": "Family2", "genus": "Genus2", "familyName": "StarshipFamily2"},
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


class TestHandleTaxaAndFamilySearchCallback:
    """Test the handle_taxa_and_family_search callback function"""
    
    def test_handle_search_reset(self):
        """Test search reset functionality"""
        original_data = [{"name": "Species1", "familyName": "Family1"}]
        
        result = handle_taxa_and_family_search(
            search_clicks=0,
            reset_clicks=1,
            taxa_search_value="",
            family_search_value="",
            original_data=original_data
        )
        
        assert result is None
    
    def test_handle_search_taxa_filter(self):
        """Test taxa search filtering"""
        original_data = [
            {"name": "Species1", "familyName": "Family1"},
            {"name": "Species2", "familyName": "Family2"},
        ]
        
        result = handle_taxa_and_family_search(
            search_clicks=1,
            reset_clicks=0,
            taxa_search_value="Species1",
            family_search_value="",
            original_data=original_data
        )
        
        assert len(result) == 1
        assert result[0]["name"] == "Species1"
    
    def test_handle_search_family_filter(self):
        """Test family search filtering"""
        original_data = [
            {"name": "Species1", "familyName": "Family1"},
            {"name": "Species2", "familyName": "Family2"},
        ]
        
        result = handle_taxa_and_family_search(
            search_clicks=1,
            reset_clicks=0,
            taxa_search_value="",
            family_search_value="Family1",
            original_data=original_data
        )
        
        assert len(result) == 1
        assert result[0]["familyName"] == "Family1"
    
    def test_handle_search_no_original_data(self):
        """Test search with no original data"""
        with pytest.raises(dash.exceptions.PreventUpdate):
            handle_taxa_and_family_search(
                search_clicks=1,
                reset_clicks=0,
                taxa_search_value="test",
                family_search_value="",
                original_data=None
            )
    
    def test_handle_search_no_clicks(self):
        """Test search with no button clicks"""
        original_data = [{"name": "Species1"}]
        
        result = handle_taxa_and_family_search(
            search_clicks=0,
            reset_clicks=0,
            taxa_search_value="test",
            family_search_value="",
            original_data=original_data
        )
        
        assert result is None


class TestUpdateSearchSunburstCallback:
    """Test the update_search_sunburst callback function"""
    
    def test_update_sunburst_with_meta_data(self):
        """Test sunburst update with meta data"""
        meta_data = [
            {"name": "Species1", "family": "Family1", "curated_status": "curated"},
            {"name": "Species2", "family": "Family2", "curated_status": "uncurated"},
        ]
        
        with patch('src.pages.wiki.create_sunburst_plot') as mock_plot:
            mock_figure = MagicMock()
            mock_plot.return_value = mock_figure
            
            result = update_search_sunburst(None, meta_data, False, False)
            
            assert result is not None
            mock_plot.assert_called_once()
    
    def test_update_sunburst_with_filtered_data(self):
        """Test sunburst update with filtered data"""
        meta_data = [
            {"name": "Species1", "family": "Family1", "curated_status": "curated"},
            {"name": "Species2", "family": "Family2", "curated_status": "uncurated"},
        ]
        filtered_data = [
            {"name": "Species1", "family": "Family1", "curated_status": "curated"},
        ]
        
        with patch('src.pages.wiki.create_sunburst_plot') as mock_plot:
            mock_figure = MagicMock()
            mock_plot.return_value = mock_figure
            
            result = update_search_sunburst(filtered_data, meta_data, False, False)
            
            assert result is not None
            mock_plot.assert_called_once()
    
    def test_update_sunburst_no_data(self):
        """Test sunburst update with no data"""
        result = update_search_sunburst(None, None, False, False)
        
        assert "No data available" in str(result)
    
    def test_update_sunburst_empty_data(self):
        """Test sunburst update with empty data"""
        result = update_search_sunburst([], None, False, False)
        
        assert "No results to display" in str(result)


class TestGenerateDownloadCallbacks:
    """Test the download generation callback functions"""
    
    def test_generate_download_all_success(self):
        """Test successful download all generation"""
        table_data = [
            {"accession_tag": "SSA000001", "familyName": "Family1"},
            {"accession_tag": "SSA000002", "familyName": "Family2"},
        ]
        
        with patch('src.pages.wiki.generate_download_helper') as mock_helper:
            mock_helper.return_value = (
                {"content": ">SSA000001\nATGC", "filename": "test.fasta", "type": "text/plain"},
                2
            )
            
            download_data, notification = generate_download_all(
                dl_all_clicks=1,
                table_data=table_data,
                curated=False,
                dereplicate=False
            )
            
            assert download_data is not None
            assert "content" in download_data
            assert "Success" in str(notification)
    
    def test_generate_download_all_no_clicks(self):
        """Test download all with no clicks"""
        with pytest.raises(dash.exceptions.PreventUpdate):
            generate_download_all(
                dl_all_clicks=0,
                table_data=[],
                curated=False,
                dereplicate=False
            )
    
    def test_generate_download_selected_success(self):
        """Test successful download selected generation"""
        table_data = [
            {"accession_tag": "SSA000001", "familyName": "Family1"},
            {"accession_tag": "SSA000002", "familyName": "Family2"},
        ]
        selected_rows = [
            {"accession_tag": "SSA000001", "familyName": "Family1"},
        ]
        
        with patch('src.pages.wiki.generate_download_helper') as mock_helper:
            mock_helper.return_value = (
                {"content": ">SSA000001\nATGC", "filename": "test.fasta", "type": "text/plain"},
                1
            )
            
            download_data, notification = generate_download_selected(
                dl_select_clicks=1,
                table_data=table_data,
                selected_rows=selected_rows,
                curated=False,
                dereplicate=False
            )
            
            assert download_data is not None
            assert "content" in download_data
            assert "Success" in str(notification)
    
    def test_generate_download_selected_no_clicks(self):
        """Test download selected with no clicks"""
        with pytest.raises(dash.exceptions.PreventUpdate):
            generate_download_selected(
                dl_select_clicks=0,
                table_data=[],
                selected_rows=[],
                curated=False,
                dereplicate=False
            )
    
    def test_generate_download_selected_no_selection(self):
        """Test download selected with no selected rows"""
        with pytest.raises(dash.exceptions.PreventUpdate):
            generate_download_selected(
                dl_select_clicks=1,
                table_data=[{"accession_tag": "SSA000001"}],
                selected_rows=[],
                curated=False,
                dereplicate=False
            )


# Integration tests for callback workflows
class TestWikiCallbackIntegration:
    """Integration tests for Wiki page callbacks"""
    
    def test_complete_data_loading_workflow(self):
        """Test complete data loading workflow"""
        with patch('src.pages.wiki.cache') as mock_cache, \
             patch('src.pages.wiki.fetch_meta_data') as mock_fetch_meta, \
             patch('src.pages.wiki.fetch_paper_data') as mock_fetch_paper:
            
            # Setup mocks
            mock_cache.get.return_value = None
            mock_meta_df = pd.DataFrame({"test": [1, 2, 3]})
            mock_paper_df = pd.DataFrame({"test": [4, 5, 6]})
            mock_fetch_meta.return_value = mock_meta_df
            mock_fetch_paper.return_value = mock_paper_df
            
            # Test meta data loading
            meta_result = load_meta_data("http://example.com")
            assert meta_result == mock_meta_df.to_dict("records")
            
            # Test paper data loading
            paper_result = load_paper_data("http://example.com")
            assert paper_result == mock_paper_df.to_dict("records")
    
    def test_search_and_filter_workflow(self):
        """Test complete search and filter workflow"""
        original_data = [
            {"name": "Species1", "familyName": "Family1", "curated_status": "curated"},
            {"name": "Species2", "familyName": "Family2", "curated_status": "uncurated"},
        ]
        
        # Test search
        filtered_data = handle_taxa_and_family_search(
            search_clicks=1,
            reset_clicks=0,
            taxa_search_value="Species1",
            family_search_value="",
            original_data=original_data
        )
        
        assert len(filtered_data) == 1
        
        # Test search results creation
        with patch('src.pages.wiki.make_dl_table') as mock_table:
            mock_table.return_value = html.Div("Mock Table")
            
            results = create_search_results(filtered_data, original_data, False, False)
            assert results is not None
    
    def test_download_workflow(self):
        """Test complete download workflow"""
        table_data = [
            {"accession_tag": "SSA000001", "familyName": "Family1"},
        ]
        
        with patch('src.pages.wiki.fetch_ships') as mock_fetch:
            mock_df = pd.DataFrame({
                "accession_tag": ["SSA000001"],
                "familyName": ["Family1"],
                "sequence": ["ATGCATGC"],
            })
            mock_fetch.return_value = mock_df
            
            # Test download helper
            from src.pages.wiki import generate_download_helper
            download_data, count = generate_download_helper(table_data, False, False)
            
            assert download_data is not None
            assert count == 1
            assert "ATGCATGC" in download_data["content"]
