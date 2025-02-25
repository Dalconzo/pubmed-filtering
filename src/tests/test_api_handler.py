import pytest
from unittest.mock import patch, MagicMock
from src.api_handler import PubMedAPI

@pytest.fixture
def api():
    return PubMedAPI(email="test@example.com", api_key="test_key")

@pytest.fixture
def mock_pubmed_response():
    return {
        'IdList': ['12345', '67890'],
        'Count': '2',
        'RetMax': '2'
    }

@pytest.fixture
def mock_article_response():
    return {
        'PubmedArticle': [{
            'MedlineCitation': {
                'Article': {
                    'ArticleTitle': 'Test Title',
                    'Abstract': {'AbstractText': ['Test abstract']},
                    'AuthorList': [
                        {'LastName': 'Smith', 'ForeName': 'John'},
                        {'LastName': 'Doe', 'ForeName': 'Jane'}
                    ],
                    'Journal': {'Title': 'Test Journal'}
                },
                'PMID': '12345'
            }
        }]
    }

@patch('Bio.Entrez.esearch')
@patch('Bio.Entrez.read')
def test_search_pubmed(mock_read, mock_esearch, api, mock_pubmed_response):
    # Set up mock for esearch
    mock_handle = MagicMock()
    mock_esearch.return_value = mock_handle
    
    # Set up mock for Entrez.read
    mock_read.return_value = mock_pubmed_response
    
    results = api.search_pubmed("test query", max_results=2)
    
    # Verify results
    assert len(results) == 2
    assert '12345' in results
    mock_esearch.assert_called_with(db="pubmed", term="test query", retmax=2)

@patch('Bio.Entrez.efetch')
@patch('Bio.Entrez.read')
def test_fetch_article_details(mock_read, mock_efetch, api, mock_article_response):
    # Set up mock for efetch
    mock_handle = MagicMock()
    mock_efetch.return_value = mock_handle
    
    # Set up mock for Entrez.read
    mock_read.return_value = mock_article_response
    
    article = api.fetch_article_details('12345')
    
    # Verify results
    assert article['PubmedArticle'][0]['MedlineCitation']['PMID'] == '12345'
    mock_efetch.assert_called_with(db="pubmed", id='12345', rettype="xml")

@patch('Bio.Entrez.efetch')
@patch('Bio.Entrez.read')
def test_batch_fetch_articles(mock_read, mock_efetch, api, mock_article_response):
    # Set up mock for efetch
    mock_handle = MagicMock()
    mock_efetch.return_value = mock_handle
    
    # Set up mock for Entrez.read to return the expected data
    mock_read.return_value = [mock_article_response]
    
    articles = api.batch_fetch_articles(['12345', '67890'], batch_size=2)
    
    # Verify results
    assert len(articles) == 1
    assert articles[0]['PubmedArticle'][0]['MedlineCitation']['PMID'] == '12345'
    
    # Verify correct API call
    mock_efetch.assert_called_with(db="pubmed", id=['12345', '67890'], rettype="xml")

def test_get_article_metadata(api, mock_article_response):
    metadata = api.get_article_metadata(mock_article_response)
    assert metadata['title'] == 'Test Title'
    assert len(metadata['authors']) == 2
    assert metadata['journal'] == 'Test Journal'
    assert metadata['pmid'] == '12345'
