import pytest
from src.api_handler import PubMedAPI
from src.ai_extractor import GeminiExtractor
from src.config import Config

@pytest.mark.live
class TestLiveAPIs:
    @pytest.fixture
    def config(self):
        return Config()
    
    @pytest.fixture
    def pubmed_api(self, config):
        credentials = config.get_pubmed_credentials()
        return PubMedAPI(email=credentials['email'], api_key=credentials['api_key'])
    
    @pytest.fixture
    def gemini_api(self, config):
        api_key = config.get_gemini_key()
        if not api_key:
            pytest.skip("Gemini API key not configured")
        return GeminiExtractor(api_key=api_key)
    
    def test_live_pubmed_search(self, pubmed_api):
        results = pubmed_api.search_pubmed("cancer treatment 2023", max_results=5)
        assert len(results) > 0
        assert all(isinstance(pmid, str) for pmid in results)
    
    def test_live_article_fetch(self, pubmed_api):
        article = pubmed_api.fetch_article_details("34567890")
        assert article is not None
        assert 'PubmedArticle' in article

    def test_live_gemini_extraction(self, gemini_api):
        test_text = """
        Study Population: 100 patients with diabetes.
        Methods: Randomized controlled trial.
        Results: 45% improvement in treatment group.
        Conclusion: Treatment was effective.
        """
        results = gemini_api.extract_sections(test_text)
        assert isinstance(results, dict)
        assert 'population' in results
        assert 'methods' in results
        assert 'results' in results
        assert 'conclusion' in results
        
    def test_live_gemini_complex_text(self, gemini_api):
        complex_text = """
        This study examined 250 participants across multiple centers.
        The methodology involved a double-blind protocol.
        Our findings indicate significant improvements.
        In conclusion, the treatment shows promise.
        """
        results = gemini_api.extract_sections(complex_text)
        assert isinstance(results, dict)
        assert any(results.values())
