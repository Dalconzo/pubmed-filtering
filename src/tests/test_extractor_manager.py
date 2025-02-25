import pytest
from src.extractor_manager import ExtractorManager
from src.config import Config

@pytest.fixture
def manager():
    config = Config()
    return ExtractorManager(config.get_gemini_key())

@pytest.fixture
def clear_text():
    return """
    Study Population: 100 diabetic patients
    Methods: Double-blind study
    Results: 45% improvement
    Conclusion: Treatment effective
    """

@pytest.fixture
def ambiguous_text():
    return """
    The research involved 100 patients.
    We conducted several tests.
    Improvements were noted.
    This suggests positive outcomes.
    """

def test_pattern_extraction(manager, clear_text):
    results = manager.extract_with_fallback(clear_text)
    assert results['population'] == '100 diabetic patients'
    assert results['methods'] == 'Double-blind study'
    assert results['results'] == '45% improvement'
    assert results['conclusion'] == 'Treatment effective'

def test_ai_fallback(manager, ambiguous_text):
    results = manager.extract_with_fallback(ambiguous_text)
    assert any(results.values()), "AI should extract at least some sections"
    assert isinstance(results, dict)
    assert all(key in results for key in ['population', 'methods', 'results', 'conclusion'])

def test_mixed_extraction(manager):
    mixed_text = """
    Study Population: 200 participants
    The study methodology varied by center.
    Results: Significant improvement observed.
    The findings suggest further research is needed.
    """
    results = manager.extract_with_fallback(mixed_text)
    assert results['population'] == '200 participants'
    assert 'Significant improvement observed' in results['results']
