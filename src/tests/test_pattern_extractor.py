import pytest
from src.pattern_extractor import PatternExtractor

@pytest.fixture
def extractor():
    return PatternExtractor()

@pytest.fixture
def sample_text():
    return """
    Study Population: The study included 100 patients with diabetes.
    
    Methods: We conducted a randomized controlled trial over 12 months.
    
    Results: Treatment group showed 45% improvement compared to control.
    
    Conclusion: The intervention was effective for diabetes management.
    """

def test_extract_population(extractor, sample_text):
    result = extractor.extract_section(sample_text, 'population')
    assert result is not None
    assert "100 patients" in result

def test_extract_all_sections(extractor, sample_text):
    results = extractor.extract_all_sections(sample_text)
    assert len(results) == 4
    assert all(key in results for key in ['population', 'methods', 'results', 'conclusion'])
    assert "randomized controlled trial" in results['methods']

def test_invalid_input(extractor):
    assert extractor.extract_section("", "population") is None
    assert extractor.extract_section(None, "population") is None

def test_missing_section(extractor):
    text = "This text has no clear section markers"
    results = extractor.extract_all_sections(text)
    assert all(value is None for value in results.values())

def test_pattern_variations(extractor):
    variations = [
        "Study Population: 50 subjects",
        "Participants: 50 subjects",
        "Subject Population: 50 subjects",
        "Sample Size: 50 subjects"
    ]
    for text in variations:
        result = extractor.extract_section(text, 'population')
        assert "50 subjects" in result
