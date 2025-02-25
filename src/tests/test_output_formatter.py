import pytest
import os
from datetime import datetime
from src.output_formatter import OutputFormatter

@pytest.fixture
def formatter():
    return OutputFormatter()

@pytest.fixture
def sample_metadata():
    return {
        'title': 'Test Article',
        'authors': ['Smith, John', 'Doe, Jane'],
        'journal': 'Test Journal',
        'pmid': '12345',
        'year': '2023'
    }

@pytest.fixture
def sample_sections():
    return {
        'population': 'Test population data',
        'methods': 'Test methodology',
        'results': 'Test results',
        'conclusion': 'Test conclusion'
    }

def test_format_article(formatter, sample_metadata, sample_sections):
    formatted = formatter.format_article(sample_metadata, sample_sections)
    assert 'Test Article' in formatted
    assert 'Smith, John' in formatted
    assert 'Test Journal' in formatted
    assert all(section in formatted for section in sample_sections.values())

def test_format_article_compact(formatter, sample_metadata, sample_sections):
    compact = formatter.format_article_compact(sample_metadata, sample_sections)
    assert '[12345]' in compact
    assert len(compact.split('\n')) < len(formatter.format_article(sample_metadata, sample_sections).split('\n'))

def test_validate_output_path(formatter, tmp_path):
    valid_path = tmp_path / "output.txt"
    assert formatter.validate_output_path(str(valid_path))
    
    invalid_path = "/nonexistent/directory/output.txt"
    assert not formatter.validate_output_path(invalid_path)

def test_validate_article_data(formatter, sample_metadata, sample_sections):
    valid_article = {
        'metadata': sample_metadata,
        'extracted_sections': sample_sections
    }
    assert formatter.validate_article_data(valid_article)
    
    invalid_article = {'metadata': sample_metadata}
    assert not formatter.validate_article_data(invalid_article)

def test_write_output(formatter, tmp_path, sample_metadata, sample_sections):
    output_file = tmp_path / "test_output.txt"
    articles = [{
        'metadata': sample_metadata,
        'extracted_sections': sample_sections
    }]
    
    formatter.write_output(
        articles=articles,
        query="test query",
        output_path=str(output_file),
        include_metadata=True
    )
    
    assert output_file.exists()
    content = output_file.read_text()
    assert "PubMed Search Results" in content
    assert "test query" in content
    assert "Test Article" in content

def test_export_summary(formatter, tmp_path, sample_metadata):
    summary_file = tmp_path / "summary.txt"
    articles = [{'metadata': sample_metadata}] * 3
    
    formatter.export_summary(articles, str(summary_file))
    
    assert summary_file.exists()
    content = summary_file.read_text()
    assert "SUMMARY STATISTICS" in content
    assert "Test Journal: 3 articles" in content
