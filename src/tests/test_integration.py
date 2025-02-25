from src.config import Config
from src.api_handler import PubMedAPI
from src.extractor_manager import ExtractorManager
from src.output_formatter import OutputFormatter
import pytest
import os


def test_full_workflow():
    # Initialize components
    config = Config()
    pubmed_api = PubMedAPI(config.get_pubmed_credentials())
    extractor_manager = ExtractorManager(config.get_gemini_key())
    output_formatter = OutputFormatter()
    # Search and fetch articles
    print("\n=== STARTING PUBMED SEARCH ===")
    pmids = pubmed_api.search_pubmed("cancer treatment 2023", max_results=2)
    print(f"\nFound {len(pmids)} articles with PMIDs:", pmids)
    articles = []
    
    # Process each article
    for pmid in pmids:
        print(f"\n=== PROCESSING ARTICLE {pmid} ===")
        article_data = pubmed_api.fetch_article_details(pmid)
        
        # Get metadata first
        metadata = pubmed_api.get_article_metadata(article_data)
        
        # Use abstract from metadata
        article_text = metadata['abstract']
        
        print("\nARTICLE ABSTRACT:")
        print("-" * 80)
        print(article_text)
        print("-" * 80)
        
        extracted_sections = extractor_manager.extract_with_fallback(article_text)
        
        print("\nEXTRACTED SECTIONS:")
        print("=" * 80)
        for section, content in extracted_sections.items():
            print(f"\n{section.upper()}:")
            print(content if content else "Not found")
            print("-" * 40)
        
        articles.append({
            'metadata': metadata,
            'extracted_sections': extracted_sections
        })
    output_path = os.path.expanduser("~/Documents/pubmed_filtering_tests/test_output.txt")
    print(f"\n=== WRITING OUTPUT TO: {output_path} ===")
    
    # Generate output with full path
    output_formatter.write_output(
        articles=articles,
        query="cancer treatment 2023",
        output_path=output_path,
        include_metadata=True
    )    
    # Add verification
    assert os.path.exists(output_path)
    with open(output_path, 'r') as f:
        content = f.read()
        assert "cancer treatment 2023" in content
        assert len(content) > 0 