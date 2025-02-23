from datetime import datetime
from typing import Dict, List, Optional
import os
from collections import Counter

class OutputFormatter:
    def __init__(self):
        self.section_order = ['population', 'methods', 'results', 'conclusion']
        self.section_headers = {
            'population': 'STUDY POPULATION',
            'methods': 'METHODOLOGY',
            'results': 'RESULTS',
            'conclusion': 'CONCLUSIONS'
        }

        self.separator_char = '='
        self.subseparator_char = '-'

    def set_separator_style(self, main_char: str, sub_char: str) -> None:
        self.separator_char = main_char
        self.subseparator_char = sub_char

    def format_article(self, metadata: Dict, extracted_sections: Dict) -> str:
        lines = []
        lines.append(self.separator_char * 80)
        lines.append(f"TITLE: {metadata.get('title', 'N/A')}")
        lines.append(f"AUTHORS: {', '.join(metadata.get('authors', ['N/A']))}")
        lines.append(f"JOURNAL: {metadata.get('journal', 'N/A')}")
        lines.append(f"PMID: {metadata.get('pmid', 'N/A')}")
        lines.append(self.subseparator_char * 80)
        
        for section in self.section_order:
            content = extracted_sections.get(section)
            lines.append(f"\n{self.section_headers[section]}:")
            lines.append(content if content else "Not found")
            lines.append(self.subseparator_char * 40)
        
        return '\n'.join(lines)

    def format_article_compact(self, metadata: Dict, extracted_sections: Dict) -> str:
        lines = []
        lines.append(f"[{metadata.get('pmid', 'N/A')}] {metadata.get('title', 'N/A')}")
        
        for section in self.section_order:
            content = extracted_sections.get(section)
            if content:
                lines.append(f"{self.section_headers[section]}: {content[:200]}...")
        
        return '\n'.join(lines)

    def export_summary(self, articles: List[Dict], output_path: str) -> None:
        journals = Counter(article['metadata'].get('journal') for article in articles)
        years = Counter(article['metadata'].get('year') for article in articles)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("SUMMARY STATISTICS\n")
            f.write(f"Total Articles: {len(articles)}\n\n")
            
            f.write("Journal Distribution:\n")
            for journal, count in journals.most_common():
                f.write(f"{journal}: {count} articles\n")
                
            f.write("\nPublication Years:\n")
            for year, count in sorted(years.items()):
                f.write(f"{year}: {count} articles\n")

    def validate_output_path(self, path: str) -> bool:
        dir_path = os.path.dirname(path) or '.'
        return os.path.exists(dir_path) and os.access(dir_path, os.W_OK)

    def validate_article_data(self, article: Dict) -> bool:
        if not isinstance(article, dict):
            return False
        required_keys = {'metadata', 'extracted_sections'}
        return all(key in article for key in required_keys)

    def write_output(self, 
                    articles: List[Dict],
                    query: str,
                    output_path: str,
                    include_metadata: bool = True) -> None:
        if not self.validate_output_path(output_path):
            raise ValueError(f"Invalid output path: {output_path}")

        with open(output_path, 'w', encoding='utf-8') as f:
            if include_metadata:
                f.write(f"PubMed Search Results\n")
                f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Query: {query}\n")
                f.write(f"Articles found: {len(articles)}\n\n")
                f.write(self.separator_char * 80 + '\n\n')
            
            for article in articles:
                if self.validate_article_data(article):
                    formatted_article = self.format_article(
                        article['metadata'],
                        article['extracted_sections']
                    )
                    f.write(formatted_article + '\n\n')