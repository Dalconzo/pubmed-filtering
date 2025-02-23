import re
from typing import Dict, Optional

class PatternExtractor:
    def __init__(self):
        # Define section markers with common variations
        self.patterns = {
            'population': r'(?i)(study population|participants|subjects|sample size|cohort)[\s:]*(.*?)(?=\n\n|\.[A-Z])',
            'methods': r'(?i)(methods|methodology|study design)[\s:]*(.*?)(?=\n\n|\.[A-Z])',
            'results': r'(?i)(results|findings|outcomes)[\s:]*(.*?)(?=\n\n|\.[A-Z])',
            'conclusion': r'(?i)(conclusion|conclusions|summary)[\s:]*(.*?)(?=\n\n|\.[A-Z])'
        }
        
    def extract_section(self, text: str, section_type: str) -> Optional[str]:
        pattern = self.patterns.get(section_type)
        if not pattern:
            return None
            
        match = re.search(pattern, text, re.DOTALL)
        return match.group(2).strip() if match else None

    def extract_all_sections(self, text: str) -> Dict[str, str]:
        return {
            section: self.extract_section(text, section)
            for section in self.patterns.keys()
        }
