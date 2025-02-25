import re
from typing import Dict, Optional

class PatternExtractor:
    def __init__(self):
        self.patterns = {
            'population': r'(?i)(study population|participants|subject population|sample size|cohort)[\s:]+(.+)',
            'methods': r'(?i)(methods|methodology|study design)[\s:]+(.+)',
            'results': r'(?i)(results|findings|outcomes)[\s:]+(.+)',
            'conclusion': r'(?i)(conclusion|conclusions|summary)[\s:]+(.+)'
        }

    def extract_section(self, text: str, section_type: str) -> Optional[str]:
        if not text or not isinstance(text, str):
            return None
            
        pattern = self.patterns.get(section_type)
        if not pattern:
            return None
            
        match = re.search(pattern, text)
        if match:
            extracted = match.group(2).strip()
            return extracted.rstrip('.')
        return None
    def extract_all_sections(self, text: str) -> Dict[str, str]:
        return {
            section: self.extract_section(text, section)
            for section in self.patterns.keys()
        }
