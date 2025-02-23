
from src.ai_extractor import GeminiExtractor
from src.pattern_extractor import PatternExtractor
from typing import Dict

class ExtractorManager:
    def __init__(self, gemini_api_key: str):
        self.pattern_extractor = PatternExtractor()
        self.ai_extractor = GeminiExtractor(gemini_api_key)
        
    def extract_with_fallback(self, text: str) -> Dict[str, str]:
        # Try pattern matching first
        pattern_results = self.pattern_extractor.extract_all_sections(text)
        
        # For any missing sections, try AI extraction
        missing_sections = [k for k, v in pattern_results.items() if v is None]
        if missing_sections:
            ai_results = self.ai_extractor.extract_sections(text)
            # Merge results, preferring pattern matches
            pattern_results.update({k: ai_results[k] for k in missing_sections})
            
        return pattern_results
