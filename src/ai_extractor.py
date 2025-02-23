import google.generativeai as genai
from typing import Dict

class GeminiExtractor:
    def __init__(self, api_key: str):
        genai.configure(api_key=api_key)
        self.model = genai.GenerativeModel('gemini-pro')
        
    def extract_sections(self, text: str) -> Dict[str, str]:
        prompt = """
        From the following research paper text, extract these sections:
        1. Population/Sample Size
        2. Methods
        3. Results
        4. Conclusions
        
        Format as JSON with these keys: population, methods, results, conclusion
        If a section is not found, return null for that key.
        
        Paper text:
        {text}
        """
        
        response = self.model.generate_content(prompt.format(text=text))
        return response.text  # Parse JSON response in actual implementation
