import google.generativeai as genai
import json
from typing import Dict

class GeminiExtractor:
    def __init__(self, api_key: str):
        genai.configure(api_key=api_key)
        self.model = genai.GenerativeModel('gemini-pro')
        
    def extract_sections(self, text: str) -> Dict[str, str]:
        prompt = """
        From the following research paper text, extract this data:
        1. Population/Sample Size
        2. Methods
        3. Results
        4. Conclusions
        
        Format as JSON with these keys: population, methods, results, conclusion
        If a section is not found, return null for that key.
        Always return valid JSON even if no sections are found.
        
        Paper text:
        {text}
        """
        
        try:
            response = self.model.generate_content(prompt.format(text=text))
            # Extract JSON from markdown code block if present
            json_str = response.text.strip('`json\n`')
            return json.loads(json_str)
        except json.JSONDecodeError:
            # Return empty sections if JSON parsing fails
            return {
                'population': None,
                'methods': None,
                'results': None,
                'conclusion': None
            }