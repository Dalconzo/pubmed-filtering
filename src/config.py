import os
from dotenv import load_dotenv

class Config:
    def __init__(self):
        load_dotenv()  # Load environment variables from .env file
        
        self.PUBMED_API_KEY = os.getenv('PUBMED_API_KEY')
        self.GEMINI_API_KEY = os.getenv('GEMINI_API_KEY')
        self.EMAIL = os.getenv('EMAIL')  # Required for PubMed API
        
    def get_pubmed_credentials(self):
        return {
            'email': self.EMAIL,
            'api_key': self.PUBMED_API_KEY
        }
        
    def get_gemini_key(self):
        return self.GEMINI_API_KEY
