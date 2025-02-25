from Bio import Entrez
from typing import List, Dict
import time

class PubMedAPI:
    def __init__(self, email: str, api_key: str = None):
        """Initialize PubMed API handler"""
        Entrez.email = email
        self.api_key = api_key
        if api_key:
            Entrez.api_key = api_key
        self.delay = 0.34  # Rate limit compliance (3 requests per second)

    def search_pubmed(self, query: str, max_results: int = 100) -> List[str]:
        """Search PubMed and return list of PMIDs"""
        time.sleep(self.delay)
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]

    def fetch_article_details(self, pmid: str) -> Dict:
        """Fetch full article details for a given PMID"""
        time.sleep(self.delay)
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
        article = Entrez.read(handle)
        handle.close()
        return article

    def batch_fetch_articles(self, pmids: List[str], batch_size: int = 50) -> List[Dict]:
        """Fetch multiple articles in batches"""
        articles = []
        for i in range(0, len(pmids), batch_size):
            batch = pmids[i:i + batch_size]
            time.sleep(self.delay)
            handle = Entrez.efetch(db="pubmed", id=batch, rettype="xml")
            articles.extend(Entrez.read(handle))
            handle.close()
        return articles

    def get_article_metadata(self, article: Dict) -> Dict:
        """Extract key metadata from article"""
        print("\nRAW PUBMED RESPONSE:")
        print("=" * 80)
        print(article)
        print("=" * 80)
        
        print("\nARTICLE STRUCTURE:")
        for key in article:
            print(f"\nKey: {key}")
            if isinstance(article[key], list):
                for item in article[key]:
                    print("\nSubkeys:")
                    for subkey in item:
                        print(f"- {subkey}")
        
        metadata = {
            'title': '',
            'authors': [],
            'abstract': '',
            'publication_date': '',
            'journal': '',
            'pmid': ''
        }
        
        try:
            article_data = article['PubmedArticle'][0]['MedlineCitation']['Article']
            print("\nARTICLE DATA STRUCTURE:")
            for key in article_data:
                print(f"- {key}")
                
            metadata['title'] = article_data['ArticleTitle']
            metadata['authors'] = [
                f"{author['LastName']} {author['ForeName']}" 
                for author in article_data.get('AuthorList', [])
            ]
            
            if 'Abstract' in article_data:
                abstract_texts = article_data['Abstract']['AbstractText']
                metadata['abstract'] = ' '.join(abstract_texts)
                    
            metadata['journal'] = article_data['Journal']['Title']
            metadata['pmid'] = article['PubmedArticle'][0]['MedlineCitation']['PMID']
            
        except (KeyError, IndexError) as e:
            print(f"Error parsing article metadata: {e}")
            
        return metadata
    
def create_api_instance(email: str, api_key: str = None) -> PubMedAPI:
    """Factory function to create PubMedAPI instance"""
    return PubMedAPI(email, api_key)
