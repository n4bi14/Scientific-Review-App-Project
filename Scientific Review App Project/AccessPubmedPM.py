from Bio import Entrez
import GetFullText as full_text

#Program Description:
#Allows access of Pubmed articles via E-Utilities (the common way of accessing Pubmed).
#This method allows the use of Pubmed's relevance scoring system, which is beneficial
#in identifying articles that are relevant to the user's search.

Entrez.email = "jacobg2071@gmail.com"

class AccessPubmedPM:
    
    def __init__(self, query):
         self.query = query
         self.valid_articles = []

    def findArticles(self):
        count = 0
        article_ids = ""
        
        handle = Entrez.esearch(db = "pmc", retmax = 3000, term = self.query)
        records = Entrez.read(handle)
        handle.close()
        
        data = self.getArticleInfo(records['IdList'])

        for record in data:
            
            if('body' in record):
                if('sec' in record['body']):
                    article_text = full_text.fetch_full_text_xml_single(record['body']['sec'])
                    self.valid_articles.append(article_text)
                    ++count
                
                    if(count >= 100):
                       break
            
    def getArticleInfo(self,articles):
        fetchHandle=Entrez.efetch(db="pmc",id=articles,rettype="full",retmode="xml")
        xml_data = Entrez.read(fetchHandle)
            
        return xml_data
