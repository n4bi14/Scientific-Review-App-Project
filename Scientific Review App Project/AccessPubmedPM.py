from Bio import Entrez

#Program Description:
#Allows access of Pubmed articles via E-Utilities (the common way of accessing Pubmed).
#This method allows the use of Pubmed's relevance scoring system, which is beneficial
#in identifying articles that are relevant to the user's search.

Entrez.email = "jacobg2071@gmail.com"

class AccessPubmedPM:
    
    def __init__(self, query):
         self.query = query

    def findArticles(self):
        handle = Entrez.esearch(db = "pmc", retmax = 3000, term = self.query)
        records = Entrez.read(handle)
        handle.close()
        
        for record in records:
            self.getArticleIds(records)
            
    def getArticleIds(records):
        
