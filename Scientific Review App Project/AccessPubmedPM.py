from Bio import Entrez
import GetFullText as full_text
import csv

#Program Description:
#Allows access of Pubmed articles via E-Utilities (the common way of accessing Pubmed).
#This method allows the use of Pubmed's relevance scoring system, which is beneficial
#in identifying articles that are relevant to the user's search.

Entrez.email = "jacob.t.galyean@gmail.com"

class AccessPubmedPM:
    
    def __init__(self, query):
         self.query = query
         self.valid_articles = []
         self.csv_file_path = 'PMC-ids.csv'

    def findArticles(self):
        count = 0
        article_ids = []
        
        handle = Entrez.esearch(db = "pubmed", retmax = 3000, term = self.query)
        records = Entrez.read(handle)
        handle.close()
        pmidMap = self.create_pmid_to_pmcid_mapping()
        
        for record in records["IdList"]:
            if(record in pmidMap):
                article_ids.append(pmidMap.get(record)) 
        
        data = self.getArticleInfo(article_ids)
        
        for record in data:
            
            if('body' in record):
                if('sec' in record['body']):
                    article_title = record['front']['article-meta']['title-group']['article-title']
                    article_text = full_text.fetch_full_text_xml_single(record['body']['sec'])
                    
                    self.valid_articles.append('Title: ' + article_title + 'Body:' + article_text)
                    count += 1
                
                    if(count >= 100):
                       break
                   
        return self.valid_articles 
            
    def getArticleInfo(self,articles):
        fetchHandle=Entrez.efetch(db="pmc",id=articles,rettype="full",retmode="xml")
        xml_data = Entrez.read(fetchHandle, validate=False)
        
        return xml_data

    
    # Function to create PMC ID to PubMed ID mapping from a CSV file
    def create_pmid_to_pmcid_mapping(self):
        pmid_to_pmcid_mapping = {}  # Initialize an empty dictionary

        with open(self.csv_file_path, mode='r', newline='') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=',')  # Using comma as the delimiter
            for row in reader:
                pmcid = row['PMCID']  # Get PMC ID from the 'PMCID' column
                pmid = row['PMID']    # Get PubMed ID from the 'PMID' column
                if pmcid and pmid:    # Check if both PMC ID and PubMed ID are present
                    pmid_to_pmcid_mapping[pmid] = pmcid  # Add mapping to the dictionary

        return pmid_to_pmcid_mapping