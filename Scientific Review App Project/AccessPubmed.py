import re
import requests
from Bio import Entrez
import ReadXML as xml
import csv
import nltk
from nltk.tokenize import word_tokenize

nltk.download('punkt') 

Entrez.email = "jacob.t.galyean@gmail.com"

class AccessPubmed:
    def __init__(self, query):
        self.query = query
        self.queryTerms = re.findall(r'[^"\s]+|"[^"]*"', self.query)
        self.validArticles = []
        self.validArticleAbs = []
        
        self.csv_file_path = 'PMC-ids.csv'
        
    #Function to find articles from the Open Access Subset of Pubmed using the OA Web Service API.
    def findArticles(self):
        article_ids = []
        count = 0
        total_count = 9999
        #Get the link, combining what the user entered.
        oaLink = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?term=" + self.query 
        print(oaLink)
        
        while(count < 10):
            result = requests.get(oaLink)
        
            if (result.status_code == 200):
                new_article_ids, total_count, oaLink = xml.getArticleIds(result.text)
                print(total_count)
                print(oaLink)
                
                article_ids.extend(new_article_ids)
                
            else:
                print("Error:", result.status_code)

            count += 1
            print(count)
            
        if article_ids:
            pmcid_to_pmid_mapping = self.create_pmcid_to_pmid_mapping()
            for article in article_ids:
                article = pmcid_to_pmid_mapping[article]

            self.extractInfo(article_ids)
               
        else:
            print("No article IDs found")
               
    def extractInfo(self,article_ids):   
            fetchHandle=Entrez.efetch(db="pubmed",id=article_ids,rettype="null",retmode="xml")
            records = Entrez.read(fetchHandle)
                
            for record in records['PubmedArticle']:  # Loop through 'PubmedArticle' records
                
                if 'MedlineCitation' in record and 'Article' in record['MedlineCitation']:
                    article_data = record['MedlineCitation']['Article']
                    
                    if 'Abstract' in article_data:
                        abstract = article_data['Abstract']['AbstractText'][0]

                        if (self.checkValidity(abstract)):
                            print("Passed!")
                            self.validArticles.append(record['MedlineCitation']['PMID'])
                            self.validArticleAbs.append(abstract)
            
            fetchHandle.close()
            
    def checkValidity(self, abstract):
        abstract_lower = abstract.lower()  # Normalize abstract to lowercase
        query_terms_lower = [term.lower() for term in self.queryTerms]
        
        abstract_tokens = word_tokenize(abstract_lower)
        
        for term in query_terms_lower:
            if term in abstract_tokens:
                print("NICE!")
                return True
            
        return False
    
    # Function to create PMC ID to PubMed ID mapping from a CSV file
    def create_pmcid_to_pmid_mapping(self):
        pmcid_to_pmid_mapping = {}  # Initialize an empty dictionary

        with open(self.csv_file_path, mode='r', newline='') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=',')  # Using comma as the delimiter
            for row in reader:
                pmcid = row['PMCID']  # Get PMC ID from the 'PMCID' column
                pmid = row['PMID']    # Get PubMed ID from the 'PMID' column
                if pmcid and pmid:    # Check if both PMC ID and PubMed ID are present
                    pmcid_to_pmid_mapping[pmcid] = pmid  # Add mapping to the dictionary

        return pmcid_to_pmid_mapping
        

        