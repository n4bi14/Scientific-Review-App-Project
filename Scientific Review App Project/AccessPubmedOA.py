import re
import requests
from Bio import Entrez
import ReadXML as xml
import csv

#Program Description:
#Allows the access of articles via the OA Web Service. This method returns only open access articles, at the cost of not
#being pre-relevance scored, returning any article that contains even one instance of each of the user's terms.

#After using the OA Web Service, it converts the PMC id to Pubmed ID and then asks E-Utilities for the abstract, 
#filters out the articles that don't have one of the user's terms in it's abstract, and then stores the 'valid articles'.

Entrez.email = "jacob.t.galyean@gmail.com"

class AccessPubmed:
    def __init__(self, query):
        self.query = query
        self.queryTerms = re.findall(r'[^"\s]+|"[^"]*"', self.query)
        self.article_ids = []
        self.validArticles = []
        self.validArticleAbs = []
        
        self.csv_file_path = 'PMC-ids.csv'
        
    #Function to find articles from the Open Access Subset of Pubmed using the OA Web Service API.
    def findArticles(self):
        pubmed_article_ids = []
        count = 0
        total_count = 9999
        #Get the link, combining what the user entered.
        oaLink = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?term=" + self.query 
        print(oaLink)
        
        #While loop that determines how many times articles are retrieved, currently set to 10 for testing purposes.
        while(count*1000 < total_count / 2 and count < 10):
            
            #Accesses the link using the requests module.
            result = requests.get(oaLink)
        
            #Checks if the result succeeds or not.
            if (result.status_code == 200):
                
                #Use the getArticleIds function in order to retrieve the following variables:
                #   - new_article_ids: Article ids that have not been retrieved yet.
                #   - total_count: The amount of articles the user's search results in.
                #   - oaLink: The link that goes to the next set of article ids.
                new_article_ids, total_count, oaLink = xml.getArticleIds(result.text)
                print(oaLink)
                
                self.article_ids.extend(new_article_ids)
                
                #Safeguard that allows the user to narrow down the search further by adding more terms if 40k articles or more are returned.
                if (total_count >= 40000 and count == 0):
                    print("Your search results in " + str(total_count) + " articles, you may want to add more terms or phrases to your search to narrow it down.")
                    print("The more articles there are, the longer the process will be and the less accurate the results may be.")
                    add_to_query = input("Please enter in other terms and phrases, if you would not like to add anymore just press ENTER: ")
                    
                    if (add_to_query != ""):
                        #If they choose to add queries, update the varaibles and link.
                        self.query = self.query + " " + add_to_query
                        self.queryTerms = re.findall(r'[^"\s]+|"[^"]*"', self.query)
                    
                        oaLink = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?term=" + self.query 
                        count -= 1
                        self.article_ids = []

            #If there is an error in fetching the article ids, print the error.
            else:
                print("Error:", result.status_code)

            count += 1
        
        if len(self.article_ids) > 0:
            #Change article_ids to go from a PMC id to a Pubmed id.
            pmcid_to_pmid_mapping = self.create_pmcid_to_pmid_mapping()
            
            for article in self.article_ids:
                pubmed_article_ids.append(pmcid_to_pmid_mapping.get(article, 'Not Found'))       

            self.extractInfo(pubmed_article_ids)
               
        else:
            print("No article IDs found")
    
            
    def extractInfo(self,article_ids):
            count = 0
            fetchHandle=Entrez.efetch(db="pubmed",id=article_ids,rettype="null",retmode="xml")
            records = Entrez.read(fetchHandle)
                
            for record in records['PubmedArticle']:  # Loop through 'PubmedArticle' records
                count += 1
                
                if 'MedlineCitation' in record and 'Article' in record['MedlineCitation']:
                    article_data = record['MedlineCitation']['Article']
                    
                    if 'Abstract' in article_data:
                        abstract = article_data['Abstract']['AbstractText'][0]

                        if (self.checkAbstract(abstract)):
                            self.validArticles.append(self.article_ids[count])
                            self.validArticleAbs.append(abstract)
            
            fetchHandle.close()
            
    def checkAbstract(self, abstract):
        #abstract_lower = abstract.lower()  # Normalize abstract to lowercase
        #query_terms_lower = [term.lower() for term in self.queryTerms]
        
        #abstract_tokens = word_tokenize(abstract_lower)
        
        for term in self.queryTerms:
            if term in abstract:
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
        

        