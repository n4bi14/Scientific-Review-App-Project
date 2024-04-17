import re
from Bio import Entrez

#Program Descritpion:
#When provided the Pubmed Central (PMC) ids, this program is able to retrieve
#the open access of all articles provided. It does this using E-Utilities.

Entrez.email = "jacob.t.galyean@gmail.com"

def fetch_full_text_xml_single(xml_body):
    xml_str = ''
    article_text = []
    
    for paragraph in xml_body:
        xml_str += str(paragraph['p'])

    #Remove XML tags
    clean_text = re.sub(r'<.*?>', '', xml_str)
    
    # Remove Unicode characters
    clean_text = re.sub(r'\\xa0', ' ', clean_text)
    
    # Remove unwanted punctuation
    clean_text = re.sub(r'[\[\]]', '', clean_text)
    
    # Remove leading and trailing whitespaces
    clean_text = clean_text.strip()
        
    article_text.append(clean_text)
    
    return article_text