import re
from Bio import Entrez

# Set your email address (required by Entrez)
Entrez.email = "jacob.t.galyean@gmail.com"

def fetch_full_text_xml(pmc_id):
    xml_str = ''
    article_texts = []
    
    handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="full", retmode="xml")
    xml_data = Entrez.read(handle)
    
    for record in xml_data:
        xml_str = ''
        xml_body = record['body']['sec']
    
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
        
        article_texts.append(clean_text)
    
    handle.close()
    
    return article_texts