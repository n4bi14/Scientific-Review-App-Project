from Bio import Entrez
import ReadXML

# Set your email address (required by Entrez)
Entrez.email = "jacob.t.galyean@gmail.com"

def fetch_full_text_xml(pmc_id):
    handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="full", retmode="xml")
    xml_data = handle.read()
    handle.close()
    
    xml_data = ReadXML.cleanText(xml_data)
    
    return xml_data