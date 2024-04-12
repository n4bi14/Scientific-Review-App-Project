from Bio import Entrez
import ReadXML

# Set your email address (required by Entrez)
Entrez.email = "jacob.t.galyean@gmail.com"

def fetch_full_text_xml(pmc_id):
    handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="full", retmode="xml")
    xml_data = Entrez.read(handle)
    
    xml_body = xml_data[0]['body']
    print(xml_body)
    p_sections = xml_body['sec'][0]['p']
    combined_text = ' '.join(p_sections)

    # Replace special characters with their appropriate representations
    xml_str_cleaned = combined_text.replace('\xa0', ' ')
    print(xml_str_cleaned)
    ReadXML.cleanText(xml_str_cleaned)

    handle.close()
    
    return xml_data