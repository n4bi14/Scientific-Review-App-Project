import xml.etree.ElementTree as et
    
#Function getArticles takes an XML document and creates an array of article ids.
#It only returns the article ids that are considered 'valid' by the checkValidity function.
def getArticleIds(response):
    root = et.fromstring(response)

    # Extract article IDs, how many results came through and the resumptionTokenLink.
    #The resumptionTokenLink is used to get the next 1000 results.
    article_ids = [record.get('id') for record in root.findall('.//record')]
    
    resumption_link = root.find('.//resumption/link').attrib.get('href') if root.find('.//resumption/link') is not None else None

    # Extract total count of records returned
    total_count = int(root.find('.//records').attrib.get('total-count') if root.find('.//records') else 0)
        
    return article_ids, total_count, resumption_link
    
    