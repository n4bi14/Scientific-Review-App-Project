from AccessPubmedOA import AccessPubmedOA
from AccessPubmedPM import AccessPubmedPM
import GetFullText as oa
import SummaryGeneration as gemini

def main():
    chosenArticles = []

    user_query = input("Please enter in your keywords. If it includes a phrase, put quotation marks around the phrase: ")
    pm = AccessPubmedPM(user_query)
    
    chosen_articles = pm.findArticles()
    
    chosenArticles.append((AccessPubmedOA.validArticles))
    
    open_texts = oa.fetch_full_text_xml_single(chosenArticles)
    
    summaryGenerator = gemini.SummaryGeneration()
    summaryGenerator.generate_summaries(open_texts)
    
main()