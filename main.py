from AccessPubmedOA import AccessPubmedOA
from AccessPubmedPM import AccessPubmedPM
import GetFullText as oa
import random
import SummaryGeneration as gemini

def main():
    chosenArticles = []

    user_query = input("Please enter in your keywords. If it includes a phrase, put quotation marks around the phrase: ")
    pm = AccessPubmedPM(user_query)
    
    chosen_articles = pm.findArticles()
    
    #for i in range(0,100):
       #chosenArticles.append(random.choice(accessPubmedOA.validArticles))
    
    #open_texts = oa.fetch_full_text_xml(chosenArticles)
    
    #summaryGenerator = gemini.summaryGeneration()
    #summaryGenerator.generate_summaries(open_texts)
    
main()