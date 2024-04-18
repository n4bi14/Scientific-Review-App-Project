from AccessPubmedOA import AccessPubmedOA
from AccessPubmedPM import AccessPubmedPM
import GetFullText as oa
import SummaryGeneration as gemini

def main():
    chosenArticles = []

    user_query = input("Please enter in your keywords. If it includes a phrase, put quotation marks around the phrase: ")
    pm = AccessPubmedPM(user_query)
    
    chosen_articles = pm.findArticles()
    
    summaryGenerator = gemini.summaryGeneration()
    summaryGenerator.generate_summaries(chosen_articles)
    
main()