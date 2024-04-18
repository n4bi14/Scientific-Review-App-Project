from AccessPubmedOA import AccessPubmedOA
from AccessPubmedPM import AccessPubmedPM
import GetFullText as oa
#import summaryGeneration as gemini

def main():
    chosen_articles = []
    count = 0

    user_query = input("Please enter in your keywords. If it includes a phrase, put quotation marks around the phrase: ")
    pm = AccessPubmedPM(user_query)
    
    chosen_articles = pm.findArticles()
    
    #summaryGenerator = gemini.SummaryGeneration()
    #summaryGenerator.generate_summaries(open_texts)
    
main()