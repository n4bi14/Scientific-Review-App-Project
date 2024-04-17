from AccessPubmedOA import AccessPubmedOA
import GetFullText as oa
import random
import summaryGeneration as gemini

def main():
    chosenArticles = []

    user_query = input("Please enter in your keywords. If it includes a phrase, put quotation marks around the phrase: ")
    accessPubmedOA = AccessPubmedOA(user_query)
    
    accessPubmedOA.findArticles()
    for i in range(0,100):
       chosenArticles.append(random.choice(accessPubmedOA.validArticles))
    
    open_texts = oa.fetch_full_text_xml(chosenArticles)
    
    summaryGenerator = gemini.summaryGeneration()
    summaryGenerator.generate_summaries(open_texts)
    
main()