from AccessPubmed import AccessPubmed
import GetOpenAccess as oa
import random

def main():
    chosenArticles = []

    user_query = input("Please enter in your keywords. If it includes a phrase, put quotation marks around the phrase: ")
    accessPubmed = AccessPubmed(user_query)
    
    accessPubmed.findArticles()
    for i in range(0,100):
        chosenArticles.append(random.choice(accessPubmed.article_ids))
    
    open_texts = oa.fetch_full_text_xml(chosenArticles)
    
    print(open_texts)
    
main()