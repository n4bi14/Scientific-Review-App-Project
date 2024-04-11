from AccessPubmed import AccessPubmed
import GetOpenAccess as oa

def main():
    user_query = input("Please enter in your keywords. If it includes a phrase, put quotation marks around the phrase: ")
    accessPubmed = AccessPubmed(user_query)
    
    accessPubmed.findArticles()
    
main()