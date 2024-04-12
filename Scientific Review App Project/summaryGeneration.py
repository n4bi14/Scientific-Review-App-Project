import os
import google.generativeai as genai 
from AccessPubmed import AccessPubmed
from GetFullText import GetFullText

# Set up API key and model
os.environ["API_KEY"] = ""
genai.configure(api_key=os.environ["API_KEY"])
model = genai.GenerativeModel('gemini-pro')

class summaryGeneration:
    def __init__(self, access_pubmed, GetFullText):
        self.titles = access_pubmed.validArticles
        self.abstracts = access_pubmed.validArticleAbs
        self.full_text = GetFullText.article_texts

    def generate_summaries(self):
        for index, (title, abstract, article_texts) in enumerate(zip(self.titles, self.abstracts, self.full_text), start=1):
            responseAbs = model.generate_content(f"Please generate a 1-2 sentence max summary of the abstract for the article titled '{title}': {abstract}")
            responseFullText = model.generate_content(f"Please generate a one to two paragraph max summary of the paper titled '{title}': {article_texts}")

            summaryAbs = responseAbs.text
            summaryFullText = responseFullText.text

            print(f"Article {index}:")
            print(f"Title: {title}")
            print(f"Abstract Summary: {summaryAbs}")
            print(f"Paper Summary: {summaryFullText}")
            print()