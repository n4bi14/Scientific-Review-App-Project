import os
import google.generativeai as genai 
from AccessPubmedPM import AccessPubmedPM

#Program Description:
#This program utilizes the Gemini AI model to generate summaries of academic articles retrieved from PubMed. 
#It takes relevant articles based on the user's query and returns a summary of each article. 
#An API key is required for access to the generative model.


# Set up API key and model
os.environ["API_KEY"] = ""
genai.configure(api_key=os.environ["API_KEY"])
model = genai.GenerativeModel('gemini-pro')

class SummaryGeneration:
    @staticmethod
    def generate_summaries(articles):
        for index, article in enumerate(articles, start=1):
            # Split the article into title and body
            article_parts = article.split('Body:')
            article_title = article_parts[0].replace('Title: ', '')  # Extract title
            article = article_parts[1]  # Extract body
            
            responseFullText = model.generate_content(f"Please generate a one to two paragraph summary of the paper. Paragraphs only. Do not include any special formatting.': {article}")
            summaryFullText = responseFullText.text
            
            print(f"{index}. {article_title} ")
            print(f"Summary: {summaryFullText}")
            print() 
            