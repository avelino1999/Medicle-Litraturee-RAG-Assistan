from Bio import Entrez
import json
import time

Entrez.email = "avelinomonteiro02@gmail.com"

def fetch_pubmed_articles(query, retmax=100):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    record = Entrez.read(handle)
    ids = record["IdList"]
    handle.close()

    # Fetch full article summaries (title, abstract, etc.)
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="xml")
    papers = Entrez.read(handle)
    handle.close()

    results = []
    for article in papers["PubmedArticle"]:
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", [""])[0]
        results.append({"title": title, "abstract": abstract})
    return results

if __name__ == "__main__":
    query = "Alzheimer's disease AND machine learning"
    articles = fetch_pubmed_articles(query)
    with open("Healthcare-RAG-Assistant-main\ML\data\pubmed_articles.json", "w", encoding="utf-8") as f:
        json.dump(articles, f, indent=2, ensure_ascii=False)
    print(f"✅ Saved {len(articles)} articles to data/pubmed_articles.json")
 
  