from Bio import Entrez
from tqdm import tqdm
import os
import time
    
    # 1️⃣ Set your email and (optional) API key
Entrez.email = "your_email@example.com"  # Replace with your email (required by NCBI)
Entrez.api_key = None  # Optional: replace with your API key for faster requests
    
    # 2️⃣ Create output directory
OUTPUT_DIR = "Healthcare-RAG-Assistant-main\ML\scripts\data\pubmed_articles"
os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # 3️⃣ Function to fetch PubMed article IDs based on a search term
def fetch_article_ids(query, max_results=10000):
        """Fetch up to max_results PubMed article IDs by paginating through API."""
        ids = []
        retstart = 0
        batch = 200  # NCBI allows 200 per request
        
        print(f"Fetching up to {max_results} article IDs for query: '{query}'")
        while retstart < max_results:
            handle = Entrez.esearch(
                db="pubmed", 
                term=query, 
                retmax=batch, 
                retstart=retstart
            )
            record = Entrez.read(handle)
            handle.close()
            
            batch_ids = record["IdList"]
            if not batch_ids:
                break
            ids.extend(batch_ids)
            
            retstart += batch
            time.sleep(0.5)  # be polite to NCBI
    
        print(f"Total article IDs fetched: {len(ids)}")
        return ids[:max_results]
    
    
    # 4️⃣ Function to fetch article details (title, abstract, etc.)
def fetch_article_details(id_list):
        ids = ",".join(id_list)
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
        data = handle.read()
        handle.close()
        return data
    
    # 5️⃣ Download articles in batches
def download_pubmed_articles(query="diabetes treatment", total_articles=10000, batch_size=100):
        print(f"Searching for PubMed articles about: {query}")
        all_ids = fetch_article_ids(query, max_results=total_articles)
        print(f"Found {len(all_ids)} articles. Downloading...")
    
        for i in tqdm(range(0, len(all_ids), batch_size)):
            batch_ids = all_ids[i:i+batch_size]
            data = fetch_article_details(batch_ids)
            filename = os.path.join(OUTPUT_DIR, f"pubmed_{i//batch_size + 1}.txt")
            with open(filename, "w", encoding="utf-8") as f:
                f.write(data)
            time.sleep(0.5)  # be polite to the API
    
        print(f"✅ Download complete! Files saved in: {OUTPUT_DIR}")
    
if __name__ == "__main__":
        # You can change the topic below
        download_pubmed_articles(query="Alzheimer disease treatment", total_articles=10000)
        