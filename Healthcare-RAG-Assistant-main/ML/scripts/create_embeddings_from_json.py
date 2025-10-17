# ML/scripts/create_embeddings_from_json.py
import os
import json
from tqdm import tqdm
from sentence_transformers import SentenceTransformer
import faiss
import numpy as np
import pickle

# -----------------------------
# 1) Locate repo and JSON file
# -----------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))        # .../ML/scripts
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))  # .../Healthcare-RAG-Assistant-main/Healthcare-RAG-Assistant-main

print("SCRIPT_DIR:", SCRIPT_DIR)
print("REPO_ROOT:", REPO_ROOT)

# Candidate default location
candidate = os.path.join(REPO_ROOT, "ML", "data", "pubmed_articles.json")

found_path = None
if os.path.exists(candidate):
    found_path = candidate
else:
    # fallback: walk the repo to find the file
    print("Default path not found, searching repo for pubmed_articles.json ... (this may take a few seconds)")
    for root, dirs, files in os.walk(REPO_ROOT):
        if "pubmed_articles.json" in files:
            found_path = os.path.join(root, "pubmed_articles.json")
            break

if not found_path:
    print("\nERROR: pubmed_articles.json not found under repo root.")
    print("Checked default candidate:", candidate)
    # Helpful diagnostics:
    ml_data_dir = os.path.join(REPO_ROOT, "ML", "data")
    print("Does ML/data exist?:", os.path.exists(ml_data_dir))
    if os.path.exists(ml_data_dir):
        print("ML/data listing (first 50):", os.listdir(ml_data_dir)[:50])
    raise FileNotFoundError("pubmed_articles.json not found. Please ensure file is at ML/data/pubmed_articles.json")

print("Using pubmed JSON at:", found_path)

# -----------------------------
# 2) Load JSON
# -----------------------------
with open(found_path, "r", encoding="utf-8") as f:
    papers = json.load(f)

print(f"Loaded {len(papers)} articles from JSON.")

# -----------------------------
# 3) Combine title+abstract and create plain texts list
# -----------------------------
texts = []
for p in papers:
    title = p.get("title", "")
    abstract = p.get("abstract", "")
    combined = f"Title: {title}\n\nAbstract: {abstract}"
    texts.append(combined)

print("Prepared", len(texts), "documents for chunking/embedding.")

# -----------------------------
# 4) Simple chunking (split long strings into chunks of ~500 words)
# -----------------------------
def chunk_text(text, chunk_words=500, overlap=50):
    words = text.split()
    chunks = []
    i = 0
    while i < len(words):
        chunk = words[i:i+chunk_words]
        chunks.append(" ".join(chunk))
        i += chunk_words - overlap
    return chunks

all_chunks = []
for t in texts:
    ch = chunk_text(t, chunk_words=500, overlap=100)
    all_chunks.extend(ch)

print("Total chunks created:", len(all_chunks))

# -----------------------------
# 5) Generate embeddings locally (sentence-transformers)
# -----------------------------
print("Loading SentenceTransformer model (this runs locally)...")
model = SentenceTransformer("all-MiniLM-L6-v2")  # small, fast, free

# encode in batches
batch_size = 64
embeddings = []
for i in tqdm(range(0, len(all_chunks), batch_size)):
    batch_texts = all_chunks[i:i+batch_size]
    emb = model.encode(batch_texts, show_progress_bar=False)
    embeddings.append(emb)
embeddings = np.vstack(embeddings).astype("float32")
print("Embeddings shape:", embeddings.shape)

# -----------------------------
# 6) Save FAISS index + metadata
# -----------------------------
vector_dir = os.path.join(REPO_ROOT, "ML", "data", "vector_db")
os.makedirs(vector_dir, exist_ok=True)

index = faiss.IndexFlatL2(embeddings.shape[1])
index.add(embeddings)
faiss_index_path = os.path.join(vector_dir, "index.faiss")
faiss.write_index(index, faiss_index_path)

# Save metadata (chunks -> to map indices back to text)
meta_path = os.path.join(vector_dir, "metadata.pkl")
with open(meta_path, "wb") as f:
    pickle.dump({"chunks": all_chunks}, f)

print(f"Saved FAISS index to: {faiss_index_path}")
print(f"Saved metadata to: {meta_path}")
print("✅ Done.")
