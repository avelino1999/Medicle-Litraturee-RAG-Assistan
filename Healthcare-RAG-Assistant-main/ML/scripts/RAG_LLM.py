import os
import pickle
import faiss
import numpy as np
import google.generativeai as genai
from sentence_transformers import SentenceTransformer

# ===========================
# 📁 PATHS
# ===========================
BASE_DIR = os.path.dirname(os.path.dirname(__file__))  # ML/
VECTOR_DB_DIR = os.path.join(BASE_DIR, "data", "vector_db")
INDEX_PATH = os.path.join(VECTOR_DB_DIR, "index.faiss")
METADATA_PATH = os.path.join(VECTOR_DB_DIR, "metadata.pkl")

# ===========================
# 📦 LOAD FAISS INDEX + METADATA
# ===========================
if os.path.exists(INDEX_PATH):
    index = faiss.read_index(INDEX_PATH)
    print(f"✅ Loaded FAISS index from: {INDEX_PATH}")
else:
    raise FileNotFoundError(f"❌ FAISS index not found at {INDEX_PATH}")

if os.path.exists(METADATA_PATH):
    with open(METADATA_PATH, "rb") as f:
        metadata = pickle.load(f)
    print(f"✅ Loaded metadata from: {METADATA_PATH}")
else:
    raise FileNotFoundError(f"❌ Metadata file not found at {METADATA_PATH}")

# ===========================
# 🔑 CONFIGURE GEMINI API
# ===========================
genai.configure(api_key=os.environ["GOOGLE_API_KEY"])

# Main Gemini model for answers
LLM_MODEL = "models/gemini-2.5-flash"
# Embedding model for query encoding
EMBED_MODEL = "all-MiniLM-L6-v2"

# Initialize the SBERT model globally
sbert_model = SentenceTransformer(EMBED_MODEL)

print("✅ Gemini 2.5 Flash model loaded successfully.")


# ===========================
# 🧠 HELPER: Generate Embedding
# ===========================
def get_embedding(text: str):
    """Generate text embedding using the SentenceTransformer model."""
    # OLD: result = genai.embed_content(model=EMBED_MODEL, content=text)
    # OLD: return np.array(result["embedding"], dtype=np.float32)
    
    # NEW: Use the SBERT model for consistent dimensions
    return sbert_model.encode(text, convert_to_numpy=True, normalize_embeddings=True, show_progress_bar=False).astype(np.float32)


# ===========================
# 🧩 MAIN RAG FUNCTION
# ===========================
def ask_rag(query: str, top_k: int = 5):
    """Retrieve top similar chunks from FAISS and ask Gemini for answer."""
    # Create embedding for query
    query_emb = get_embedding(query).reshape(1, -1)

    # Search in FAISS
    distances, indices = index.search(query_emb, top_k)

    # Retrieve top context
    chunks = []
    for idx in indices[0]:
        if idx < len(metadata):
            chunks.append(metadata[idx]["text"])
    context = "\n\n".join(chunks)

    # Build prompt
    prompt = f"""
You are a helpful medical AI assistant. Use the context below to answer the question clearly and accurately.

Context:
{context}

Question: {query}
Answer:
"""

    # Generate response using Gemini
    model = genai.GenerativeModel(LLM_MODEL)
    response = model.generate_content(prompt)

    return response.text


# ===========================
# 🚀 MAIN SCRIPT
# ===========================
if __name__ == "__main__":
    question = input("Ask your question: ")
    answer = ask_rag(question)
    print("\nAnswer:", answer)
