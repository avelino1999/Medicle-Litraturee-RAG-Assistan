import pickle

METADATA_PATH = r"C:\Users\aveli\Downloads\Healthcare-RAG-Assistant-main\Healthcare-RAG-Assistant-main\ML\data\vector_db\metadata.pkl"

with open(METADATA_PATH, "rb") as f:
    metadata_old = pickle.load(f)

if "chunks" in metadata_old:
    chunks = metadata_old["chunks"]
    metadata_new = {i: {"text": text} for i, text in enumerate(chunks)}
    
    with open(METADATA_PATH, "wb") as f:
        pickle.dump(metadata_new, f)
    print(f"✅ Fixed metadata format. Total entries: {len(metadata_new)}")
else:
    print("⚠️ No 'chunks' key found — metadata already in correct format.")
