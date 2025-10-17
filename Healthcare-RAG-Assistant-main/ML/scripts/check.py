import pickle

with open(r"C:\Users\aveli\Downloads\Healthcare-RAG-Assistant-main\Healthcare-RAG-Assistant-main\ML\data\vector_db\metadata.pkl", "rb") as f:
    metadata = pickle.load(f)

print("Number of entries in metadata:", len(metadata))
print("Example entry:", list(metadata.items())[:1])
