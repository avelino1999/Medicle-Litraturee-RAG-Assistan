import spacy

# Load the scispaCy model
nlp = spacy.load("en_core_sci_md")

def extract_medical_entities(text: str):
    doc = nlp(text)
    entities = [{"text": ent.text, "label": ent.label_} for ent in doc.ents]
    return entities

if __name__ == "__main__":
    sample = "Patients with Alzheimer's disease were treated using donepezil and memantine."
    print(extract_medical_entities(sample))
