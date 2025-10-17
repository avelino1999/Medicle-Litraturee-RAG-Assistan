def normalize_query(query: str) -> str:
    query = query.lower().strip()

    synonyms = {
        "heart attack": "myocardial infarction",
        "alzheimers": "alzheimer's disease",
        "high blood pressure": "hypertension",
    }

    for k, v in synonyms.items():
        if k in query:
            query = query.replace(k, v)

    return query

if __name__ == "__main__":
    sample_query = "heart attack prevention"
    print("Normalized:", normalize_query(sample_query))
