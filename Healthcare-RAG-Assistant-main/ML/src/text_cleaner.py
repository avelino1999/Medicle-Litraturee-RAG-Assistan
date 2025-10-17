import re

def clean_text(text: str) -> str:
    if not text:
        return ""
    text = re.sub(r"<.*?>", " ", text)  # remove HTML tags
    text = re.sub(r"\[[0-9]+\]", " ", text)  # remove references like [1]
    text = re.sub(r"\s+", " ", text)  # normalize whitespace
    text = re.sub(r"[^a-zA-Z0-9.,;:()\-– ]", " ", text)  # remove special chars
    text = text.strip()
    return text

if __name__ == "__main__":
    sample = "We studied Alzheimer's disease [1] using <b>deep learning</b> methods."
    print("Before:", sample)
    print("After:", clean_text(sample))
