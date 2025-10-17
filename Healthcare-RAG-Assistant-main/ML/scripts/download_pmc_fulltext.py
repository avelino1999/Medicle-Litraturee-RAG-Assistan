# ml/scripts/download_pmc_fulltext_fixed.py
import os
import re
import json
import time
from Bio import Entrez
from tqdm import tqdm
import xml.etree.ElementTree as ET

# config
Entrez.email = "avelinomonteiro02@gmail.com"   # <-- set your email
Entrez.api_key = "275b72c6dfbfc360e7293cf0696df2c40e09"                    # <-- optional, set your NCBI API key if you have one

# Where to save (absolute path based on this script file)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))  # two levels up to repo root
OUT_DIR = os.path.join(PROJECT_ROOT, "ml", "data", "processed")
OUT_FILE = os.path.join(OUT_DIR, "pmc_fulltext.json")
DEBUG_DIR = os.path.join(OUT_DIR, "debug_xml")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(DEBUG_DIR, exist_ok=True)

def remove_namespace(xml_str: str) -> str:
    # remove default namespace declarations and prefixed namespaces to make parsing easier
    xml_str = re.sub(r'\sxmlns(:\w+)?="[^"]+"', '', xml_str)
    xml_str = re.sub(r'<(/?)[a-zA-Z0-9]+:', r'<\1', xml_str)
    return xml_str

def parse_pmc_xml(xml_str: str):
    """
    Returns list of dicts: {'title':..., 'full_text':...}
    """
    papers = []
    if not xml_str or not xml_str.strip():
        return papers

    xml_str = remove_namespace(xml_str)
    try:
        root = ET.fromstring(xml_str)
    except Exception as e:
        return papers

    # Find all <article> nodes
    for article in root.findall('.//article'):
        # Title: article-title (may be nested)
        title_elem = article.find('.//article-title')
        title = ""
        if title_elem is not None:
            title = "".join(title_elem.itertext()).strip()

        # Body paragraphs: look for <body> then <p> tags
        paras = []
        for p in article.findall('.//body//p'):
            text = "".join(p.itertext()).strip()
            if text:
                paras.append(text)

        full_text = "\n\n".join(paras).strip()
        if title and full_text:
            papers.append({"title": title, "full_text": full_text})
    return papers

def fetch_pmc_fulltext(query="diabetes treatment", total_articles=200, batch_size=20, sleep=0.4):
    print(f"Searching for '{query}' in PMC...")
    search_handle = Entrez.esearch(db="pmc", term=query, retmax=total_articles)
    search_results = Entrez.read(search_handle)
    id_list = search_results.get("IdList", [])
    print(f"Found {len(id_list)} open-access articles (PMC IDs).")

    papers = []
    batch_count = 0
    for start in tqdm(range(0, len(id_list), batch_size)):
        batch = id_list[start:start + batch_size]
        batch_count += 1
        ids_str = ",".join(batch)
        try:
            handle = Entrez.efetch(db="pmc", id=ids_str, rettype="full", retmode="xml")
            xml_str_bytes = handle.read()
            xml_str = xml_str_bytes.decode("utf-8")
            handle.close()
        except Exception as e:
            print(f"Warning: efetch failed for batch {batch_count}: {e}")
            time.sleep(sleep)
            continue

        # parse xml -> articles
        parsed = parse_pmc_xml(xml_str)
        if not parsed:
            # Save raw xml for debugging so you can inspect why parsing failed
            debug_file = os.path.join(DEBUG_DIR, f"pmc_batch_{batch_count}.xml")
            with open(debug_file, "w", encoding="utf-8") as f:
                f.write(xml_str if isinstance(xml_str, str) else xml_str.decode("utf-8", errors="ignore"))
        else:
            papers.extend(parsed)

        time.sleep(sleep)

    # Save results
    with open(OUT_FILE, "w", encoding="utf-8") as f:
        json.dump(papers, f, indent=2, ensure_ascii=False)

    print(f"✅ Saved {len(papers)} full-text papers to: {OUT_FILE}")
    if os.listdir(DEBUG_DIR):
        print(f"⚠️  Some batches couldn't be parsed; debug XML saved to: {DEBUG_DIR}")
    return papers

if __name__ == "__main__":
    # Example: increase total_articles or split queries by year if you want more
    fetched = fetch_pmc_fulltext(query="Alzheimer disease biomarkers", total_articles=200, batch_size=20)
