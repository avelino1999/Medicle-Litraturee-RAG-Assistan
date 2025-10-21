🧠 Healthcare RAG Assistant (Gemini-Powered)

A Retrieval-Augmented Generation (RAG) system designed to serve as a Healthcare AI Assistant, combining Google Gemini’s reasoning capabilities with local vector embeddings for contextual retrieval.
This project allows you to ask medical and healthcare-related queries with relevant document-based context.

🚀 Features

Uses Google Gemini API as the main LLM.

Retrieves relevant context using FAISS vector search.

Designed for healthcare-related RAG queries.

Fully open-source and extendable for medical NLP research.

🧩 Project Structure
Healthcare-RAG-Assistant/
│
├── ML/
│   ├── scripts/
│   │   ├── create_embeddings.py      # (Optional) Script to generate embeddings if needed
│   │   ├── RAG_LLM.py                # Main RAG pipeline (uses Gemini)
│   │   ├── check_api_key_ver.py      # Verify Gemini API key setup
│   │
│   └── vector_db/
│       ├── index.faiss               # FAISS index for vector search
│       ├── metadata.json             # Metadata mapping documents to embeddings
│
├── requirements.txt
├── README.md
└── .env                              # Environment file for Gemini API key

⚙️ Setup Instructions
1. Clone the Repository
git clone https://github.com/<your-username>/Healthcare-RAG-Assistant.git
cd Healthcare-RAG-Assistant

2. Create a Conda Environment
conda create -n healthcare_rag python=3.10
conda activate healthcare_rag

3. Install Dependencies
pip install -r requirements.txt

4. Setup Gemini API Key

Create a .env file in the project root with your Gemini API key:

GOOGLE_API_KEY=your_google_gemini_api_key_here


Or set it manually in your terminal:

setx GOOGLE_API_KEY "your_google_gemini_api_key_here"  # Windows
export GOOGLE_API_KEY="your_google_gemini_api_key_here"  # macOS/Linux

5. Verify Your API Key

Run the following script to confirm Gemini is accessible:

python ML/scripts/check_api_key_ver.py


✅ If it prints available Gemini models, your setup is correct.

6. Run the RAG Assistant

You can now run the main script that uses Gemini and your FAISS vector database:

python ML/scripts/RAG_LLM.py


If everything is configured properly, it will:

Load your existing vector_db/index.faiss and metadata.json.

Accept a query (e.g., “What are common treatments for diabetes?”).

Retrieve relevant documents from your local database.

Generate an informed response using Gemini.

7. (Optional) Regenerate Embeddings

If you modify your data or want to rebuild the FAISS database:

python ML/scripts/create_embeddings.py


This will recreate index.faiss, metadata.json, and optionally embeddings.pkl.

🧠 Example Output

User Input:

“What are the symptoms of hypertension?”

Model Response:

“Common symptoms include headaches, dizziness, and blurred vision. However, many individuals may remain asymptomatic until complications arise.”

📚 Technologies Used

Google Gemini API (LLM)

LangChain for RAG orchestration

FAISS for vector similarity search

Python-dotenv for environment management

💡 Future Enhancements

Add Gradio-based web interface

Support PDF document uploads

Enable local caching for offline use

Add summarization and medical citation support

🧾 License

This project is open-source and available under the MIT License.
