
# 🧬 GenoSync — Evolutionary DNA Matcher

GenoSync is a bioinformatics web application that computes the evolutionary similarity between multiple organisms by comparing their gene marker from DNA sequences. Fetching live genetic data from the NCBI GenBank and Taxonomy databases, the application performs both global and local sequence alignments to generate similarity matrices, heatmaps, and neighbor-joining phylogenetic trees. 

WEBSITE: [GenoSync Streamlit App](https://genosync-dnamatcher.streamlit.app/)

## 🏗️ Core Architecture & Data Flow

GenoSync is built with a modular architecture separating the user interface, data fetching, and computational algorithms:

1. **Query Resolution:** User inputs (common names) are passed through `dna_matcher/core.py` to NCBI's Taxonomy API to resolve exact scientific names, reducing biological ambiguity.
2. **Marker Selection:** The system iterates through a prioritized biological marker chain (e.g., COI, 16S for animals; rbcL for plants) to find a shared genetic marker among all selected species.
3. **Data Retrieval & Caching:** Sequences are fetched via Biopython/Entrez. To mitigate NCBI API rate limits and timeouts, GenoSync implements a local SQLite caching layer (`dna_matcher_cache.sqlite3`). 
4. **Alignment Engine:** Sequences are capped (default 1000bp) and routed to the algorithmic core, which computes Longest Common Subsequence (LCS), Hirschberg (global), and Smith-Waterman (local) alignments.
5. **Visualization:** Streamlit renders the results into interactive pairwise metrics, heatmaps, and clustering trees.

---

## ✨ Detailed Features

### 🖥️ User Interface (`app.py`)
* **Dynamic Input Handling:** Supports flexible input of 2 to 5 species.
* **Interactive Visualizations:** Generates matplotlib-based heatmaps for global similarity and renders text-based neighbor-joining trees.
* **Pairwise Explorer:** Allows users to drill down into specific species pairs to view exact matching genetic segments (strings) and numerical scores.
* **Interview/Demo Safe Mode:** Automatically falls back to the pre-populated SQLite database for species selection if the live NCBI API times out or becomes unresponsive during live evaluations.

### 🧬 Data Engineering & Persistence (`dna_matcher/core.py`, `dna_matcher/db.py`, `dna_matcher/fetch.py`)
* **Scientific Name Resolution:** Cross-references user inputs with internal search caches and live NCBI taxonomy databases to enforce strict scientific naming conventions.
* **Marker Priority Chain:** Intelligently selects the most appropriate gene marker based on organism type (Animals: COI, 16S, cytb | Plants: rbcL, matK | Fungi: ITS).
* **Two-Tier Caching System:** 
  * *Search Cache:* Stores mappings of common names to scientific names.
  * *Sequence Cache:* Stores fetched FASTA sequences to bypass redundant API calls, guaranteeing high availability and performance.

### 🧮 Alignment Algorithms
* **Hirschberg Algorithm (Global Alignment - `dna_matcher/alg_hirschberg.py`):** A divide-and-conquer approach to sequence alignment that computes the global similarity percentage. It is space-efficient, operating in linear space while maintaining standard time complexity.
* **Smith-Waterman (Local Alignment - `dna_matcher/alg_smith_waterman.py`):** Identifies the most highly conserved local sub-regions between two sequences, providing a localized match score and percentage.
* **Longest Common Subsequence (LCS - `alg_lcs.py`):** Foundational helper functions used to construct the scoring matrices for the alignment algorithms.

---

## 🚀 Quick Start

### 1. Install Dependencies
Ensure you have Python 3.9+ installed, then run:
```bash
pip install -r requirements.txt

```

### 2. Environment Configuration

Create a `.env` file in the root directory to responsibly identify your requests to NCBI:

```env
NCBI_EMAIL=your.email@example.com
NCBI_API_KEY=your_optional_ncbi_api_key

```

### 3. Run the Application

```bash
streamlit run app.py

```

---

## 📂 Project Structure

```text
├── dna_matcher/
│   ├── __init__.py
│   ├── alg_hirschberg.py       # Global alignment implementation
│   ├── alg_smith_waterman.py   # Local alignment implementation
│   ├── algorithms.py           # Algorithm routing and execution
│   ├── compare.py              # Sequence comparison logic and tree building
│   ├── core.py                 # Orchestration: caching, marker routing, name resolution
│   ├── db.py                   # SQLite schema and database interaction logic
│   └── fetch.py                # NCBI Entrez/Biopython network calls
├── .env                        # Environment variables for NCBI API
├── .gitignore                  # Git ignore rules
├── alg_lcs.py                  # Longest Common Subsequence utilities
├── app.py                      # Main Streamlit application and UI
├── dna_matcher_cache.sqlite3   # Pre-populated cache for API fail-safes
├── README.md                   # Project documentation
├── requirements.txt            # Python dependencies
└── TEST.py                     # Testing scripts

```

---

## ⚠️ Notes & Limitations

* **Sequence Truncation:** To ensure reasonable processing times on standard hardware, sequences are strictly capped at 1000 base pairs during algorithmic comparison.
* **Single-Marker Constraint:** Evolutionary distance is calculated based on a single shared marker (e.g., Mitochondrial COI) rather than whole-genome comparisons.
* **Nomenclature Ambiguity:** While the app attempts to resolve common names (e.g., "bear" -> *Ursus arctos*), users should verify the resolved scientific name in the UI to ensure biological accuracy.
* **Evaluation Reliability:** The pre-cached SQLite database ensures reliable performance and minimizes API timeouts and busy server errors during evaluations or demonstrations.

---

## 📚 References & Acknowledgments

* **NCBI GenBank & Taxonomy:** [National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/)
* **Biopython Entrez:** [Biopython Documentation](https://biopython.org/docs/1.75/api/Bio.Entrez.html)
* **Algorithms:**
* Hirschberg, D. S. (1975). *A linear space algorithm for computing maximal common subsequences.* Communications of the ACM.
* Smith, T. F., & Waterman, M. S. (1981). *Identification of common molecular subsequences.* Journal of Molecular Biology.



```

```