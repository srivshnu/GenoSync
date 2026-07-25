```markdown
# 🧬 GenoSync — Evolutionary DNA Matcher

GenoSync compares DNA sequences fetched from NCBI to estimate similarity between two organisms using global and local alignment methods.

WEBSITE: [GenoSync Streamlit App](https://genosync-dnamatcher.streamlit.app/)

## ✨ Features

* **UI (`app.py`)**: Streamlit interface for entering species and viewing results — makes the tool approachable.
* **Name resolution** (`core.py`, `dna_matcher/fetch.py`): Resolves common names to scientific names via NCBI Taxonomy — reduces ambiguity.
* **Marker priority chain**: Picks shared gene markers (COI, 16S, ITS, etc.) in order — ensures biologically valid comparisons.
* **Sequence fetching** (Entrez / Biopython): Retrieves GenBank records on demand — uses real, public data.
* **Caching**: Local SQLite caching layer (`dna_matcher_cache.sqlite3`) and short-term caches to speed repeated queries, reduce API load, and bypass external API instability during live demonstrations.
* **Scientific Name Resolution**: Automatically resolves common names to precise scientific classifications using local search caches and NCBI taxonomy data.
* **Compare-length cap**: Trims comparisons (max 1000bp) and uses the smaller sequence length — bounds runtime and makes percentages fair.
* **LCS helpers** (`alg_lcs.py`): Building blocks for similarity measures and Hirschberg.
* **Hirschberg (global)** (`alg_hirschberg.py`): Space-efficient full-sequence alignment → global similarity %.
* **Smith–Waterman (local)** (`alg_smith_waterman.py`): Finds best-matching sub-regions → local match and score.
* **Orchestration** (`dna_matcher/compare.py`, `core.py`): Decides markers, runs algorithms, formats results — centralizes logic.
* **Persistence** (`dna_matcher/db.py` / sqlite): Lightweight storage for caches, pre-populated marker data, and metadata bundled for Streamlit Cloud deployment.
* **Dependencies** (`requirements.txt`): Lists required Python packages for reproducible setup.

---

## 🚀 Quick Start

1. Install dependencies:
   ```bash
   pip install -r requirements.txt

```

2. Create a `.env` with your NCBI contact and optional key:
```env
NCBI_EMAIL=you@example.com
NCBI_API_KEY=your_ncbi_api_key_here

```


3. Run the app:
```bash
streamlit run app.py

```



---

## 🔍 How It Works

* Resolve common names to scientific names via NCBI Taxonomy (with fallback to local search results).
* Try markers in priority order until all selected species share the same marker.
* Fetch sequences (cached locally in SQLite), trim to 1000bp, then compute:
* **Hirschberg**: Global similarity percentage and sequence alignment matching segments.
* **Smith–Waterman**: Best local match, alignment score, and sub-region matching percentage.
* **Neighbor-Joining**: Evolutionary clustering approximation tree.



---

## 📂 Project Layout

```text
app.py                  # Streamlit UI
core.py                 # Fetching, caching, marker routing, and scientific name resolution
alg_lcs.py              # LCS helpers
alg_hirschberg.py       # Global alignment (Hirschberg)
alg_smith_waterman.py   # Local alignment (Smith–Waterman)
requirements.txt        # Dependencies

```

---

## 🧬 Common Markers

* **Animals**: COI, 16S, 12S, cytb
* **Plants**: rbcL, matK
* **Fungi**: ITS, LSU, SSU

---

## ⚠️ Notes & Limitations

* Single-marker comparisons only (not whole-genome).
* Sequences capped at 1000bp for performance.
* Common names can be ambiguous—verify the scientific name shown.
* Pre-cached SQLite database ensures reliable performance and minimize API timeouts and busy server during evaluations.

---

## 📚 References

* NCBI GenBank, Biopython Entrez, Hirschberg (1975), Smith & Waterman (1981)

```

```