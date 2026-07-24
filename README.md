# GenoSync — DNA Sequence Matching

GenoSync compares DNA sequences fetched from NCBI to estimate similarity between two organisms using global and local alignment methods.

<div style="font-family: Inter, -apple-system, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif; color: #333333;">
<h2 style="color:#00796B; margin-top:0;">Features (brief)</h2>
<ul>
    <li><strong>UI (<code>app.py</code>)</strong>: Streamlit interface for entering species and viewing results — makes the tool approachable.</li>
    <li><strong>Name resolution</strong> (<code>core.py</code>, <code>dna_matcher/fetch.py</code>): resolves common names to scientific names via NCBI Taxonomy — reduces ambiguity.</li>
    <li><strong>Marker priority chain</strong>: picks shared gene markers (COI, 16S, ITS, etc.) in order — ensures biologically valid comparisons.</li>
    <li><strong>Sequence fetching</strong> (Entrez / Biopython): retrieves GenBank records on demand — uses real, public data.</li>
    <li><strong>Caching</strong>: short-term cache of fetched sequences — speeds repeated queries and reduces API load.</li>
    <li><strong>Compare-length cap</strong>: trims comparisons (max 1000bp) and uses the smaller sequence length — bounds runtime and makes percentages fair.</li>
    <li><strong>LCS helpers</strong> (<code>alg_lcs.py</code>): building blocks for similarity measures and Hirschberg.</li>
    <li><strong>Hirschberg (global)</strong> (<code>alg_hirschberg.py</code>): space-efficient full-sequence alignment → global similarity %.</li>
    <li><strong>Smith–Waterman (local)</strong> (<code>alg_smith_waterman.py</code>): finds best-matching sub-regions → local match and score.</li>
    <li><strong>Orchestration</strong> (<code>dna_matcher/compare.py</code>, <code>core.py</code>): decides markers, runs algorithms, formats results — centralizes logic.</li>
    <li><strong>Optional persistence</strong> (<code>dna_matcher/db.py</code> / sqlite): lightweight storage for caches or metadata.</li>
    <li><strong>Dependencies</strong> (<code>requirements.txt</code>): lists required Python packages for reproducible setup.</li>
</ul>
</div>

Quick start

1. Install dependencies:

```bash
pip install -r requirements.txt
```

2. Create a `.env` with your NCBI contact and optional key:

```
NCBI_EMAIL=you@example.com
NCBI_API_KEY=your_ncbi_api_key_here
```

3. Run the app:

```bash
streamlit run app.py
```

How it works (overview)
- Resolve common names to scientific names via NCBI Taxonomy.
- Try markers in priority order until both species have the same marker.
- Fetch sequences (cached), trim to 1000bp, then compute:
    - Hirschberg: global similarity percentage
    - Smith–Waterman: best local match and score

Project layout

```
app.py                  # Streamlit UI
core.py                 # fetching, caching, marker routing
alg_lcs.py              # LCS helpers
alg_hirschberg.py       # global alignment (Hirschberg)
alg_smith_waterman.py   # local alignment (Smith–Waterman)
requirements.txt        # dependencies
```

Common markers
- Animals: COI, 16S, 12S, cytb
- Plants: rbcL, matK
- Fungi: ITS, LSU, SSU

Notes & limitations
- Single-marker comparisons only (not whole-genome).
- Sequences capped at 1000bp for performance.
- Common names can be ambiguous—verify the scientific name shown.
- Cache is per Streamlit session; use an API key to avoid rate limits.

References
- NCBI GenBank, Biopython Entrez, Hirschberg (1975), Smith & Waterman (1981)

See source files for implementation details and parameters.
</div>
