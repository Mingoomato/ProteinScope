# ProteinScope

**Comprehensive Protein Information Aggregator**

ProteinScope is a command-line tool that aggregates detailed protein information
from multiple public databases into a single structured report — including
sequences, function, clinical data, interaction domains, 3D structures, and
citable references, suitable for academic use.

---

## Data Sources

| Data | Source |
|---|---|
| Protein identity, function, sequences | UniProtKB REST API |
| DNA / CDS sequence | NCBI E-utilities (Gene, Nucleotide) |
| Clinical variants | NCBI ClinVar |
| Disease / deficiency | OMIM (API key required) |
| Temperature & pH optima | BRENDA (credentials required) |
| AlphaFold2 structure | AlphaFold DB |
| Experimental structures | RCSB PDB |
| References | PubMed / NCBI Efetch |
| Cross-species orthologs | UniProt BLAST + OrthoDB |

---

## Installation

### 1. Clone / copy the project

```bash
cd proteinscope
```

### 2. Create and activate a virtual environment

```bash
python -m venv .venv

# Windows
.venv\Scripts\activate

# macOS / Linux
source .venv/bin/activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

> **Note:** WeasyPrint (for PDF output) requires additional system libraries
> on Windows and Linux. See https://doc.courtbouillon.org/weasyprint/stable/first_steps.html

---

## API Key Setup

Copy the example env file and fill in your keys:

```bash
cp .env.example .env
```

Edit `.env`:

```
# NCBI (optional — increases rate limit to 10 req/sec)
NCBI_API_KEY=your_key_here

# OMIM (required for disease data)
# Register at: https://www.omim.org/api
OMIM_API_KEY=your_key_here

# BRENDA (required for temperature / pH optima)
# Register at: https://www.brenda-enzymes.org/register.php
BRENDA_USERNAME=your_email@example.com
BRENDA_PASSWORD=your_brenda_password
```

UniProt, AlphaFold DB, RCSB PDB, and PubMed do **not** require API keys.

---

## Usage

```bash
python main.py [OPTIONS] QUERY
```

`QUERY` can be:
- A UniProt accession: `P68871`
- A gene symbol: `HBB`
- A protein name: `"hemoglobin beta"`

### Options

| Option | Default | Description |
|---|---|---|
| `--organism TEXT` | — | Filter by organism, e.g. `"Homo sapiens"` |
| `--format` | `markdown` | Output format: `json`, `markdown`, or `pdf` |
| `--output / -o` | `./report` | Output file path (no extension) |
| `--cross-species` | off | Include cross-species binding domain comparison |
| `--no-cache` | off | Bypass local cache and always fetch fresh |
| `--threshold FLOAT` | `80.0` | Identity % threshold for cross-species compatibility |

---

## CLI Examples

```bash
# Basic search by UniProt accession — Markdown report
python main.py P68871 --format markdown --output ./reports/hemoglobin

# Search by protein name with organism filter — PDF
python main.py "EGFR" --organism "Homo sapiens" --format pdf --output ./reports/egfr

# Full report with cross-species binding domain comparison
python main.py "insulin receptor" --format pdf --cross-species --output ./reports/insr

# JSON dump for programmatic use
python main.py P04637 --format json --output ./reports/tp53

# Skip cache (always fetch fresh)
python main.py P00533 --no-cache --format markdown
```

---

## Output Formats

### Markdown (`.md`)
Structured human-readable report with all sections from the data model.

### JSON (`.json`)
Full `ProteinRecord` serialization with a metadata header (query time,
data sources, version). Suitable for downstream programmatic processing.

### PDF (`.pdf`)
Academic-style PDF rendered via WeasyPrint + Jinja2. Includes:
- Cover page (protein name, gene, organism, date)
- Numbered sections with page numbers
- Formatted APA reference list
- BibTeX block
- AlphaFold confidence score and direct links to RCSB / AlphaFold Server

---

## Report Sections

1. Identity
2. Function & Biological Role
3. Sequences (canonical AA, CDS DNA, isoforms)
4. Cofactors & Kinetics
5. Subcellular Location
6. Clinical Information (diseases, deficiency, ClinVar)
7. Interaction Domains (binding sites, active sites, domain map)
8. Cross-Species Applicability
9. 3D Structure (AlphaFold + experimental PDB)
10. References (APA, BibTeX)

---

## Running Tests

```bash
pytest
```

Tests use `pytest-httpx` fixtures — no live API calls are made in the test suite.

Test proteins used as fixtures:
- `P68871` — Hemoglobin subunit beta (diseases, isoforms, cofactors)
- `P00533` — EGFR receptor (binding sites, cross-species)
- `P04637` — TP53 tumor suppressor (clinical variants, domains)

---

## Caching

API responses are cached locally in `.cache/proteinscope.db` (SQLite) with a
7-day TTL. Use `--no-cache` to force fresh fetches.

---

## Notes & Limitations

- **AlphaFold:** The public DB hosts AF2 models. For AF3, use https://alphafoldserver.com
- **BRENDA:** Covers enzymes only. Non-enzyme proteins may have no temperature data.
- **ClinVar:** Human variants only.
- **OMIM rate limit:** Free academic key allows ~10,000 queries/day.
- **Cross-species BLAST:** Asynchronous; implements exponential backoff polling.
