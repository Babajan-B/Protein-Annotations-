# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Protein Sequence Annotation Workbench** — a Flask web application for analyzing human protein missense mutations. Given a gene symbol, protein name, or UniProt accession + a missense variant (e.g., `R248W`):

1. Resolves the canonical human protein via UniProt REST API
2. Validates the submitted amino-acid substitution
3. Runs a configured pipeline of sequence-based prediction tools
4. Generates a comprehensive report with structural figures, tables, and plots

**Key principle:** The app is a *sequence-first* tool. It does not require 3D structures and works from primary sequence alone. Structural annotations come from sequence-based predictors (PSIPRED, DMVFL-RSA, AlphaFold pLDDT, etc.) and experimental PDB cross-references.

---

## Quick Start

### Development Environment

```bash
# Activate venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Run Flask dev server
./.venv/bin/python run.py
# Server runs on http://127.0.0.1:5000
```

### Run Tests

```bash
# Run all tests
pytest tests/ -v

# Run single test file
pytest tests/test_variant_parser.py -v

# Run specific test function
pytest tests/test_analysis.py::test_some_function -v
```

### Database and Resource Setup

The app uses external public APIs and optional local prediction tools:

- **UniProt REST** — always available, no setup needed
- **EBI Job Dispatcher** — used for remote HMMER, InterProScan, Phobius; no auth required
- **AlphaFold EBI API** — free public API, no setup needed
- **Local predictors** (PSIPRED, DMVFL-RSA, NetSurfP-3.0) — optional; require binary or environment setup

See `true_analysis_resources.md` and `resources/tool-downloads/` for detailed tool installation guides.

---

## Architecture Overview

### Request Flow

```
Browser Request
    ↓
routes/main.py (_build_report_from_request_args)
    ├─→ VariantParser → parse_missense_variant()
    ├─→ HumanProteinResolver → resolve protein from UniProt
    ├─→ SequenceAnalysisService → build report (WT vs mutant metrics)
    ├─→ TrueAnalysisOrchestrator → run enabled predictors (PSIPRED, DMVFL-RSA, etc.)
    ├─→ AnnotationAggregator → (optional) curated fallback from public APIs
    ├─→ AlphaFoldLookup → fetch pLDDT scores from EBI
    ├─→ PdbLookup → find experimental PDB structures via UniProt xrefs
    ├─→ GenomicInsightService → map to chromosome/transcript
    └─→ Jinja2 Template
        ├─ result.html (full report)
        └─ structure_annotations.html (SS + AlphaFold + PDB detail)
```

### Service Architecture

**`app/services/`** — core analysis logic:

| Module | Purpose |
|---|---|
| `variant_parser.py` | Parse and normalize missense variants (e.g., `R248W` → position, wild-type, mutant) |
| `protein_lookup.py` | Resolve human proteins from UniProt (accession, sequence, metadata) |
| `analysis.py` | Build main report: WT/mutant metrics, structural verdict, amphipathic helix detection, secondary structure diagrams |
| `true_analysis.py` | Orchestrate configured predictor wrappers; run and normalize outputs |
| `alphafold_lookup.py` | Fetch AlphaFold pLDDT and AlphaMissense scores from EBI; classify confidence levels |
| `pdb_lookup.py` | Resolve UniProt → PDB cross-references; map to variant position; rank by resolution and coverage |
| `annotation_sources.py` | (Optional curated fallback) Query EBI, Ensembl, NCBI, MaveDB for functional annotations |
| `genomic_insights.py` | Map protein position back to chromosome, transcript exon structure |
| `tool_registry.py` | Manage predictor specs (name, env vars, binary paths, auto-detection) |
| `exporters.py` | Generate JSON/CSV downloads of reports |

**Predictor wrappers** live in `wrappers/`. Each wrapper:

- Accepts placeholders: `{input_fasta}`, `{output_json}`, `{variant_position}`, `{wild_type}`, `{mutant}`, etc.
- Runs the real tool (local binary or remote EBI API call)
- Writes normalized JSON to `{output_json}` following the contract in `true_analysis_resources.md`

Example:
```bash
python run_psipred_single.py \
  --input-fasta /tmp/seq.fa \
  --output-json /tmp/psipred.json \
  --variant-position 248
```

The app reads the JSON and merges results into the report.

### Report Structure

The `analysis.py:SequenceAnalysisService.build_report()` method returns a dictionary with:

- `protein` — ProteinRecord (accession, gene, sequence, length, etc.)
- `variant` — MissenseVariant (position, WT, mutant)
- `mutant_sequence` — full sequence with substitution applied
- `summary_cards` — top-level findings
- `metric_deltas` — WT vs mutant physicochemical shifts (GRAVY, pI, MW, etc.)
- `secondary_structure_panel` — PSIPRED SS prediction; **now includes `diagram_rows`** for PDBsum-style layout
- `accessibility_panel` — DMVFL-RSA RSA/ASA predictions
- `structural_verdict` — integrated SS + RSA verdict card + domain context
- `amphipathic_panel` — helix amphipathicity analysis (if applicable, μH ≥ 0.35)
- `alphafold_panel` — per-residue pLDDT, AlphaMissense score, confidence distribution
- `pdb_panel` — experimental PDB entries ranked by coverage and resolution
- `interval_features` — domain/motif intervals from all sources
- `site_features` — exact-position features (PTMs, binding sites, etc.)
- `evidence_rows` — source annotations
- `architecture_rows` — protein architecture track for main report
- ... and many others

---

## Key Patterns and Conventions

### Predictor Registration

Predictors are registered in `tool_registry.py:PREDICTOR_SPECS`:

```python
PredictorSpec(
    key="psipred",
    label="PSIPRED",
    category="structure",
    summary="Secondary structure prediction...",
    command_env="PSIPRED_COMMAND_TEMPLATE",
    docs_url="https://...",
    resource_hint="Defaults to envs/bio-tools/bin/psipred if present"
)
```

The app will:
1. Check the `PSIPRED_COMMAND_TEMPLATE` environment variable
2. Fall back to a default command template if the variable is unset
3. Skip the predictor if the binary/resource path doesn't exist

Enable/disable predictors via the `ENABLED_PREDICTOR_KEYS` config (comma-separated, default: `"psipred,dmvfl_rsa"`).

### Environment Variables

Predictors are configured via environment variables. The app loads `.env.predictors` if it exists:

```bash
# .env.predictors
PSIPRED_COMMAND_TEMPLATE=python /path/to/run_psipred_single.py --input-fasta {input_fasta} --output-json {output_json} --variant-position {variant_position}
DMVFL_RSA_COMMAND_TEMPLATE=python /path/to/run_dmvfl_rsa.py ...
```

Or set them directly:
```bash
export PSIPRED_COMMAND_TEMPLATE="..."
python run.py
```

### Secondary Structure Diagram Layout

**Recent change:** The `secondary_structure_panel` now includes a `diagram_rows` key (list of dicts) that subdivide all residues into lines (default 100 per line) and group consecutive residues by SS state for PDBsum-style rendering:

```python
# In analysis.py:_secondary_structure_panel()
diagram_rows: list[dict] = _structure_diagram_rows(rows, max_per_line=100)
```

Each `diagram_row`:
```python
{
    "line_number": 1,
    "start": 1,
    "end": 100,
    "segments": [
        {
            "state_code": "C",
            "class_name": "coil",
            "label": "Coil",
            "residues": [...],
            "length": 30
        },
        {
            "state_code": "H",
            "class_name": "helix",
            "label": "Alpha helix",
            "residues": [...],
            "length": 70
        },
        ...
    ]
}
```

The template renders each segment as:
- **Helix**: rounded pill with diagonal stripes (red spring-like shape)
- **Strand**: arrow block with clip-path (blue arrow pointing right)
- **Coil**: thin bottom border line (gray, with residues rendered above)

---

## Testing Strategy

Tests are in `tests/`:

- **Unit tests** — parser, service logic, exporters
- **Integration tests** — end-to-end report building (mocked API calls)
- **Template tests** — Jinja2 rendering

Mock predictor outputs in `tests/fixtures/mock_predictor.py`.

**Run tests before committing:**
```bash
pytest tests/ -v
```

---

## Configuration (config.py)

Key settings:

| Setting | Default | Purpose |
|---|---|---|
| `UNIPROT_BASE_URL` | `https://rest.uniprot.org` | UniProt REST API |
| `ENABLED_PREDICTOR_KEYS` | `psipred,dmvfl_rsa` | Comma-separated list of active predictors |
| `ANALYSIS_WORKDIR` | `.analysis-cache/` | Temp working directory for predictor runs |
| `PREDICTOR_TIMEOUT_SECONDS` | 240 | Timeout for each predictor job |
| `REQUEST_TIMEOUT_SECONDS` | 20 | Timeout for HTTP requests |
| `ENABLE_CURATED_FALLBACK` | False | Enable optional curated annotation sources (slower) |

Override via environment:
```bash
export ENABLED_PREDICTOR_KEYS="psipred,dmvfl_rsa,netsurfp3"
export PREDICTOR_TIMEOUT_SECONDS=300
python run.py
```

---

## Recent Additions (Structural Annotations)

### AlphaFold pLDDT Integration

- **File:** `app/services/alphafold_lookup.py`
- **API:** `https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}`
- **Output:** Per-residue pLDDT fetched from PDB B-factors
- **Rendering:** Colour-coded strip + AlphaMissense badge in report
- **pLDDT levels:**
  - Very high (>90): #0053D6 (dark blue)
  - High (70–90): #65CBF3 (light blue)
  - Low (50–70): #FFDB13 (yellow)
  - Very low (<50): #FF7D45 (orange)

### PDB Cross-References

- **File:** `app/services/pdb_lookup.py`
- **Source:** UniProt PDB cross-references
- **Output:** Ranked list of PDB entries covering (or near) the variant position
- **Ranking:** Experimental structures first (X-ray, cryo-EM, NMR), sorted by resolution; then NMR
- **Table:** Shows method, resolution, chain, residue range, coverage

### Structural Verdict Card

- **File:** `app/services/analysis.py:_structural_verdict()`
- **Combines:** PSIPRED SS class + DMVFL-RSA exposure + domain overlap
- **Output:** Single sentence summary + pills for each signal
- **Location:** Main report after RSA panel, before AlphaFold

### Amphipathic Helix Detection

- **File:** `app/services/analysis.py:_amphipathic_helix_panel()`
- **Method:** Eisenberg hydrophobic moment (μH) for helical segments ≥7 aa
- **Threshold:** μH ≥ 0.35 (weak); ≥0.45 (moderate); ≥0.60 (highly amphipathic)
- **Interpretation:** Suggests membrane-inserting or protein-interaction helix
- **Location:** Rendered in structural verdict card if detected

### NetSurfP-3.0 Wrapper

- **File:** `wrappers/run_netsurfp3_remote.py`
- **Input:** Local NetSurfP-3.0 installation
- **Output:** Per-residue SS, RSA, disorder, and backbone angles
- **Registration:** `tool_registry.py` (auto-detect in `envs/bio-tools/bin/netsurfp3`)

---

## Common Development Tasks

### Add a New Predictor

1. **Create wrapper** in `wrappers/run_<tool>.py`
   - Read `true_analysis_resources.md` for the JSON contract
   - Use `predictor_utils.py` helpers (`interval_feature()`, `site_feature()`, `table()`)
   - Test: `python wrappers/run_<tool>.py --input-fasta test.fa --output-json out.json --variant-position 100`

2. **Register in tool_registry.py**
   ```python
   PredictorSpec(
       key="mytool",
       label="MyTool",
       category="structure",  # or "domains", "topology", "conservation", etc.
       summary="...",
       command_env="MYTOOL_COMMAND_TEMPLATE",
       docs_url="https://...",
       resource_hint="..."
   )
   # Add auto-detection in _default_command_template():
   "mytool": _wrapper_command(python_bin, "run_mytool.py", "...args...")
   ```

3. **Enable in config.py**
   ```python
   ENABLED_PREDICTOR_KEYS = ["psipred", "dmvfl_rsa", "mytool"]
   ```

4. **Test** via `pytest tests/`

### Modify Secondary Structure Rendering

The per-residue diagram is generated in `analysis.py:_structure_diagram_rows()`. It produces `diagram_rows` keyed by `{state_code, class_name, label, residues[], length}`.

Template rendering in `result.html` and `structure_annotations.html` uses:
- CSS classes: `.ss-helix`, `.ss-strand`, `.ss-coil` for shape styling
- Residue cells: `.ss-res` inside each segment, with `.ss-variant` highlight
- Position ticks: every 10 residues

To change layout: adjust `max_per_line` parameter (default 100) in the call to `_structure_diagram_rows()`.

### Generate Report for Testing

```python
from app import create_app
from app.services.variant_parser import parse_missense_variant
from app.services.protein_lookup import HumanProteinResolver

app = create_app()
with app.app_context():
    resolver = HumanProteinResolver(
        base_url="https://rest.uniprot.org",
        timeout_seconds=20
    )
    protein = resolver.resolve("TP53")  # Gene symbol
    variant = parse_missense_variant("R248W")
    
    from app.services.analysis import SequenceAnalysisService
    service = SequenceAnalysisService(app.config["ANALYSIS_SECTIONS"])
    report = service.build_report(protein=protein, variant=variant)
    
    # Inspect report keys
    print(report.keys())
    print(report["secondary_structure_panel"]["diagram_rows"])
```

---

## Known Limitations & TODOs

- **Batch analysis** — not supported; one query at a time
- **Async jobs** — currently synchronous Flask; RQ + Redis deferred for Phase 3
- **3D structures** — not integrated; uses sequence-based predictions only
- **Cache persistence** — `.analysis-cache/` is ephemeral; no database
- **User accounts** — not implemented; academic prototype only

---

## Useful References

- **Product Requirements:** `prd.md`
- **Implementation Plan:** `plan.md`
- **Predictor Contracts:** `true_analysis_resources.md`
- **Wrapper README:** `wrappers/README.md`

---

## Git Workflow

- Main branch: `main`
- Create feature branches: `git checkout -b feature/description`
- Commit messages should be concise: `Add AlphaFold pLDDT panel to reports`
- Run tests before pushing: `pytest tests/ -v`
- No force-pushes to `main`
