# Protein Sequence Annotation Workbench Plan

## Goal

Build a Flask-based academic web application for human protein mutation analysis from a gene symbol, protein name, or UniProt accession plus a required protein missense variant. The app should resolve the canonical human protein, validate the submitted amino-acid substitution, aggregate sequence-based annotations, and generate figures, plots, and tables for interpretation.

## Product Scope

### V1 Scope

- Human proteins only.
- Academic use only.
- Single amino-acid substitutions only, for example `R248W` or `p.Arg248Trp`.
- Inputs:
  - required `query`: gene symbol, protein name, or UniProt accession
  - required `variant`: missense substitution
  - optional `position`: must match the position encoded in the variant if supplied
  - optional `isoform`: explicit UniProt isoform or transcript override
- Outputs:
  - resolved protein metadata
  - variant validation summary
  - wild-type versus mutant physicochemical comparison
  - protein architecture figure
  - residue-centric context panel
  - annotation tables
  - downloadable JSON/CSV

### Deferred

- Indels, nonsense, splice, frameshift, and synonymous variants
- Multi-variant haplotypes
- Batch uploads
- User accounts
- Commercial deployment hardening

## Analysis Modules

### V1 Implementation Order

1. Resolver and mutation parser
2. Wild-type versus mutant sequence metrics
3. Protein overview figure and residue context panel
4. True-analysis predictor orchestration
5. Optional curated fallback and background jobs

### Canonical Sequence Source

- `UniProt REST`: sequence resolution, accession metadata, canonical protein information

### Planned True Predictors

- `InterProScan` for domains and motifs
- `HMMER / HHblits` for homolog search and alignment seeding
- `DMVFL-RSA` for RSA and ASA accessibility profiles
- `PSIPRED` for secondary structure
- `DeepTMHMM / SignalP 6.0 / DeepLoc 2.x` for topology and localization
- `IUPred3 / DISOPRED3` for disorder
- `MusiteDeep` for PTM prediction
- `EVcouplings / ConSurf` for conservation and coupling context

## Figures, Plots, and Tables

### Figures

- Protein architecture track with variant marker
- Residue neighborhood view around the mutated site
- Domain and motif overlay

### Plots

- Wild-type versus mutant physicochemical delta chart
- Conservation score track with highlighted residue
- Disorder score track
- Solvent accessibility track
- PTM and evidence lollipop track

### Tables

- Protein summary
- Variant summary
- Residue context summary
- Feature interval table
- Evidence provenance table
- Predictor output table

## Architecture

### App Layout

```text
app/
  __init__.py
  routes/
    main.py
  services/
    analysis.py
    protein_lookup.py
    variant_parser.py
  templates/
  static/
config.py
run.py
tests/
```

### Service Boundaries

- `variant_parser.py`
  - parse and normalize missense variants
  - validate syntax and residue alphabet
- `protein_lookup.py`
  - resolve human protein metadata and canonical sequence from UniProt
- `analysis.py`
  - validate variant against the resolved wild-type sequence
  - build the mutant sequence
  - compute baseline deltas
  - emit normalized payloads that the UI can render
- `true_analysis.py`
  - orchestrate configured predictor wrappers
  - merge normalized true-analysis outputs into report sections

### Data Model Direction

Normalize all outputs into a single report structure with:

- `protein_metadata`
- `variant_summary`
- `interval_features`
- `site_features`
- `residue_scores`
- `metric_deltas`
- `source_provenance`

## Execution Model

### Immediate

- Synchronous Flask requests while the app does lookup, validation, and any configured predictor calls

### Next

- `RQ + Redis` background jobs once the heavier predictor wrappers and databases are connected

## Delivery Phases

### Phase 1

- Scaffold Flask application
- Add search form and validation
- Add human protein resolver against UniProt
- Add missense parser and wild-type residue validation
- Add first result page with summary cards, tables, and placeholder plots

### Phase 2

- Add predictor registry and wrapper contracts
- Add resource audit and wrapper inventory UI
- Add first full protein architecture figure
- Add CSV/JSON export

### Phase 3

- Run configured true-analysis predictors and merge normalized outputs
- Add async jobs and persistent cache

### Phase 4

- Add optional curated fallback only where it helps interpretation and does not replace true predictors

## Risks

- Most advanced predictors are not pip-installable libraries; they require binaries, models, academic downloads, or sequence databases.
- Protein names are ambiguous; gene symbol plus explicit isoform gives cleaner resolution.
- Mutation claims should be clearly labeled as true predictor output versus missing module or optional fallback.

## Current Implementation Slice

This repository will start with:

- `plan.md`
- Flask scaffold
- variant parsing and validation
- human UniProt lookup
- baseline wild-type versus mutant metrics
- true-analysis tool contracts and wrapper inventory
- automatic default wrappers for:
  - `iprscan5` via the official EBI Job Dispatcher service
  - local `DMVFL-RSA` for RSA and ASA prediction
  - `HMMER phmmer` via the official EBI Job Dispatcher service
  - local `PSIPRED` single-sequence mode from the bundled bio-tools environment
  - `Phobius` via the official EBI Job Dispatcher service
- real architecture, conservation, topology, and structure outputs without requiring manual environment exports
