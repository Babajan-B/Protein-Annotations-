# Manual Downloads Needed

These tools cannot be treated as normal pip dependencies. Some require academic download forms, some require large standalone packages, and some require model or database bundles.

## Academic or Manual Packages

### IUPred3
- Status: academic form download
- Official: https://iupred3.elte.hu/download_new
- What to do:
  - request the academic package
  - unpack it under `resources/tool-downloads/iupred3`
  - create a wrapper script and point `IUPRED3_COMMAND_TEMPLATE` at it

### SignalP 6.0
- Status: academic package / model download
- Official: https://services.healthtech.dtu.dk/services/SignalP-6.0/
- What to do:
  - obtain the official package
  - unpack it under `resources/tool-downloads/signalp6`
  - create a wrapper script and point `SIGNALP6_COMMAND_TEMPLATE` at it

### DeepTMHMM
- Status: standalone download
- Official: https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0
- What to do:
  - obtain the official package
  - unpack it under `resources/tool-downloads/deeptmhmm`
  - create a wrapper script and point `DEEPTMHMM_COMMAND_TEMPLATE` at it

### DeepLoc 2.x
- Status: software download
- Official: https://services.healthtech.dtu.dk/services/DeepLoc-2.0/
- What to do:
  - obtain the official package
  - unpack it under `resources/tool-downloads/deeploc`
  - create a wrapper script and point `DEEPLOC_COMMAND_TEMPLATE` at it

### InterProScan
- Status: large standalone download
- Official: https://interpro-documentation.readthedocs.io/en/latest/download.html
- What to do:
  - install InterProScan under `resources/tool-downloads/interproscan`
  - provision the required member-database data
  - create a wrapper script and point `INTERPROSCAN_COMMAND_TEMPLATE` at it
  - note: the app can also use the official EBI `iprscan5` service by default for a no-install fallback

### DISOPRED3
- Status: source install plus dependencies
- Official: https://github.com/psipred/disopred
- What to do:
  - install it into the local bio env or another controlled location
  - create a wrapper script and point `DISOPRED3_COMMAND_TEMPLATE` at it

### MusiteDeep
- Status: source install with model/runtime requirements
- Official: https://github.com/duolinwang/MusiteDeep
- What to do:
  - install it into an isolated env
  - point `MUSITEDEEP_COMMAND_TEMPLATE` at the wrapper

### ConSurf
- Status: no standard local one-command install path
- Official: https://consurf.tau.ac.il/
- What to do:
  - use an internal/local conservation workflow or an academic installation
  - point `CONSURF_COMMAND_TEMPLATE` at the wrapper

## Openly Installable Core

These are the tools I can usually install directly:
- HMMER
- HH-suite
- PSIPRED
- BLAST
- Clustal Omega
- MAFFT
- EVcouplings
- Phobius / InterProScan / HMMER remote services via EBI Job Dispatcher

Use `scripts/bootstrap_bio_tools.sh` for that first wave.
