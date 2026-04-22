from __future__ import annotations

import os


class Config:
    SECRET_KEY = os.environ.get("SECRET_KEY", "dev-secret-key")
    UNIPROT_BASE_URL = os.environ.get("UNIPROT_BASE_URL", "https://rest.uniprot.org")
    ENABLE_CURATED_FALLBACK = os.environ.get("ENABLE_CURATED_FALLBACK", "").lower() in {"1", "true", "yes"}
    ANALYSIS_WORKDIR = os.environ.get(
        "ANALYSIS_WORKDIR", os.path.join(os.getcwd(), ".analysis-cache")
    )
    PREDICTOR_TIMEOUT_SECONDS = float(os.environ.get("PREDICTOR_TIMEOUT_SECONDS", "240"))
    ENABLED_PREDICTOR_KEYS = [
        key.strip()
        for key in os.environ.get("ENABLED_PREDICTOR_KEYS", "psipred").split(",")
        if key.strip()
    ]
    INTERPRO_API_BASE_URL = os.environ.get("INTERPRO_API_BASE_URL", "https://www.ebi.ac.uk/interpro/api")
    EBI_PROTEINS_API_BASE_URL = os.environ.get("EBI_PROTEINS_API_BASE_URL", "https://www.ebi.ac.uk/proteins/api")
    ENSEMBL_REST_BASE_URL = os.environ.get("ENSEMBL_REST_BASE_URL", "https://rest.ensembl.org")
    NCBI_DATASETS_BASE_URL = os.environ.get("NCBI_DATASETS_BASE_URL", "https://api.ncbi.nlm.nih.gov/datasets/v2")
    NCBI_EUTILS_BASE_URL = os.environ.get(
        "NCBI_EUTILS_BASE_URL", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    )
    MAVEDB_API_BASE_URL = os.environ.get("MAVEDB_API_BASE_URL", "https://api.mavedb.org/api/v1")
    HUMAN_TAXON_ID = 9606
    REQUEST_TIMEOUT_SECONDS = float(os.environ.get("REQUEST_TIMEOUT_SECONDS", "20"))
    USER_AGENT = os.environ.get("USER_AGENT", "ProteinMutationWorkbench/0.1 academic prototype")
    ANALYSIS_SECTIONS = [
        {
            "key": "resolver",
            "label": "Human protein resolution",
            "status": "implemented",
            "summary": "Resolve the canonical protein sequence for the submitted human target.",
        },
        {
            "key": "mutation",
            "label": "Variant parsing and validation",
            "status": "implemented",
            "summary": "Parse missense substitutions and verify that the wild-type residue matches the resolved sequence.",
        },
        {
            "key": "physicochemical",
            "label": "Wild-type versus mutant metrics",
            "status": "implemented",
            "summary": "Compare baseline sequence composition and ProtParam-style metrics after the substitution.",
        },
        {
            "key": "structure",
            "label": "Secondary structure prediction",
            "status": "implemented",
            "summary": "PSIPRED-derived helix, strand, and coil calls are rendered as intervals, tracks, and residue summaries.",
        },
        {
            "key": "accessibility",
            "label": "RSA and ASA prediction",
            "status": "implemented",
            "summary": "DMVFL-RSA-driven relative and absolute solvent accessibility is rendered with residue-level graphics and tables when the local model repo is present.",
        },
        {
            "key": "conservation",
            "label": "Primate conservation profiling",
            "status": "implemented",
            "summary": "Ensembl primate orthologues are summarized at the submitted residue and rendered as a clean comparative tree plus residue conservation table.",
        },
        {
            "key": "phylogeny",
            "label": "Phylogenetic analysis",
            "status": "implemented",
            "summary": "The report now adds an explicit primate phylogenetic analysis layer with UPGMA clustering, pairwise identity matrix, nearest and farthest orthologues, and Newick export.",
        },
        {
            "key": "genomic_mapping",
            "label": "Chromosome and exon mapping",
            "status": "implemented",
            "summary": "The submitted protein residue is mapped back onto chromosome coordinates, transcript exon structure, and a split exon zoom view using Ensembl plus NCBI-backed fallback when one source is incomplete.",
        },
        {
            "key": "alphafold",
            "label": "AlphaFold2 structural confidence",
            "status": "implemented",
            "summary": "Per-residue pLDDT confidence scores are fetched from the AlphaFold EBI API and rendered as a colour-coded strip with variant-site spotlight. AlphaMissense pathogenicity score is extracted when available.",
        },
        {
            "key": "pdb_xrefs",
            "label": "Experimental PDB structure cross-references",
            "status": "implemented",
            "summary": "UniProt PDB cross-references are resolved to identify experimental structures (X-ray, cryo-EM, NMR) that cover the submitted protein, ranked by resolution and variant-site coverage.",
        },
        {
            "key": "structural_context",
            "label": "Integrated structural context",
            "status": "implemented",
            "summary": "PSIPRED secondary structure, DMVFL-RSA solvent accessibility, and domain overlap are combined into a single variant-site structural verdict. Amphipathic helix analysis is included when the variant falls in a predicted helical segment.",
        },
        {
            "key": "stability",
            "label": "Mutation stability prediction",
            "status": "implemented",
            "summary": "Empirical ΔΔG estimate from BLOSUM62 evolutionary compatibility, hydrophobicity delta, and residue volume change. Classifies the substitution as stabilizing, neutral, or destabilizing.",
        },
    ]
