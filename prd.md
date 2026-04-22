Here is the focused list of bioinformatics methodologies and the specific sequence-based tools/algorithms we will use for each annotation category. Since we are strictly avoiding 3D structure inputs, all tools listed here predict these features directly from the **1D amino acid sequence**:

### 1. Functional Domains & Motifs
*   **InterProScan:** The master tool for sequence feature mapping. It aggregates multiple databases to find functional domains and sites.
    *   **Pfam / CDD:** For identifying large functional domains.
    *   **PROSITE:** For small catalytic sites, zinc fingers, and highly conserved binding motifs.
    *   **SMART:** For signaling and extracellular domains.

### 2. Secondary Structure Prediction (Sequence-Based)
*   **PSIPRED:** The gold standard for predicting alpha-helices, beta-sheets, and coils directly from the sequence using neural networks.
*   **NetSurfP-3.0:** An advanced deep learning tool that predicts secondary structure directly from the linear sequence.

### 3. Accessible Surface Area (ASA) / Solvent Accessibility
*(Predicting if the mutated residue is buried in the core or exposed to water, purely from sequence)*
*   **NetSurfP-3.0 / SPIDER3:** These tools calculate Relative Solvent Accessibility (RSA) solely from the primary sequence. It will tell us if a mutation flips a residue from "buried" (hydrophobic core) to "exposed" (surface).
*   **SABLE:** Another excellent sequence-based predictor for relative solvent accessibility.

### 4. Evolutionary & Conservation Analysis
*   **HMMER / HHblits:** Used to search sequence databases (like UniRef) to build Hidden Markov Models (HMMs) and Multiple Sequence Alignments (MSAs).
*   **ConSurf (Sequence Mode):** While often used with structures, ConSurf can take a pure sequence, build the MSA, and calculate a strict conservation score (1-9) for every single amino acid position to see if the mutated residue is evolutionarily locked.
*   **EVcouplings (Co-evolution):** Identifies pairs of residues that co-evolve across species, helping predict if mutating one residue will break a necessary interaction with another.

### 5. Transmembrane Topology & Subcellular Localization
*   **DeepTMHMM / Phobius:** Predicts transmembrane helices and the inside/outside topology of the sequence.
*   **SignalP 6.0:** Predicts the presence and exact cleavage site of signal peptides.
*   **DeepLoc 2.0:** Predicts where the protein lives in the cell (nucleus, cytoplasm, mitochondria, etc.) based purely on sequence sorting signals.

### 6. Intrinsic Disorder & Flexibility
*   **IUPred3:** Predicts intrinsically disordered regions (IDRs)—segments of the protein that lack a fixed structure and act as flexible linkers or binding interfaces.
*   **DISOPRED3:** Another highly accurate sequence-based predictor for unstructured regions and protein-binding sites within disordered regions.

### 7. Post-Translational Modifications (PTMs)
*   **MusiteDeep:** A deep-learning tool for predicting sequence-specific PTMs (phosphorylation, glycosylation, ubiquitination, methylation) to see if the mutation destroys a modification site.

### 8. Physicochemical & Thermodynamic Properties
*   **BioPython (ProtParam):** Calculates baseline biophysical shifts caused by the mutation.
    *   **GRAVY Score:** Grand Average of Hydropathicity (shifts in hydrophobicity).
    *   **Isoelectric Point (pI) & Charge:** To see if the mutation shifts the sequence's net charge at physiological pH.
    *   **Aliphatic Index / Instability Index:** General thermal stability metrics based on sequence composition.

Sources
