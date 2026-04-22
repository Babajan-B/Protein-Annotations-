# Wrapper Scripts

Bundled wrappers now exist for:

- `run_psipred_single.py`
- `run_hmmer_phmmer_remote.py`
- `run_phobius_remote.py`
- `run_iprscan5_remote.py`
- `run_dmvfl_rsa.py`

Each predictor wrapper should:

1. Accept the placeholders passed through the corresponding `*_COMMAND_TEMPLATE`
2. Run the real tool
3. Convert the tool output into the normalized JSON contract in `true_analysis_resources.md`

The app already knows how to consume that JSON contract. What is still missing for each tool is the tool-specific wrapper implementation once the real binary or package is present.
