# True Analysis Resources

The app now expects real predictor runs to come from configured tool wrappers rather than curated web annotations. Each predictor is enabled by setting a command-template environment variable that writes normalized JSON output for the app.

## Wrapper Contract

Each configured command receives the following placeholders:

- `{input_fasta}`
- `{output_dir}`
- `{output_json}`
- `{sequence_id}`
- `{gene_symbol}`
- `{variant}`
- `{variant_position}`
- `{wild_type}`
- `{mutant}`
- `{organism}`
- `{workdir}`

Use a plain executable-plus-arguments command template. The app tokenizes the command and does not run it through a shell, so shell redirection or shell-specific constructs should live inside your wrapper script rather than inside the environment variable itself.

The command must create `{output_json}` with a payload like:

```json
{
  "summary": "Completed predictor run",
  "interval_features": [
    {
      "label": "DNA-binding domain",
      "category": "domains",
      "start": 100,
      "end": 290
    }
  ],
  "site_features": [
    {
      "label": "Predicted phosphosite",
      "category": "ptm",
      "position": 248
    }
  ],
  "residue_tracks": [
    {
      "key": "disorder",
      "label": "Disorder score",
      "values": [
        {"position": 1, "value": 0.12},
        {"position": 2, "value": 0.18}
      ]
    }
  ],
  "evidence_rows": [],
  "tables": [],
  "notices": []
}
```

## Environment Variables

- `INTERPROSCAN_COMMAND_TEMPLATE`
- `HMMER_COMMAND_TEMPLATE`
- `HHBLITS_COMMAND_TEMPLATE`
- `PSIPRED_COMMAND_TEMPLATE`
- `DMVFL_RSA_COMMAND_TEMPLATE`
- `PHOBIUS_COMMAND_TEMPLATE`
- `DEEPTMHMM_COMMAND_TEMPLATE`
- `SIGNALP6_COMMAND_TEMPLATE`
- `DEEPLOC_COMMAND_TEMPLATE`
- `IUPRED3_COMMAND_TEMPLATE`
- `DISOPRED3_COMMAND_TEMPLATE`
- `MUSITEDEEP_COMMAND_TEMPLATE`
- `EVCOUPLINGS_COMMAND_TEMPLATE`
- `CONSURF_COMMAND_TEMPLATE`

## Practical Notes

- The app does not fake these predictions. If a command template is missing, that predictor is reported as unavailable.
- Some predictors can now auto-configure to bundled local wrappers or official remote services. Environment variables still take precedence when you want to override the defaults.
- For heavy tools that require models or databases, the wrapper command should already know where those resources live.
- Use small shell wrapper scripts if the native tool CLI is complex. Point the corresponding `*_COMMAND_TEMPLATE` variable at that wrapper.
