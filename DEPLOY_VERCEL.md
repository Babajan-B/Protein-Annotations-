# Vercel Deployment

## What is already configured

- `main.py` exposes the Flask `app` instance in a Vercel-detectable entrypoint.
- `.python-version` pins Python `3.12` for Vercel builds.
- `vercel.json` sets a longer function duration and excludes local-only files from the Python bundle.
- On Vercel, the app uses `/tmp/protein-mutation-workbench` for writable analysis cache data.
- On Vercel, predictor execution is disabled by default unless `ENABLED_PREDICTOR_KEYS` is set explicitly.

## Import the project into Vercel

1. Open Vercel and import the GitHub repository.
2. Keep the project root at the repository root.
3. Let Vercel detect the Python/Flask setup automatically.

## Environment variables

Use `.env.example` as the source for your Vercel project variables.

Recommended minimum:

- `SECRET_KEY`
- `USER_AGENT`
- `EBI_JOB_EMAIL`

Optional:

- `ENABLED_PREDICTOR_KEYS=interproscan,hmmer,phobius`
- `ENABLE_CURATED_FALLBACK=true`

## Automatic deployments

Once the repository is connected in Vercel, pushes to the production branch deploy automatically and pull requests get preview deployments.

## Important runtime note

Local heavy predictors such as PSIPRED and DMVFL-RSA depend on binaries and model assets that are not appropriate for a Vercel serverless deployment. Keep those for local or non-serverless environments.
