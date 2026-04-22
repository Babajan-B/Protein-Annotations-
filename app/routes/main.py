from __future__ import annotations

from urllib.parse import urlencode

from flask import Blueprint, Response, current_app, jsonify, make_response, redirect, render_template, request, url_for

from app.services.alphafold_lookup import AlphaFoldLookup
from app.services.analysis import AnalysisError, SequenceAnalysisService
from app.services.annotation_sources import AnnotationAggregator
from app.services.docx_exporter import DocxExportError, report_to_docx
from app.services.pdf_exporter import PdfExportError, report_to_pdf
from app.services.exporters import report_to_csv, serialize_report, structure_annotations_to_csv
from app.services.zip_exporter import report_to_zip
from app.services.genomic_insights import GenomicInsightService
from app.services.pdb_lookup import PdbLookup
from app.services.protein_lookup import HumanProteinResolver, ProteinLookupError
from app.services.open_targets import OpenTargetsLookup
from app.services.protein_atlas import ProteinAtlasLookup
from app.services.stability_prediction import StabilityPredictor
from app.services.tool_registry import build_resource_cards
from app.services.true_analysis import TrueAnalysisError, TrueAnalysisOrchestrator
from app.services.variant_parser import VariantParseError, parse_missense_variant

bp = Blueprint("main", __name__)

_FAVICON_SVG = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">
  <rect width="64" height="64" rx="14" fill="#16324f"/>
  <path d="M18 18h12c9 0 16 6.3 16 14.2S39 46.5 30 46.5H18z" fill="#7fd1b9"/>
  <path d="M30 24c4.8 0 8.5 3.5 8.5 8.2S34.8 40.5 30 40.5h-5.8V24z" fill="#f5fbff"/>
  <circle cx="46" cy="18" r="5" fill="#f2b134"/>
</svg>"""


@bp.route("/", methods=["GET"])
def index():
    form_data = {
        "query": request.args.get("query", "").strip(),
        "variant": request.args.get("variant", "").strip(),
        "position": request.args.get("position", "").strip(),
        "isoform": request.args.get("isoform", "").strip(),
    }
    return render_template(
        "index.html",
        error=request.args.get("error"),
        form_data=form_data,
        analysis_sections=current_app.config["ANALYSIS_SECTIONS"],
        resource_cards=build_resource_cards(set(current_app.config["ENABLED_PREDICTOR_KEYS"])),
    )


@bp.route("/favicon.ico", methods=["GET"])
def favicon():
    return Response(_FAVICON_SVG, mimetype="image/svg+xml")


@bp.route("/analyze", methods=["POST"])
def analyze():
    form_data = _read_form_data(request.form)
    query_params = {key: value for key, value in form_data.items() if value}
    return redirect(url_for("main.report_view") + "?" + urlencode(query_params))


@bp.route("/report", methods=["GET"])
def report_view():
    try:
        report, query_args = _build_report_from_request_args()
        return render_template("result.html", report=report, query_args=query_args)
    except (VariantParseError, ProteinLookupError, AnalysisError, TrueAnalysisError) as exc:
        return redirect(url_for("main.index", error=str(exc), **_read_form_data(request.args)))


@bp.route("/structure-annotations", methods=["GET"])
def structure_annotations_view():
    try:
        report, query_args = _build_report_from_request_args()
        return render_template("structure_annotations.html", report=report, query_args=query_args)
    except (VariantParseError, ProteinLookupError, AnalysisError, TrueAnalysisError) as exc:
        return redirect(url_for("main.index", error=str(exc), **_read_form_data(request.args)))


@bp.route("/structure-annotations.csv", methods=["GET"])
def structure_annotations_csv():
    try:
        report, _ = _build_report_from_request_args()
        response = make_response(structure_annotations_to_csv(report))
        response.headers["Content-Type"] = "text/csv; charset=utf-8"
        response.headers["Content-Disposition"] = "attachment; filename=structure-annotations.csv"
        return response
    except (VariantParseError, ProteinLookupError, AnalysisError, TrueAnalysisError) as exc:
        response = make_response(str(exc), 400)
        response.headers["Content-Type"] = "text/plain; charset=utf-8"
        return response


@bp.route("/report.docx", methods=["GET"])
def report_docx():
    try:
        report, _ = _build_report_from_request_args()
        docx_bytes = report_to_docx(report)
        accession = report["protein"].accession
        variant = report["variant"].short_notation.replace(" ", "_")
        filename = f"report_{accession}_{variant}.docx"
        response = make_response(docx_bytes)
        response.headers["Content-Type"] = (
            "application/vnd.openxmlformats-officedocument.wordprocessingml.document"
        )
        response.headers["Content-Disposition"] = f"attachment; filename={filename}"
        return response
    except DocxExportError as exc:
        response = make_response(str(exc), 503)
        response.headers["Content-Type"] = "text/plain; charset=utf-8"
        return response
    except (VariantParseError, ProteinLookupError, AnalysisError, TrueAnalysisError) as exc:
        response = make_response(str(exc), 400)
        response.headers["Content-Type"] = "text/plain; charset=utf-8"
        return response


@bp.route("/report.pdf", methods=["GET"])
def report_pdf():
    try:
        report, _ = _build_report_from_request_args()
        pdf_bytes = report_to_pdf(report)
        accession = report["protein"].accession
        variant = report["variant"].short_notation.replace(" ", "_")
        filename = f"report_{accession}_{variant}.pdf"
        response = make_response(pdf_bytes)
        response.headers["Content-Type"] = "application/pdf"
        response.headers["Content-Disposition"] = f"attachment; filename={filename}"
        return response
    except PdfExportError as exc:
        response = make_response(str(exc), 503)
        response.headers["Content-Type"] = "text/plain; charset=utf-8"
        return response
    except (VariantParseError, ProteinLookupError, AnalysisError, TrueAnalysisError) as exc:
        response = make_response(str(exc), 400)
        response.headers["Content-Type"] = "text/plain; charset=utf-8"
        return response


@bp.route("/report-bundle.zip", methods=["GET"])
def report_bundle_zip():
    try:
        report, _ = _build_report_from_request_args()
        zip_bytes = report_to_zip(report)
        accession = report["protein"].accession
        variant = report["variant"].short_notation.replace(" ", "_")
        filename = f"report_bundle_{accession}_{variant}.zip"
        response = make_response(zip_bytes)
        response.headers["Content-Type"] = "application/zip"
        response.headers["Content-Disposition"] = f"attachment; filename={filename}"
        return response
    except (VariantParseError, ProteinLookupError, AnalysisError, TrueAnalysisError) as exc:
        response = make_response(str(exc), 400)
        response.headers["Content-Type"] = "text/plain; charset=utf-8"
        return response


@bp.route("/phylogeny.newick", methods=["GET"])
def phylogeny_newick():
    try:
        report, _ = _build_report_from_request_args()
        panel = report.get("phylogenetic_panel") or {}
        newick = panel.get("newick", "").strip()
        if not newick:
            raise AnalysisError("No phylogenetic tree is available for this report.")
        response = make_response(newick + ("\n" if not newick.endswith("\n") else ""))
        response.headers["Content-Type"] = "text/plain; charset=utf-8"
        response.headers["Content-Disposition"] = "attachment; filename=primate-phylogeny.newick"
        return response
    except (VariantParseError, ProteinLookupError, AnalysisError, TrueAnalysisError) as exc:
        response = make_response(str(exc), 400)
        response.headers["Content-Type"] = "text/plain; charset=utf-8"
        return response


@bp.route("/report.json", methods=["GET"])
def report_json():
    try:
        report, _ = _build_report_from_request_args()
        return jsonify(serialize_report(report))
    except (VariantParseError, ProteinLookupError, AnalysisError, TrueAnalysisError) as exc:
        return jsonify({"error": str(exc)}), 400


@bp.route("/report.csv", methods=["GET"])
def report_csv():
    try:
        report, _ = _build_report_from_request_args()
        response = make_response(report_to_csv(report))
        response.headers["Content-Type"] = "text/csv; charset=utf-8"
        response.headers["Content-Disposition"] = "attachment; filename=protein-mutation-report.csv"
        return response
    except (VariantParseError, ProteinLookupError, AnalysisError, TrueAnalysisError) as exc:
        response = make_response(str(exc), 400)
        response.headers["Content-Type"] = "text/plain; charset=utf-8"
        return response


def _build_report_from_request_args():
    form_data = _read_form_data(request.args)
    variant = parse_missense_variant(form_data["variant"])
    submitted_position = _parse_optional_position(form_data["position"])
    if submitted_position and submitted_position != variant.position:
        raise VariantParseError(
            f"Optional position {submitted_position} does not match variant position {variant.position}."
        )

    resolver = HumanProteinResolver(
        base_url=current_app.config["UNIPROT_BASE_URL"],
        timeout_seconds=current_app.config["REQUEST_TIMEOUT_SECONDS"],
    )
    protein = resolver.resolve(query=form_data["query"], isoform=form_data["isoform"] or None)

    orchestrator = TrueAnalysisOrchestrator(
        workdir=current_app.config["ANALYSIS_WORKDIR"],
        timeout_seconds=current_app.config["PREDICTOR_TIMEOUT_SECONDS"],
        enabled_predictors=set(current_app.config["ENABLED_PREDICTOR_KEYS"]),
    )
    annotations = orchestrator.collect(protein=protein, variant=variant)

    if current_app.config["ENABLE_CURATED_FALLBACK"]:
        aggregator = AnnotationAggregator(
            interpro_base_url=current_app.config["INTERPRO_API_BASE_URL"],
            ebi_proteins_base_url=current_app.config["EBI_PROTEINS_API_BASE_URL"],
            ensembl_base_url=current_app.config["ENSEMBL_REST_BASE_URL"],
            ncbi_base_url=current_app.config["NCBI_EUTILS_BASE_URL"],
            mavedb_base_url=current_app.config["MAVEDB_API_BASE_URL"],
            timeout_seconds=current_app.config["REQUEST_TIMEOUT_SECONDS"],
            user_agent=current_app.config["USER_AGENT"],
        )
        curated = aggregator.collect(protein=protein, variant=variant)
        annotations = _merge_annotation_bundles(annotations, curated)

    analysis_service = SequenceAnalysisService(current_app.config["ANALYSIS_SECTIONS"])
    report = analysis_service.build_report(protein=protein, variant=variant, annotations=annotations)
    insight_service = GenomicInsightService(
        ensembl_base_url=current_app.config["ENSEMBL_REST_BASE_URL"],
        ncbi_datasets_base_url=current_app.config["NCBI_DATASETS_BASE_URL"],
        timeout_seconds=current_app.config["REQUEST_TIMEOUT_SECONDS"],
        user_agent=current_app.config["USER_AGENT"],
    )
    report.update(insight_service.collect(protein=protein, variant=variant))

    af_lookup = AlphaFoldLookup(
        timeout_seconds=current_app.config["REQUEST_TIMEOUT_SECONDS"],
        user_agent=current_app.config["USER_AGENT"],
    )
    report.update(
        af_lookup.collect(
            accession=protein.accession,
            variant_position=variant.position,
            wild_type=variant.wild_type,
            mutant=variant.mutant,
            sequence_length=protein.length,
        )
    )

    pdb_lookup = PdbLookup(
        base_url=current_app.config["UNIPROT_BASE_URL"],
        timeout_seconds=current_app.config["REQUEST_TIMEOUT_SECONDS"],
        user_agent=current_app.config["USER_AGENT"],
    )
    report.update(
        pdb_lookup.collect(
            accession=protein.accession,
            variant_position=variant.position,
        )
    )

    stability_predictor = StabilityPredictor()
    report.update(
        stability_predictor.collect(
            wild_type=variant.wild_type,
            mutant=variant.mutant,
            position=variant.position,
        )
    )

    gene_sym = protein.gene_symbol or ""
    atlas_lookup = ProteinAtlasLookup(
        timeout_seconds=current_app.config["REQUEST_TIMEOUT_SECONDS"],
        user_agent=current_app.config["USER_AGENT"],
    )
    report.update(atlas_lookup.collect(gene_symbol=gene_sym))

    ot_lookup = OpenTargetsLookup(
        timeout_seconds=current_app.config["REQUEST_TIMEOUT_SECONDS"],
        user_agent=current_app.config["USER_AGENT"],
    )
    report.update(ot_lookup.collect(gene_symbol=gene_sym))

    return report, form_data


def _read_form_data(source) -> dict[str, str]:
    return {
        "query": source.get("query", "").strip(),
        "variant": source.get("variant", "").strip(),
        "position": source.get("position", "").strip(),
        "isoform": source.get("isoform", "").strip(),
    }


def _parse_optional_position(raw_position: str) -> int | None:
    if not raw_position:
        return None

    if not raw_position.isdigit():
        raise VariantParseError("Optional amino-acid position must be a positive integer.")

    position = int(raw_position)
    if position <= 0:
        raise VariantParseError("Optional amino-acid position must be a positive integer.")
    return position


def _merge_annotation_bundles(primary: dict, secondary: dict) -> dict:
    merged = {
        "interval_features": primary.get("interval_features", []) + secondary.get("interval_features", []),
        "site_features": primary.get("site_features", []) + secondary.get("site_features", []),
        "evidence_rows": primary.get("evidence_rows", []) + secondary.get("evidence_rows", []),
        "residue_tracks": primary.get("residue_tracks", []) + secondary.get("residue_tracks", []),
        "predictor_tables": primary.get("predictor_tables", []) + secondary.get("predictor_tables", []),
        "resource_cards": primary.get("resource_cards", []),
        "architecture_rows": primary.get("architecture_rows", []) or secondary.get("architecture_rows", []),
        "feature_counts": primary.get("feature_counts", []) or secondary.get("feature_counts", []),
        "evidence_plot": primary.get("evidence_plot", []) + secondary.get("evidence_plot", []),
        "source_status": primary.get("source_status", []) + secondary.get("source_status", []),
        "notices": primary.get("notices", []) + secondary.get("notices", []),
    }
    return merged
