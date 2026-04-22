from __future__ import annotations

import io
from datetime import date
from typing import Any

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import cm, mm
    from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
        HRFlowable, PageBreak, KeepTogether
    )
    from reportlab.platypus.flowables import HRFlowable
    _RL_AVAILABLE = True
except ImportError:
    _RL_AVAILABLE = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from reportlab.platypus import Image as RLImage
    _MPL_AVAILABLE = True
except ImportError:
    _MPL_AVAILABLE = False


class PdfExportError(RuntimeError):
    pass


# ── Colour palette ────────────────────────────────────────────────────────
_GREEN      = colors.HexColor("#0b6b57")
_GREEN_LIGHT = colors.HexColor("#e8f5ef")
_BLUE       = colors.HexColor("#1e3a5f")
_BLUE_LIGHT = colors.HexColor("#eaeff8")
_AMBER      = colors.HexColor("#b06520")
_AMBER_LIGHT = colors.HexColor("#fef6ea")
_GRAY       = colors.HexColor("#555555")
_LIGHT_GRAY = colors.HexColor("#f4f7f4")
_RED        = colors.HexColor("#c0392b")
_WHITE      = colors.white
_BLACK      = colors.black
_HEADER_BG  = colors.HexColor("#16324f")


def report_to_pdf(report: dict[str, Any]) -> bytes:
    """Generate a professional PDF from a report dict. Returns raw bytes."""
    if not _RL_AVAILABLE:
        raise PdfExportError("reportlab is not installed. Run: pip install reportlab")

    protein = report["protein"]
    variant = report["variant"]
    gene    = protein.gene_symbol or protein.entry_name or protein.accession

    buf = io.BytesIO()
    doc = SimpleDocTemplate(
        buf,
        pagesize=A4,
        leftMargin=2*cm, rightMargin=2*cm,
        topMargin=2.5*cm, bottomMargin=2*cm,
        title=f"{gene} {variant.short_notation} — Mutation Analysis Report",
        author="Protein Mutation Workbench",
    )

    styles = _build_styles()
    story  = []

    # ── Cover ────────────────────────────────────────────────────────────
    story += _cover(protein, variant, gene, styles)
    story.append(PageBreak())

    # ── Table of contents hint ───────────────────────────────────────────
    toc_items = [
        "1. Executive Summary",
        "2. Physicochemical Metrics",
        "3. Residue Context",
        "4. Mutation Stability Prediction",
        "5. Secondary Structure (PSIPRED)",
        "6. Solvent Accessibility (DMVFL-RSA)",
        "7. Integrated Structural Verdict",
        "8. AlphaFold2 Structural Confidence",
        "9. Experimental PDB Structures",
        "10. Primate Conservation",
        "11. Genomic Location",
        "12. Human Protein Atlas",
        "13. Open Targets Platform",
        "14. Feature Annotations",
        "15. Evidence & Data Sources",
    ]
    story.append(_h1("Contents", styles))
    for item in toc_items:
        story.append(Paragraph(item, styles["toc"]))
    story.append(Spacer(1, 6*mm))
    story.append(PageBreak())

    # ── 1. Executive Summary ─────────────────────────────────────────────
    story.append(_h1("1. Executive Summary", styles))
    story += _kv_table([
        ("Accession",       protein.accession),
        ("Gene symbol",     gene),
        ("Protein name",    protein.protein_name or "—"),
        ("Organism",        protein.organism_name or "—"),
        ("Sequence length", f"{protein.length} aa"),
        ("Variant submitted",   variant.original_text),
        ("Variant normalized",  variant.short_notation),
        ("HGVS protein",        variant.hgvs_protein or "—"),
        ("Residue change",      f"{variant.wild_type} → {variant.mutant} at position {variant.position}"),
    ], styles)

    sc = report.get("summary_cards") or []
    if sc:
        story.append(_h2("Key Findings", styles))
        story += _kv_table([(c["label"], str(c["value"])) for c in sc], styles)

    sp = report.get("stability_panel")
    if sp:
        story.append(_callout(f"Stability verdict: {sp.get('summary','')}", styles, color=_GREEN_LIGHT))

    # ── 2. Physicochemical Metrics ────────────────────────────────────────
    story.append(_h1("2. Physicochemical Metrics (WT vs Mutant)", styles))
    md = report.get("metric_deltas") or []
    if md:
        rows = [["Metric", "Wild-type", "Mutant", "Δ (Mutant − WT)"]]
        for m in md:
            wt  = f"{m.wild_type:.4f}" if m.wild_type  is not None else "N/A"
            mt  = f"{m.mutant:.4f}"    if m.mutant      is not None else "N/A"
            dlt = f"{m.delta:+.4f}"   if m.delta       is not None else "N/A"
            rows.append([m.label, wt, mt, dlt])
        story.append(_data_table(rows, styles, col_widths=[7*cm, 3.5*cm, 3.5*cm, 3.5*cm]))

        if _MPL_AVAILABLE:
            img = _metric_chart(md)
            if img:
                story.append(RLImage(img, width=14*cm, height=7*cm))

    # ── 3. Residue Context ────────────────────────────────────────────────
    story.append(_h1("3. Residue Context", styles))
    ctx = report.get("context") or {}
    if ctx:
        story.append(Paragraph(
            f"Sequence window: residues <b>{ctx.get('window_start','?')}</b> – "
            f"<b>{ctx.get('window_end','?')}</b>",
            styles["body"]))
        story.append(Spacer(1, 3*mm))
        story.append(Paragraph(
            f"<b>WT:</b>  <font name='Courier'>{ctx.get('wild_window','')}</font>",
            styles["body"]))
        story.append(Paragraph(
            f"<b>Mut:</b> <font name='Courier'>{ctx.get('mutant_window','')}</font>",
            styles["body"]))
        wtp = ctx.get("wild_type_properties") or {}
        mtp = ctx.get("mutant_properties") or {}
        if wtp or mtp:
            story.append(_h2("Residue Properties at Variant Site", styles))
            story += _kv_table([
                ("WT amino acid",     variant.wild_type),
                ("WT hydropathy",     str(wtp.get("hydropathy","—"))),
                ("WT charge",         wtp.get("charge_class","—")),
                ("WT polarity",       wtp.get("polarity_class","—")),
                ("WT aromatic",       "Yes" if wtp.get("aromatic") else "No"),
                ("Mutant amino acid", variant.mutant),
                ("Mutant hydropathy", str(mtp.get("hydropathy","—"))),
                ("Mutant charge",     mtp.get("charge_class","—")),
                ("Mutant polarity",   mtp.get("polarity_class","—")),
                ("Mutant aromatic",   "Yes" if mtp.get("aromatic") else "No"),
            ], styles)

    # ── 4. Stability Prediction ───────────────────────────────────────────
    story.append(_h1("4. Mutation Stability Prediction", styles))
    if sp:
        ddg_color = _GREEN_LIGHT if sp.get("ddg_class") == "stabilizing" else (
            _AMBER_LIGHT if sp.get("ddg_class") == "neutral" else colors.HexColor("#fdecea"))
        story.append(_callout(
            f"ΔΔG estimate: {sp['ddg_estimate']:+.2f} kcal/mol  —  {sp.get('ddg_label','')}", styles,
            color=ddg_color, bold=True))
        story += _kv_table([
            ("BLOSUM62 score",           str(sp.get("blosum62_score","—"))),
            ("BLOSUM62 interpretation",  sp.get("blosum62_interpretation","—")),
            ("Hydrophobicity Δ (KD)",    f"{sp.get('hydro_delta',0):+.2f}"),
            ("Volume Δ (Å³)",            f"{sp.get('vol_delta',0):+.1f}"),
            ("SIFT-proxy tolerance",     f"{sp.get('sift_proxy','—')} ({sp.get('sift_tolerance','—')})"),
        ], styles)
        if sp.get("property_changes"):
            story.append(_h2("Property Changes", styles))
            rows = [["Property", "From", "To", "Impact"]]
            for c in sp["property_changes"]:
                rows.append([c["property"], c["from"], c["to"], c["impact"]])
            story.append(_data_table(rows, styles, col_widths=[4.5*cm, 4*cm, 4*cm, 5*cm]))
        story.append(Paragraph(sp.get("summary",""), styles["body"]))
    else:
        story.append(Paragraph("Stability prediction data not available.", styles["muted"]))

    # ── 5. Secondary Structure ────────────────────────────────────────────
    story.append(_h1("5. Secondary Structure Prediction (PSIPRED)", styles))
    ss = report.get("secondary_structure_panel")
    if ss:
        vr = ss.get("variant_row") or {}
        story.append(_callout(
            f"Variant {variant.wild_type}{variant.position}{variant.mutant}: "
            f"{vr.get('state','—')} with {vr.get('confidence_pct','—')}% confidence",
            styles, color=_BLUE_LIGHT))
        cards = ss.get("state_cards") or []
        if cards:
            story.append(_h2("State Distribution", styles))
            rows = [["State", "Residue count", "Percent (%)"]]
            for c in cards:
                rows.append([c["label"], str(c["count"]), f"{c['percent']}%"])
            story.append(_data_table(rows, styles, col_widths=[6*cm, 5*cm, 5*cm]))
            if _MPL_AVAILABLE:
                img = _ss_pie(cards)
                if img:
                    story.append(RLImage(img, width=8*cm, height=8*cm))
        rows_ss = ss.get("rows") or []
        if rows_ss:
            story.append(_h2("Per-residue Predictions (first 60)", styles))
            rows = [["Pos", "AA", "State", "Confidence (%)"]]
            for r in rows_ss[:60]:
                rows.append([str(r.get("position","—")), r.get("aa","—"),
                             r.get("state","—"), str(r.get("confidence_pct","—"))])
            story.append(_data_table(rows, styles, col_widths=[3*cm, 3*cm, 6*cm, 5.5*cm]))
    else:
        story.append(Paragraph(
            "PSIPRED predictor not configured. Install PSIPRED and set PSIPRED_COMMAND_TEMPLATE in .env.predictors.",
            styles["muted"]))

    # ── 6. Solvent Accessibility ──────────────────────────────────────────
    story.append(_h1("6. Solvent Accessibility (DMVFL-RSA)", styles))
    acc = report.get("accessibility_panel")
    if acc:
        vr = acc.get("variant_row") or {}
        story += _kv_table([
            ("Mean RSA (protein-wide)",     str(acc.get("mean_rsa","—"))),
            ("Mean ASA (Å²)",               str(acc.get("mean_asa","—"))),
            ("Exposed residues (RSA > 0.2)",str(acc.get("exposed_count","—"))),
            ("Buried residues (RSA ≤ 0.2)", str(acc.get("buried_count","—"))),
            (f"Variant RSA (pos. {variant.position})", str(vr.get("RSA","—"))),
            ("Variant ASA (Å²)",            str(vr.get("ASA","—"))),
            ("Variant exposure class",      str(vr.get("Exposure","—"))),
        ], styles)
    else:
        story.append(Paragraph(
            "DMVFL-RSA predictor not configured. Install it and set DMVFL_RSA_COMMAND_TEMPLATE.",
            styles["muted"]))

    # ── 7. Structural Verdict ─────────────────────────────────────────────
    story.append(_h1("7. Integrated Structural Verdict", styles))
    sv = report.get("structural_verdict")
    if sv:
        story.append(_callout(sv.get("verdict_sentence",""), styles, color=_GREEN_LIGHT, bold=True))
        if sv.get("domain_context"):
            story.append(Paragraph(f"Domain context: {sv['domain_context']}", styles["body"]))
        if sv.get("parts"):
            story.append(Paragraph("Evidence pillars: " + "  ·  ".join(sv["parts"]), styles["body"]))
    else:
        story.append(Paragraph("No structural verdict available (requires PSIPRED + DMVFL-RSA).", styles["muted"]))

    ap = report.get("amphipathic_panel")
    if ap:
        story.append(_h2("Amphipathic Helix Detected", styles))
        story += _kv_table([
            ("Segment",    f"{ap.get('segment_start','—')}–{ap.get('segment_end','—')}"),
            ("μH (moment)",str(ap.get("mu_h","—"))),
            ("Class",      ap.get("class_label","—")),
            ("Summary",    ap.get("summary","—")),
        ], styles)

    # ── 8. AlphaFold2 ─────────────────────────────────────────────────────
    story.append(_h1("8. AlphaFold2 Structural Confidence", styles))
    af = report.get("alphafold_panel")
    if af:
        story.append(Paragraph(
            f"Entry: <b>{af.get('entry_id','—')}</b>  ·  "
            f"Model date: {af.get('model_date','—')}  ·  "
            f"Coverage: residues {af.get('uniprot_start','—')}–{af.get('uniprot_end','—')}",
            styles["body"]))
        if af.get("plddt_available"):
            vr = af.get("variant_row") or {}
            story += _kv_table([
                ("Mean pLDDT",     f"{af.get('mean_plddt','—')} — {af.get('mean_plddt_label','—')}"),
                (f"Variant pLDDT (pos. {variant.position})",
                 f"{vr.get('plddt','—')} — {vr.get('plddt_label','—')}"),
            ], styles)
            conf = af.get("confidence_cards") or []
            if conf:
                story.append(_h2("pLDDT Confidence Distribution", styles))
                rows = [["Band", "Count", "Percent (%)", "Colour"]]
                for c in conf:
                    rows.append([c["label"], str(c["count"]), f"{c['percent']}%", c.get("color","—")])
                story.append(_data_table(rows, styles, col_widths=[6*cm, 3.5*cm, 3.5*cm, 4.5*cm]))
        am = af.get("alphamissense")
        if am:
            story.append(_h2("AlphaMissense Pathogenicity Score", styles))
            story += _kv_table([
                ("Score",          str(am.get("score","—"))),
                ("Class",          am.get("label","—")),
                ("Interpretation", am.get("interpretation","—")),
            ], styles)
    else:
        story.append(Paragraph("No AlphaFold entry found for this protein.", styles["muted"]))

    # ── 9. Experimental PDB Structures ───────────────────────────────────
    story.append(_h1("9. Experimental PDB Structures", styles))
    pdb = report.get("pdb_panel")
    if pdb:
        story.append(Paragraph(
            f"<b>{pdb.get('total_count',0)}</b> PDB entries found for {protein.accession}. "
            f"<b>{pdb.get('variant_covered_count',0)}</b> cover position {variant.position}.",
            styles["body"]))
        best = pdb.get("best_entry")
        if best:
            story.append(_h2("Best Structure at Variant Site", styles))
            story += _kv_table([
                ("PDB ID",       best.get("pdb_id","—")),
                ("Method",       best.get("method","—")),
                ("Resolution",   str(best.get("resolution","—"))),
                ("Chain",        best.get("chain","—")),
                ("Residue range",f"{best.get('residue_start','—')}–{best.get('residue_end','—')}"),
            ], styles)
        entries = (pdb.get("entries") or [])[:20]
        if entries:
            story.append(_h2("All PDB Entries (top 20)", styles))
            rows = [["PDB ID", "Method", "Resolution", "Chain", "Range", "Covers"]]
            for e in entries:
                rows.append([e.get("pdb_id","—"), e.get("method","—"),
                             str(e.get("resolution","—")), e.get("chain","—"),
                             f"{e.get('residue_start','—')}–{e.get('residue_end','—')}",
                             "Yes" if e.get("covers_variant") else "No"])
            story.append(_data_table(rows, styles,
                col_widths=[2.5*cm, 3.5*cm, 3*cm, 2*cm, 3.5*cm, 3*cm]))
    else:
        story.append(Paragraph("No experimental PDB cross-references found.", styles["muted"]))

    # ── 10. Primate Conservation ──────────────────────────────────────────
    story.append(_h1("10. Primate Conservation", styles))
    cons = report.get("conservation_panel")
    if cons:
        story += _kv_table([
            ("Conservation at position",    f"{cons.get('conservation_pct','—')}%"),
            ("Human residue",               cons.get("human_residue","—")),
            ("Dominant residue in primates",cons.get("dominant_residue","—")),
            ("Conserved orthologues",       f"{cons.get('conserved_count','—')}/{cons.get('comparable_count','—')}"),
            ("Orthologues matching mutant", str(cons.get("mutant_match_count","—"))),
        ], styles)
        rows_c = cons.get("species_rows") or []
        if rows_c:
            story.append(_h2("Species Table", styles))
            pos = cons.get("variant_position","?")
            rows = [[f"Residue at {pos}", "Species", "Match class"]]
            for r in rows_c:
                rows.append([r["site_residue"], r["species_label"], r["match_class"]])
            story.append(_data_table(rows, styles, col_widths=[4*cm, 8*cm, 5.5*cm]))
    else:
        story.append(Paragraph("Conservation data not available.", styles["muted"]))

    # ── 11. Genomic Location ──────────────────────────────────────────────
    story.append(_h1("11. Genomic Location", styles))
    chr_panel = report.get("chromosome_location_panel")
    if chr_panel:
        story += _kv_table([
            ("Chromosome",    f"chr{chr_panel.get('chromosome','—')}"),
            ("Gene",          chr_panel.get("gene_symbol","—")),
            ("Gene span",     chr_panel.get("gene_span_label","—")),
            ("Strand",        chr_panel.get("strand_label","—")),
            ("Transcript",    chr_panel.get("transcript_id","—")),
            ("Variant span",  chr_panel.get("variant_span_label","—")),
        ], styles)
    else:
        story.append(Paragraph("Genomic location data not available.", styles["muted"]))

    # ── 12. Human Protein Atlas ───────────────────────────────────────────
    story.append(_h1("12. Human Protein Atlas", styles))
    pa = report.get("protein_atlas_panel")
    if pa:
        story.append(Paragraph(
            f"Gene: <b>{pa.get('gene_symbol','—')}</b>  ·  "
            f"Ensembl: {pa.get('ensembl_id','—')}  ·  "
            f"Evidence: {pa.get('evidence','—')}", styles["body"]))
        if pa.get("gene_description"):
            story.append(Paragraph(pa["gene_description"], styles["body"]))
        for label, key in [
            ("Protein Class",      "protein_class"),
            ("Biological Process", "biological_process"),
            ("Molecular Function", "molecular_function"),
        ]:
            vals = pa.get(key) or []
            if vals:
                story.append(_h2(label, styles))
                story.append(Paragraph(" · ".join(vals), styles["body"]))
        locs = pa.get("subcellular_locations") or []
        if locs:
            story.append(_h2("Subcellular Localisation", styles))
            rows = [["Location", "Reliability"]]
            for l in locs:
                rows.append([l["location"], l.get("reliability","—")])
            story.append(_data_table(rows, styles, col_widths=[9*cm, 8.5*cm]))
        diseases = pa.get("disease_associations") or []
        if diseases:
            story.append(_h2("Disease Involvement", styles))
            rows = [["Disease", "Source"]]
            for d in diseases:
                rows.append([d["disease"], d.get("source","HPA")])
            story.append(_data_table(rows, styles, col_widths=[13*cm, 4.5*cm]))
        rna = pa.get("rna_fields") or []
        if rna:
            story.append(_h2("RNA Expression Summary", styles))
            story += _kv_table([(r["label"], r["value"]) for r in rna], styles)
        if pa.get("interactions"):
            story.append(Paragraph(f"Known protein interactions: {pa['interactions']}", styles["body"]))
    else:
        story.append(Paragraph("No Human Protein Atlas data retrieved for this gene.", styles["muted"]))

    # ── 13. Open Targets Platform ─────────────────────────────────────────
    story.append(_h1("13. Open Targets Platform", styles))
    ot = report.get("open_targets_panel")
    if ot:
        story.append(Paragraph(
            f"Ensembl: <b>{ot.get('ensembl_id','—')}</b>  ·  "
            f"Biotype: {ot.get('biotype','—')}  ·  "
            f"{ot.get('total_disease_count',0)} disease associations  ·  "
            f"{ot.get('total_drug_count',0)} known drugs",
            styles["body"]))
        funcs = ot.get("function_descriptions") or []
        if funcs:
            story.append(_h2("Function Descriptions", styles))
            for f in funcs:
                story.append(Paragraph(f"• {f}", styles["body"]))
        disease_rows = ot.get("disease_rows") or []
        if disease_rows:
            story.append(_h2(
                f"Disease Associations (top {len(disease_rows)} of {ot.get('total_disease_count',0)})",
                styles))
            rows = [["Disease", "Score", "Therapeutic Areas"]]
            for r in disease_rows:
                rows.append([r["disease_name"], str(r["score"]), r.get("therapeutic_areas","—")])
            story.append(_data_table(rows, styles, col_widths=[7*cm, 3*cm, 7.5*cm]))
        drug_rows = ot.get("drug_rows") or []
        if drug_rows:
            story.append(_h2(f"Known Drugs ({ot.get('total_drug_count',0)} total)", styles))
            rows = [["Drug", "Type", "Phase", "Mechanism", "Indication"]]
            for r in drug_rows:
                rows.append([r["drug_name"], r.get("drug_type","—"), r.get("phase_label","—"),
                             (r.get("mechanism") or "—")[:40], (r.get("indication") or "—")[:40]])
            story.append(_data_table(rows, styles,
                col_widths=[3.5*cm, 3*cm, 2.5*cm, 5*cm, 3.5*cm]))
        pathways = ot.get("pathways") or []
        if pathways:
            story.append(_h2("Pathway Memberships", styles))
            rows = [["Pathway", "Top-level Term"]]
            for pw in pathways:
                rows.append([pw["name"], pw.get("top_level","—")])
            story.append(_data_table(rows, styles, col_widths=[10*cm, 7.5*cm]))
    else:
        story.append(Paragraph("No Open Targets data retrieved for this gene.", styles["muted"]))

    # ── 14. Feature Annotations ───────────────────────────────────────────
    story.append(_h1("14. Feature Annotations", styles))
    ifeats = report.get("interval_features") or []
    sfeats = report.get("site_features") or []
    if ifeats:
        story.append(_h2("Interval Features (Domains / Regions)", styles))
        rows = [["Source", "Type", "Label", "Coordinates", "Variant Overlap"]]
        for f in ifeats:
            rows.append([f["source"], f["type"], f["label"],
                        f"{f['start']}-{f['end']}", "Yes" if f.get("variant_hit") else "No"])
        story.append(_data_table(rows, styles, col_widths=[3*cm, 3*cm, 5*cm, 3.5*cm, 3*cm]))
    if sfeats:
        story.append(_h2("Site Features (PTMs / Binding Sites)", styles))
        rows = [["Source", "Type", "Label", "Position", "Evidence"]]
        for f in sfeats:
            rows.append([f["source"], f["type"], f["label"],
                        str(f["position"]), (f.get("evidence") or "—")[:50]])
        story.append(_data_table(rows, styles, col_widths=[3*cm, 3*cm, 5*cm, 2.5*cm, 4*cm]))
    if not ifeats and not sfeats:
        story.append(Paragraph("No feature annotation data available.", styles["muted"]))

    # ── 15. Evidence & Sources ────────────────────────────────────────────
    story.append(_h1("15. Evidence & Data Sources", styles))
    evidence = report.get("evidence_rows") or []
    if evidence:
        story.append(_h2("Evidence Annotations", styles))
        rows = [["Source", "Label", "Position", "Summary"]]
        for e in evidence:
            rows.append([e["source"], e["label"], str(e["position"]),
                        (e.get("summary") or "—")[:60]])
        story.append(_data_table(rows, styles, col_widths=[3.5*cm, 4*cm, 2.5*cm, 7.5*cm]))

    prov = report.get("source_provenance") or []
    if prov:
        story.append(_h2("Data Sources", styles))
        rows = [["Source", "Purpose"]]
        for p in prov:
            rows.append([p["source"], (p.get("purpose") or "—")[:80]])
        story.append(_data_table(rows, styles, col_widths=[5*cm, 12.5*cm]))

    # ── Footer note ───────────────────────────────────────────────────────
    story.append(Spacer(1, 8*mm))
    story.append(HRFlowable(width="100%", thickness=0.5, color=_GRAY))
    story.append(Paragraph(
        f"Generated by Protein Mutation Workbench · {date.today().isoformat()} · "
        f"Academic use only — not a clinical diagnostic tool.",
        styles["footer"]))

    doc.build(story, onFirstPage=_page_footer, onLaterPages=_page_footer)
    return buf.getvalue()


# ── Cover page ────────────────────────────────────────────────────────────

def _cover(protein, variant, gene: str, styles) -> list:
    elems = []
    elems.append(Spacer(1, 3*cm))

    # Title block
    title_tbl = Table([[
        Paragraph(
            f"<font color='#0b6b57'><b>Protein Mutation</b></font><br/>"
            f"<font color='#1e3a5f'><b>Analysis Report</b></font>",
            styles["cover_title"])
    ]], colWidths=[17.5*cm])
    title_tbl.setStyle(TableStyle([
        ("ALIGN", (0,0), (-1,-1), "CENTER"),
        ("BOTTOMPADDING", (0,0), (-1,-1), 12),
    ]))
    elems.append(title_tbl)
    elems.append(Spacer(1, 1*cm))

    # Variant pill
    variant_tbl = Table([[
        Paragraph(f"<font size='20'><b>{gene}</b></font>  "
                  f"<font size='18' color='#c65d44'>{variant.short_notation}</font>",
                  styles["cover_variant"])
    ]], colWidths=[17.5*cm])
    variant_tbl.setStyle(TableStyle([
        ("ALIGN", (0,0), (-1,-1), "CENTER"),
        ("BACKGROUND", (0,0), (-1,-1), _GREEN_LIGHT),
        ("ROUNDEDCORNERS", [10]),
        ("TOPPADDING", (0,0), (-1,-1), 14),
        ("BOTTOMPADDING", (0,0), (-1,-1), 14),
        ("LEFTPADDING", (0,0), (-1,-1), 24),
        ("RIGHTPADDING", (0,0), (-1,-1), 24),
    ]))
    elems.append(variant_tbl)
    elems.append(Spacer(1, 0.8*cm))

    # Protein metadata
    elems.append(Paragraph(protein.protein_name or "", styles["cover_meta"]))
    elems.append(Paragraph(
        f"{protein.accession}  ·  {protein.length} aa  ·  {protein.organism_name or ''}",
        styles["cover_meta"]))
    elems.append(Spacer(1, 4*cm))
    elems.append(Paragraph(f"Generated: {date.today().isoformat()}", styles["cover_date"]))
    elems.append(Paragraph("Protein Mutation Workbench — Academic Use Only", styles["cover_date"]))
    return elems


# ── Style builder ──────────────────────────────────────────────────────────

def _build_styles() -> dict:
    base = getSampleStyleSheet()
    s = {}

    def ps(name, **kwargs):
        parent = kwargs.pop("parent", "Normal")
        p = ParagraphStyle(name, parent=base[parent], **kwargs)
        return p

    s["body"]   = ps("body",   fontSize=9,  leading=13,  textColor=colors.HexColor("#222222"))
    s["muted"]  = ps("muted",  fontSize=9,  leading=13,  textColor=_GRAY, fontName="Helvetica-Oblique")
    s["footer"] = ps("footer", fontSize=7.5, textColor=_GRAY, alignment=TA_CENTER)
    s["toc"]    = ps("toc",   fontSize=9,  leading=14, leftIndent=10)

    s["h1"] = ps("h1", parent="Heading1", fontSize=13, textColor=_HEADER_BG,
                  spaceBefore=14, spaceAfter=6, fontName="Helvetica-Bold")
    s["h2"] = ps("h2", parent="Heading2", fontSize=10, textColor=_GREEN,
                  spaceBefore=10, spaceAfter=4, fontName="Helvetica-Bold")

    s["cover_title"]   = ps("ct",  fontSize=28, alignment=TA_CENTER, leading=36,
                              fontName="Helvetica-Bold")
    s["cover_variant"] = ps("cv",  fontSize=20, alignment=TA_CENTER, leading=28,
                              fontName="Helvetica-Bold")
    s["cover_meta"]    = ps("cm",  fontSize=11, alignment=TA_CENTER, textColor=_GRAY,
                              spaceAfter=4)
    s["cover_date"]    = ps("cd",  fontSize=9,  alignment=TA_CENTER, textColor=_GRAY)

    return s


# ── Page footer (called by SimpleDocTemplate) ─────────────────────────────

def _page_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(_GRAY)
    canvas.drawString(2*cm, 1.2*cm, "Protein Mutation Workbench — Academic Use Only")
    canvas.drawRightString(A4[0] - 2*cm, 1.2*cm, f"Page {doc.page}")
    canvas.restoreState()


# ── Flowable helpers ───────────────────────────────────────────────────────

def _h1(text: str, styles) -> Paragraph:
    return Paragraph(text, styles["h1"])


def _h2(text: str, styles) -> Paragraph:
    return Paragraph(text, styles["h2"])


def _callout(text: str, styles, color=None, bold: bool = False) -> Table:
    color = color or _GREEN_LIGHT
    style = ParagraphStyle("callout", parent=getSampleStyleSheet()["Normal"],
                           fontSize=9.5, leading=14,
                           fontName="Helvetica-Bold" if bold else "Helvetica")
    tbl = Table([[Paragraph(text, style)]], colWidths=[17.5*cm])
    tbl.setStyle(TableStyle([
        ("BACKGROUND",    (0,0), (-1,-1), color),
        ("TOPPADDING",    (0,0), (-1,-1), 8),
        ("BOTTOMPADDING", (0,0), (-1,-1), 8),
        ("LEFTPADDING",   (0,0), (-1,-1), 10),
        ("RIGHTPADDING",  (0,0), (-1,-1), 10),
        ("ROUNDEDCORNERS",[6]),
    ]))
    return tbl


def _kv_table(rows: list[tuple[str, str]], styles) -> list:
    if not rows:
        return []
    data = [[Paragraph(f"<b>{k}</b>", styles["body"]),
             Paragraph(str(v) if v is not None else "—", styles["body"])]
            for k, v in rows]
    tbl = Table(data, colWidths=[6*cm, 11.5*cm])
    tbl.setStyle(TableStyle([
        ("GRID",          (0,0), (-1,-1), 0.4, colors.HexColor("#d4ddd8")),
        ("BACKGROUND",    (0,0), (0,-1),  _LIGHT_GRAY),
        ("VALIGN",        (0,0), (-1,-1), "TOP"),
        ("TOPPADDING",    (0,0), (-1,-1), 4),
        ("BOTTOMPADDING", (0,0), (-1,-1), 4),
        ("LEFTPADDING",   (0,0), (-1,-1), 6),
        ("RIGHTPADDING",  (0,0), (-1,-1), 6),
        ("ROWBACKGROUNDS",(0,0), (-1,-1), [_WHITE, _LIGHT_GRAY]),
    ]))
    return [tbl, Spacer(1, 3*mm)]


def _data_table(rows: list[list], styles, col_widths: list | None = None) -> Table:
    if not rows:
        return Paragraph("No data.", styles["muted"])
    body_style = ParagraphStyle("td", parent=getSampleStyleSheet()["Normal"], fontSize=8, leading=11)
    formatted = []
    for i, row in enumerate(rows):
        fmt_row = []
        for cell in row:
            text = str(cell) if cell is not None else "—"
            if i == 0:
                fmt_row.append(Paragraph(f"<b>{text}</b>", body_style))
            else:
                fmt_row.append(Paragraph(text, body_style))
        formatted.append(fmt_row)

    tbl = Table(formatted, colWidths=col_widths, repeatRows=1)
    tbl.setStyle(TableStyle([
        ("BACKGROUND",    (0,0), (-1,0),  _HEADER_BG),
        ("TEXTCOLOR",     (0,0), (-1,0),  _WHITE),
        ("GRID",          (0,0), (-1,-1), 0.4, colors.HexColor("#c8d4ce")),
        ("VALIGN",        (0,0), (-1,-1), "TOP"),
        ("TOPPADDING",    (0,0), (-1,-1), 4),
        ("BOTTOMPADDING", (0,0), (-1,-1), 4),
        ("LEFTPADDING",   (0,0), (-1,-1), 5),
        ("RIGHTPADDING",  (0,0), (-1,-1), 5),
        ("ROWBACKGROUNDS",(0,1), (-1,-1), [_WHITE, _LIGHT_GRAY]),
    ]))
    return tbl


# ── Chart helpers ──────────────────────────────────────────────────────────

def _metric_chart(metric_deltas) -> io.BytesIO | None:
    if not _MPL_AVAILABLE:
        return None
    items = [(m.label, m.delta) for m in metric_deltas if m.delta is not None]
    if not items:
        return None
    labels, values = zip(*items)
    colors_list = ["#e05c3a" if v > 0 else "#0b6b57" for v in values]
    fig, ax = plt.subplots(figsize=(7, max(2.5, len(labels) * 0.45)))
    bars = ax.barh(labels, values, color=colors_list, height=0.6, edgecolor="white")
    ax.axvline(0, color="#444", linewidth=0.8)
    ax.set_xlabel("Δ (mutant − wild-type)", fontsize=8)
    ax.set_title("Physicochemical metric changes", fontsize=9)
    ax.tick_params(labelsize=7.5)
    for bar, val in zip(bars, values):
        ha = "left" if val >= 0 else "right"
        ax.text(val + (0.003 if val >= 0 else -0.003),
                bar.get_y() + bar.get_height() / 2,
                f"{val:+.4f}", va="center", ha=ha, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def _ss_pie(state_cards) -> io.BytesIO | None:
    if not _MPL_AVAILABLE or not state_cards:
        return None
    cmap = {"helix": "#c65d44", "strand": "#356d9c", "coil": "#7a8a78"}
    items = [(c["label"], c["count"], cmap.get(c["class_name"], "#aaa"))
             for c in state_cards if c["count"] > 0]
    if not items:
        return None
    labels, values, colors_list = zip(*items)
    fig, ax = plt.subplots(figsize=(4, 4))
    wedges, texts, autotexts = ax.pie(
        values, labels=labels, colors=colors_list,
        autopct="%1.1f%%", startangle=90,
        textprops={"fontsize": 9}, pctdistance=0.78)
    for at in autotexts:
        at.set_fontsize(8)
    ax.set_title("Secondary structure composition", fontsize=9)
    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf
