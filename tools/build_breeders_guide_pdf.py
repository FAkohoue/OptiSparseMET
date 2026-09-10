"""Build the September 2026 production edition of the Breeder Guide PDF.

The editable source for the original guide is not retained in the working tree.
This builder preserves the original substantive pages, replaces the former
blank page 18, appends two production-engine reference pages, and refreshes the
continued contents page. It is idempotent for either the original 20-page input
or the generated 22-page edition.
"""

from __future__ import annotations

import argparse
import io
from datetime import datetime
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from reportlab.lib.colors import Color, HexColor, white
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph


PAGE_WIDTH = 612
PAGE_HEIGHT = 792
NAVY = HexColor("#0C4C80")
BLUE = HexColor("#2878C7")
TEXT = HexColor("#17324D")
MUTED = HexColor("#5F7082")


def _paragraph_style(name: str, **overrides) -> ParagraphStyle:
    values = {
        "fontName": "Helvetica",
        "fontSize": 9.1,
        "leading": 11.8,
        "textColor": TEXT,
        "alignment": TA_LEFT,
        "spaceAfter": 0,
        "allowWidows": 0,
        "allowOrphans": 0,
    }
    values.update(overrides)
    return ParagraphStyle(name, **values)


BODY = _paragraph_style("Body")
INTRO = _paragraph_style("Intro", fontSize=9.4, leading=12.2)
HEADING = _paragraph_style(
    "Heading", fontName="Helvetica-Bold", fontSize=10.9, leading=13.2,
    textColor=NAVY,
)
BULLET = _paragraph_style(
    "Bullet", leftIndent=12, firstLineIndent=-8, fontSize=8.85, leading=11.35,
)
SMALL = _paragraph_style("Small", fontSize=8.45, leading=10.7, textColor=MUTED)


def _draw_page_frame(c: canvas.Canvas, page_number: int) -> None:
    c.setFillColor(MUTED)
    c.setFont("Helvetica", 8.1)
    c.drawString(72, 756, "OptiSparseMET  |  Breeder's Guide")
    c.drawRightString(
        540, 35, f"OptiSparseMET 0.2.0   |   {page_number}"
    )


def _draw_paragraph(c: canvas.Canvas, text: str, style: ParagraphStyle,
                    x: float, y: float, width: float) -> float:
    paragraph = Paragraph(text, style)
    _, height = paragraph.wrap(width, PAGE_HEIGHT)
    paragraph.drawOn(c, x, y - height)
    return y - height


def _cover_overlay() -> bytes:
    stream = io.BytesIO()
    c = canvas.Canvas(stream, pagesize=(PAGE_WIDTH, PAGE_HEIGHT))

    # Replace only the date line while preserving the original cover design.
    c.setFillColor(white)
    c.rect(244, 350, 126, 23, fill=1, stroke=0)
    c.setFillColor(MUTED)
    c.setFont("Helvetica", 8.1)
    c.drawCentredString(PAGE_WIDTH / 2, 360.5, "Version 0.2.0  |  September 2026")
    c.save()
    return stream.getvalue()


def _environment_update_page() -> bytes:
    stream = io.BytesIO()
    c = canvas.Canvas(stream, pagesize=(PAGE_WIDTH, PAGE_HEIGHT))
    _draw_page_frame(c, 18)
    x = 72
    width = 468
    y = 711

    c.setFillColor(NAVY)
    c.setFont("Helvetica-Bold", 16)
    c.drawString(x, y, "12.5 Environmental characterization additions")
    y -= 26

    y = _draw_paragraph(
        c,
        "The environmental workflow now separates descriptive structure from "
        "validated hard grouping and exposes the evidence behind that decision. "
        "The additions below preserve the package's conservative interpretation "
        "while making historical weather and SoilGrids analyses easier to audit.",
        INTRO, x, y, width,
    ) - 13

    y = _draw_paragraph(c, "Candidate discovery and validation", HEADING, x, y, width) - 4
    bullets = [
        "- <b>infer_mega_environments()</b> keeps strict defaults for hard groups. "
        "If validation fails, <b>membership</b> remains one unpartitioned network, "
        "while the best estimable structure is retained in the <b>candidate_*</b> fields.",
        "- <b>infer_environmental_strata()</b> returns the descriptive candidate and "
        "permits singleton strata by default. These strata are not genetic mega-"
        "environments unless historical GxE responses support that interpretation.",
        "- Block-specific agreement and adjusted-Rand diagnostics show whether "
        "instability comes from weather years, soil, another modality, or conflict "
        "between evidence blocks.",
    ]
    for bullet in bullets:
        y = _draw_paragraph(c, bullet, BULLET, x, y, width) - 3
    y -= 4

    y = _draw_paragraph(c, "Request-safe weather and soil metadata", HEADING, x, y, width) - 4
    bullets = [
        "- <b>available_weather_parameters()</b> returns supported NASA POWER API "
        "codes. The catalogue now separates <b>api_code</b>, <b>output_name</b>, "
        "and advisory <b>default_aggregation</b>, so retrieval and interpretation "
        "use the correct names.",
        "- SoilGrids requests are validated as property-by-depth combinations. Most "
        "properties use the six standard depth intervals; organic carbon stock "
        "(<b>ocs</b>) is available at 0-30 cm. Unsupported combinations stop with an "
        "informative error before network retrieval.",
    ]
    for bullet in bullets:
        y = _draw_paragraph(c, bullet, BULLET, x, y, width) - 3
    y -= 4

    y = _draw_paragraph(c, "Optional within-block redundancy control", HEADING, x, y, width) - 4
    y = _draw_paragraph(
        c,
        "<b>build_environment_kernels()</b> retains <b>redundancy = \"none\"</b> "
        "as the default and adds correlation filtering, PCA, and whitening as "
        "explicit alternatives. The returned audit preserves original covariates "
        "and records retained variables or components, correlation groups, "
        "loadings, variance explained, and effective rank before and after control.",
        BODY, x, y, width,
    ) - 13

    y = _draw_paragraph(c, "One audited historical weather and soil workflow", HEADING, x, y, width) - 4
    y = _draw_paragraph(
        c,
        "<b>historical_environment_characterization()</b> connects crop-window "
        "weather histories, optional SoilGrids profiles, temporal stability, "
        "separate modality kernels, weather-soil agreement, descriptive strata, "
        "strict validation, and provenance/QC in one result. Without a historical "
        "genetic-response target, it uses a weight-free descriptive consensus; "
        "with a defensible target, the existing covariance calibration can estimate "
        "supported modality contributions.",
        BODY, x, y, width,
    ) - 13

    y = _draw_paragraph(c, "How to act on the result", HEADING, x, y, width) - 4
    actions = [
        "- <b>candidate_k &gt; 1 and hard_groups = FALSE:</b> use the partition for "
        "description and sensitivity analysis, not as a compulsory allocation rule.",
        "- <b>Stable weather but low cross-modality agreement:</b> the weather "
        "pattern is repeatable, but soil describes a different geometry; inspect "
        "the block diagnostics rather than lowering validation thresholds.",
        "- <b>hard_groups = TRUE:</b> environmental validation passed. Continue to "
        "require historical GxE evidence before calling the groups genetic mega-"
        "environments.",
    ]
    for action in actions:
        y = _draw_paragraph(c, action, BULLET, x, y, width) - 3

    # A compact provenance note; the original page header and footer remain visible.
    y -= 4
    c.setStrokeColor(Color(0.82, 0.86, 0.90))
    c.setLineWidth(0.5)
    c.line(x, y, x + width, y)
    y -= 10
    _draw_paragraph(
        c,
        "Implementation update: September 2026. See the package function reference "
        "for arguments, return fields, and reproducible examples.",
        SMALL, x, y, width,
    )

    c.save()
    return stream.getvalue()


def _production_engine_page() -> bytes:
    stream = io.BytesIO()
    c = canvas.Canvas(stream, pagesize=(PAGE_WIDTH, PAGE_HEIGHT))
    _draw_page_frame(c, 21)
    x = 72
    width = 468
    y = 711

    c.setFillColor(NAVY)
    c.setFont("Helvetica-Bold", 16)
    c.drawString(x, y, "Appendix A Production engine additions")
    y -= 26
    y = _draw_paragraph(
        c,
        "The final design calculation can now use the breeding programme's "
        "historical response covariance, target-environment priorities, local "
        "field geometry, residual precision, and physical costs in one objective.",
        INTRO, x, y, width,
    ) - 13

    y = _draw_paragraph(c, "Historical genetic covariance", HEADING, x, y, width) - 4
    y = _draw_paragraph(
        c,
        "<b>fit_historical_met()</b> fits diagonal, factor-analytic, or "
        "unstructured environment covariance by REML. Missing genotype-by-"
        "environment cells remain in the likelihood. Supply residual variance "
        "when cells are unreplicated; otherwise replicated cells estimate it. "
        "Compare convergence, AIC, BIC, blocked prediction, and biological "
        "plausibility before selecting the central covariance.",
        BODY, x, y, width,
    ) - 13

    y = _draw_paragraph(c, "Programme target and site precision", HEADING, x, y, width) - 4
    bullets = [
        "- <b>tpe_weights</b> records the relative importance of environments "
        "in the target population and is normalised to sum to one.",
        "- <b>sigma_e2</b> may differ by environment. Do not use a common value "
        "when historical trial analyses support material precision differences.",
        "- TPE weights are programme priorities. They are not environmental-"
        "kernel weights and are not inferred from the diagonal of the covariance.",
    ]
    for bullet in bullets:
        y = _draw_paragraph(c, bullet, BULLET, x, y, width) - 3
    y -= 4

    y = _draw_paragraph(c, "Joint allocation and field layout", HEADING, x, y, width) - 4
    y = _draw_paragraph(
        c,
        "<b>fieldbook_design_evaluator()</b> converts every candidate allocation "
        "into actual local fieldbooks. It returns integer replication, repeated-"
        "check overhead, site-specific cost, and full treatment-information "
        "matrices from <b>local_treatment_information()</b>. "
        "<b>optimize_design()</b> therefore accepts or rejects an allocation "
        "using the plantable layouts it produces, not a fixed efficiency guess.",
        BODY, x, y, width,
    ) - 13

    y = _draw_paragraph(c, "Large-network calculation", HEADING, x, y, width) - 4
    y = _draw_paragraph(
        c,
        "<b>met_information(solver = \"auto\")</b> uses the exact dense solve "
        "for moderate problems and matrix-free preconditioned conjugate gradients "
        "(PCG) above the declared dimension threshold. PCG avoids the full "
        "JE-by-JE matrix and estimates target PEV diagonals with deterministic "
        "probes. Record the probe count, tolerance, iterations, residual norms, "
        "and convergence fraction, and confirm design rankings are stable as the "
        "probe count increases.",
        BODY, x, y, width,
    )

    c.save()
    return stream.getvalue()


def _release_record_page() -> bytes:
    stream = io.BytesIO()
    c = canvas.Canvas(stream, pagesize=(PAGE_WIDTH, PAGE_HEIGHT))
    _draw_page_frame(c, 22)
    x = 72
    width = 468
    y = 711

    c.setFillColor(NAVY)
    c.setFont("Helvetica-Bold", 16)
    c.drawString(x, y, "Appendix B Release checks and design record")
    y -= 26
    y = _draw_paragraph(
        c,
        "The released fieldbook must reproduce the design that was scored. The "
        "package now treats allocation, realised plot counts, covariance "
        "assumptions, diagnostics, and provenance as one versioned record.",
        INTRO, x, y, width,
    ) - 13

    y = _draw_paragraph(c, "Integer replication and physical cost", HEADING, x, y, width) - 4
    y = _draw_paragraph(
        c,
        "<b>recommend_replication()</b> converts fractional p-rep targets to "
        "balanced integer plots. A target of 1.5 gives one plot to half of the "
        "occupied cells and two to the remainder. <b>design_objective()</b> "
        "accepts environment-specific cost per plot and fixed overhead such as "
        "repeated checks. Both candidate and overhead plots count against the "
        "budget.",
        BODY, x, y, width,
    ) - 13

    y = _draw_paragraph(c, "Optimisation evidence", HEADING, x, y, width) - 4
    bullets = [
        "- Review proposed, accepted, and improving moves and the acceptance rate.",
        "- Investigate infeasible proposals and fieldbook-evaluator failures.",
        "- Review unique candidate evaluations, cache hits, and every restart "
        "trajectory; retain the random seed and objective definition.",
    ]
    for bullet in bullets:
        y = _draw_paragraph(c, bullet, BULLET, x, y, width) - 3
    y -= 4

    y = _draw_paragraph(c, "Validated release object", HEADING, x, y, width) - 4
    y = _draw_paragraph(
        c,
        "<b>sparse_met_design()</b> stores the treatment-by-environment "
        "allocation, integer replication, fieldbooks, G, Sigma_E, TPE weights, "
        "variance assumptions, diagnostics, and provenance. "
        "<b>validate_sparse_met_design()</b> checks identifier alignment, "
        "covariance symmetry and positive semidefiniteness, and exact agreement "
        "between each fieldbook and its replication matrix.",
        BODY, x, y, width,
    ) - 13

    y = _draw_paragraph(c, "Release checklist", HEADING, x, y, width) - 4
    actions = [
        "- Every allocated treatment appears in the correct environment and at "
        "the exact recorded integer replication.",
        "- Check plots and other fixed overhead are included in capacity and cost.",
        "- The network-wide seed ledger remains non-negative after the reserve.",
        "- Dense or PCG information diagnostics meet the declared numerical rule.",
        "- Robustness and benchmarking cover covariance, variance, missing data, "
        "site loss, and explicit reference designs.",
        "- The validated design object, fieldbooks, seed, software version, and "
        "analysis assumptions are archived together.",
    ]
    for action in actions:
        y = _draw_paragraph(c, action, BULLET, x, y, width) - 3

    c.save()
    return stream.getvalue()


def _contents_overlay() -> bytes:
    stream = io.BytesIO()
    c = canvas.Canvas(stream, pagesize=(PAGE_WIDTH, PAGE_HEIGHT))
    c.setFillColor(white)
    c.rect(48, 54, 516, 682, fill=1, stroke=0)
    x = 72
    y = 711
    c.setFillColor(NAVY)
    c.setFont("Helvetica-Bold", 16)
    c.drawString(x, y, "Contents continued")
    y -= 25
    entries = [
        ("8.2 Equireplicate allocation", 12, False),
        ("8.3 Colmant contribution optimisation", 13, False),
        ("8.4 Genetic and environmental structure in allocation", 13, False),
        ("9. Network-wide seed accounting", 13, True),
        ("9.1 Reserve policy", 14, False),
        ("9.2 Allocation and replication sequence", 14, False),
        ("9.3 Audit fields", 14, False),
        ("10. Within-environment field design", 14, True),
        ("10.1 Block-based repeated-check designs", 14, False),
        ("10.2 Alpha row-column stream designs", 14, False),
        ("10.3 Local design checklist", 15, False),
        ("11. Statistical evaluation and optimisation", 15, True),
        ("11.1 Information-matrix criteria", 15, False),
        ("11.2 Outcome simulation", 15, False),
        ("11.3 Optimisation methods", 15, False),
        ("11.4 Robust aggregation", 16, False),
        ("12. Benchmarking, capacity, Pareto frontiers, and release", 16, True),
        ("12.1 Site-capacity selection", 16, False),
        ("12.2 Benchmarking and validation", 16, False),
        ("12.3 Pareto decisions", 17, False),
        ("12.4 Pre-release checklist", 17, False),
        ("12.5 Environmental characterization additions", 18, False),
        ("13. Decision guide", 19, True),
        ("14. Glossary", 19, True),
        ("15. References", 20, True),
        ("16. Final perspective", 20, True),
        ("Appendix A Production engine additions", 21, True),
        ("Appendix B Release checks and design record", 22, True),
    ]
    for label, page, major in entries:
        c.setFillColor(TEXT)
        c.setFont("Helvetica-Bold" if major else "Helvetica", 8.55)
        c.drawString(x, y, label)
        c.drawRightString(540, y, str(page))
        y -= 16.2
    y -= 2
    _draw_paragraph(
        c,
        "How to use this guide. Read Sections 1-3 before choosing inputs. Use "
        "Sections 4-9 while building the design. Use Sections 10-12 to evaluate "
        "and release the fieldbooks. Section 12.5 and Appendices A-B document "
        "the September 2026 production additions.",
        SMALL, x, y, 468,
    )
    c.save()
    return stream.getvalue()


def build(source: Path, output: Path) -> None:
    reader = PdfReader(source)
    if len(reader.pages) not in (20, 22):
        raise ValueError(
            f"Expected the original 20-page or generated 22-page guide, "
            f"found {len(reader.pages)} pages"
        )

    cover = PdfReader(io.BytesIO(_cover_overlay())).pages[0]
    contents = PdfReader(io.BytesIO(_contents_overlay())).pages[0]
    environment_update = PdfReader(
        io.BytesIO(_environment_update_page())
    ).pages[0]
    production_update = PdfReader(
        io.BytesIO(_production_engine_page())
    ).pages[0]
    release_update = PdfReader(io.BytesIO(_release_record_page())).pages[0]

    writer = PdfWriter()
    for index, page in enumerate(reader.pages[:17]):
        if index == 0:
            page.merge_page(cover, over=True)
        if index == 2:
            page.merge_page(contents, over=True)
        writer.add_page(page)
    writer.add_page(environment_update)
    for page in reader.pages[18:20]:
        writer.add_page(page)
    writer.add_page(production_update)
    writer.add_page(release_update)

    metadata = dict(reader.metadata or {})
    metadata.update({
        "/Title": "The OptiSparseMET Breeder's Guide",
        "/Author": "Félicien Akohoue",
        "/Subject": "Designing, evaluating, and releasing connected and seed-feasible sparse multi-environment trials",
        "/Keywords": "OptiSparseMET, sparse MET, REML, factor analytic, PCG, field design, breeding trials",
        "/Creator": "OptiSparseMET Breeder Guide builder",
        "/Producer": "pypdf and ReportLab",
        "/ModDate": datetime.now().astimezone().strftime("D:%Y%m%d%H%M%S%z"),
    })
    writer.add_metadata(metadata)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as handle:
        writer.write(handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    build(args.source, args.output)


if __name__ == "__main__":
    main()
