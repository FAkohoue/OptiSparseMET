"""Add the environmental-improvements page to the Breeder Guide PDF.

The editable Word source for the current guide is not retained in the working
tree.  This builder therefore preserves the original PDF pages and places the
September 2026 environmental update on the intentionally blank page 18.
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


def _environment_update_overlay() -> bytes:
    stream = io.BytesIO()
    c = canvas.Canvas(stream, pagesize=(PAGE_WIDTH, PAGE_HEIGHT))
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


def build(source: Path, output: Path) -> None:
    reader = PdfReader(source)
    if len(reader.pages) != 20:
        raise ValueError(f"Expected the 20-page guide, found {len(reader.pages)} pages")

    cover = PdfReader(io.BytesIO(_cover_overlay())).pages[0]
    update = PdfReader(io.BytesIO(_environment_update_overlay())).pages[0]

    writer = PdfWriter()
    for index, page in enumerate(reader.pages):
        if index == 0:
            page.merge_page(cover, over=True)
        if index == 17:
            page.merge_page(update, over=True)
        writer.add_page(page)

    metadata = dict(reader.metadata or {})
    metadata.update({
        "/Title": "The OptiSparseMET Breeder's Guide",
        "/Author": "Félicien Akohoue",
        "/Subject": "Designing connected, seed-feasible, and statistically robust sparse multi-environment trials",
        "/Keywords": "OptiSparseMET, sparse MET, environmental characterization, breeding trials",
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
