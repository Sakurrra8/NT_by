#!/usr/bin/env python3
"""Assemble all benchmark PDFs for one case into a single indexed report."""

import argparse
import csv
import io
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


SUMMARY_FIELDS = [
    ("n_D_0", "Triangle D density"),
    ("n_D2_0", "Triangle D2 density"),
    ("T_D_0", "Triangle D temperature"),
    ("T_D2_0", "Triangle D2 temperature"),
    ("B2_n_D_0", "B2 D density"),
    ("B2_n_D2_0", "B2 D2 density"),
    ("B2_T_D_0", "B2 D temperature"),
    ("B2_T_D2_0", "B2 D2 temperature"),
    ("Tri_inner_target_r30_36_n_D2_0", "Inner-target D2 density r30-36"),
    ("SourceStratum_tri_n_D_0_total", "Triangle D source inventory"),
    ("SourceStratum_tri_n_D2_0_total", "Triangle D2 source inventory"),
    ("SourceStratum_b2_n_D_0_total", "B2 D source inventory"),
    ("SourceStratum_b2_n_D2_0_total", "B2 D2 source inventory"),
]


def load_metrics(path):
    with path.open(newline="") as stream:
        return {row["field"]: row for row in csv.DictReader(stream)}


def number(row, key):
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return float("nan")


def summary_pdf(metrics, case_label, job_id):
    buffer = io.BytesIO()
    document = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        rightMargin=16 * mm,
        leftMargin=16 * mm,
        topMargin=15 * mm,
        bottomMargin=15 * mm,
        title=f"NT and EIRENE comparison - {case_label}",
    )
    styles = getSampleStyleSheet()
    story = [
        Paragraph("NT and EIRENE Full Comparison", styles["Title"]),
        Spacer(1, 5 * mm),
        Paragraph(f"Case: {case_label}", styles["Heading2"]),
        Paragraph(f"HPC job: {job_id}", styles["Normal"]),
        Paragraph("Wall-side model: EIRENE ILSIDE=1 back-side absorption", styles["Normal"]),
        Spacer(1, 5 * mm),
        Paragraph(
            "Relative bias is (NT - EIRENE) / EIRENE. Relative L1 measures "
            "the volume-weighted spatial mismatch.",
            styles["BodyText"],
        ),
        Spacer(1, 4 * mm),
    ]
    data = [["Quantity", "NT integral", "EIRENE integral", "Bias", "Rel. L1"]]
    for field, label in SUMMARY_FIELDS:
        row = metrics.get(field)
        if not row:
            continue
        code = number(row, "code_volume_integral")
        ref = number(row, "ref_volume_integral")
        bias = number(row, "volume_weighted_rel_bias")
        l1 = number(row, "volume_weighted_rel_l1")
        data.append([
            label,
            f"{code:.4e}",
            f"{ref:.4e}",
            f"{bias:+.2%}",
            f"{l1:.2%}",
        ])
    table = Table(data, colWidths=[65 * mm, 31 * mm, 31 * mm, 22 * mm, 22 * mm], repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#27364a")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (0, -1), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#edf1f5")]),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#a9b2bd")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ALIGN", (1, 1), (-1, -1), "RIGHT"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.extend([table, Spacer(1, 5 * mm)])
    story.append(Paragraph(
        "Following pages contain every comparison figure produced by benchmark_tri.py. "
        "Each map shows NT, EIRENE, and their ratio using the same case geometry.",
        styles["BodyText"],
    ))
    document.build(story)
    buffer.seek(0)
    return buffer


def figure_group(path):
    name = path.name
    if name.startswith("InputBackground"):
        return 0, "Input background parity"
    if name.endswith("_tri_compare.pdf"):
        return 1, "Triangle-mesh fields"
    if name.startswith("B2_"):
        return 2, "B2 neutral fields"
    if "Source" in name or "sink" in name or "Ion" in name:
        return 3, "Sources and sinks"
    return 4, "Additional B2 comparisons"


def add_page_number(page, page_number):
    from reportlab.pdfgen import canvas

    width = float(page.mediabox.width)
    height = float(page.mediabox.height)
    overlay_buffer = io.BytesIO()
    overlay = canvas.Canvas(overlay_buffer, pagesize=(width, height))
    overlay.setFillColor(colors.HexColor("#4f5965"))
    overlay.setFont("Helvetica", 8)
    overlay.drawRightString(width - 10 * mm, 6 * mm, str(page_number))
    overlay.save()
    overlay_buffer.seek(0)
    page.merge_page(PdfReader(overlay_buffer).pages[0])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--case-label", required=True)
    parser.add_argument("--job-id", default="unknown")
    args = parser.parse_args()

    benchmark_dir = Path(args.benchmark_dir)
    metrics_path = benchmark_dir / "benchmark_tri_metrics.csv"
    metrics = load_metrics(metrics_path)
    figures = sorted(
        (path for path in benchmark_dir.glob("*.pdf") if path.name != Path(args.output).name),
        key=lambda path: (figure_group(path)[0], path.name),
    )
    if not figures:
        raise RuntimeError(f"No comparison PDFs found in {benchmark_dir}")

    writer = PdfWriter()
    for page in PdfReader(summary_pdf(metrics, args.case_label, args.job_id)).pages:
        writer.add_page(page)
    writer.add_outline_item("Summary", 0)

    current_group = None
    for figure in figures:
        _, group = figure_group(figure)
        start_page = len(writer.pages)
        for page in PdfReader(figure).pages:
            writer.add_page(page)
        if group != current_group:
            writer.add_outline_item(group, start_page)
            current_group = group

    for index, page in enumerate(writer.pages, start=1):
        add_page_number(page, index)
    writer.add_metadata({
        "/Title": f"NT and EIRENE comparison - {args.case_label}",
        "/Subject": "Neutral transport density, temperature, velocity, source and sink comparison",
        "/Author": "NT benchmark workflow",
    })
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as stream:
        writer.write(stream)
    print(f"wrote {output} ({len(writer.pages)} pages, {len(figures)} figures)")


if __name__ == "__main__":
    main()
