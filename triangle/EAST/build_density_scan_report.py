#!/usr/bin/env python3
"""Build density-trend summary and optionally prepend it to all case reports."""

import argparse
import csv
import io
import math
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import Flowable, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


FIELDS = {
    "tri_D": "n_D_0",
    "tri_D2": "n_D2_0",
    "b2_D": "B2_n_D_0",
    "b2_D2": "B2_n_D2_0",
    "inner_D2": "Tri_inner_target_r30_36_n_D2_0",
    "tri_T_D": "T_D_0",
    "tri_T_D2": "T_D2_0",
    "tri_source_D": "SourceStratum_tri_n_D_0_total",
    "tri_source_D2": "SourceStratum_tri_n_D2_0_total",
}


def load_case(specification):
    label, raw_path = specification.split("=", 1)
    metrics_path = Path(raw_path) / "benchmark_tri_metrics.csv"
    with metrics_path.open(newline="") as stream:
        rows = {row["field"]: row for row in csv.DictReader(stream)}
    density = float(label.replace("e19", ""))
    return density, label, Path(raw_path), rows


def value(rows, field, key="volume_weighted_rel_bias"):
    row = rows.get(field)
    if not row or not row.get(key):
        return float("nan")
    return float(row[key])


def percent(entry):
    return f"{entry:+.1%}" if math.isfinite(entry) else "n/a"


class DensityTrendChart(Flowable):
    def __init__(self, cases, width=245 * mm, height=145 * mm):
        super().__init__()
        self.cases = cases
        self.width = width
        self.height = height

    def draw_panel(self, canvas, x, y, width, height, series, title):
        plot_x = x + 14 * mm
        plot_y = y + 12 * mm
        plot_width = width - 18 * mm
        plot_height = height - 23 * mm
        densities = [case[0] for case in self.cases]
        all_values = [
            100.0 * value(case[3], FIELDS[key])
            for _, key in series
            for case in self.cases
        ]
        finite = [entry for entry in all_values if entry == entry]
        lower = min(finite + [0.0])
        upper = max(finite + [0.0])
        margin = max(5.0, 0.12 * max(upper - lower, 1.0))
        lower -= margin
        upper += margin
        density_min = min(densities)
        density_max = max(densities)
        density_span = max(density_max - density_min, 1.e-12)
        value_span = max(upper - lower, 1.e-12)

        def map_x(entry):
            return plot_x + plot_width * (entry - density_min) / density_span

        def map_y(entry):
            return plot_y + plot_height * (entry - lower) / value_span

        canvas.setFillColor(colors.HexColor("#27364a"))
        canvas.setFont("Helvetica-Bold", 9)
        canvas.drawString(x, y + height - 5 * mm, title)
        canvas.setStrokeColor(colors.HexColor("#d1d7de"))
        canvas.setLineWidth(0.4)
        for tick in range(5):
            tick_value = lower + tick * value_span / 4.0
            tick_y = map_y(tick_value)
            canvas.line(plot_x, tick_y, plot_x + plot_width, tick_y)
            canvas.setFillColor(colors.HexColor("#4f5965"))
            canvas.setFont("Helvetica", 6.5)
            canvas.drawRightString(plot_x - 1.5 * mm, tick_y - 1.8, f"{tick_value:.0f}")
        canvas.setStrokeColor(colors.HexColor("#4f5965"))
        canvas.line(plot_x, plot_y, plot_x, plot_y + plot_height)
        canvas.line(plot_x, plot_y, plot_x + plot_width, plot_y)
        if lower <= 0.0 <= upper:
            canvas.setStrokeColor(colors.HexColor("#7b8490"))
            canvas.line(plot_x, map_y(0.0), plot_x + plot_width, map_y(0.0))
        for density in densities:
            tick_x = map_x(density)
            canvas.setFillColor(colors.HexColor("#4f5965"))
            canvas.setFont("Helvetica", 6.5)
            canvas.drawCentredString(tick_x, plot_y - 3 * mm, f"{density:g}")

        palette = [colors.HexColor("#087f8c"), colors.HexColor("#d1495b")]
        for index, (label, key) in enumerate(series):
            points = []
            for case in self.cases:
                entry = 100.0 * value(case[3], FIELDS[key])
                if math.isfinite(entry):
                    points.append((map_x(case[0]), map_y(entry)))
            canvas.setStrokeColor(palette[index])
            canvas.setFillColor(palette[index])
            canvas.setLineWidth(1.4)
            for first, second in zip(points, points[1:]):
                canvas.line(first[0], first[1], second[0], second[1])
            for point_x, point_y in points:
                canvas.circle(point_x, point_y, 1.5, fill=1, stroke=0)
            legend_x = plot_x + index * 34 * mm
            legend_y = y + 2.5 * mm
            canvas.line(legend_x, legend_y, legend_x + 5 * mm, legend_y)
            canvas.setFillColor(colors.HexColor("#27364a"))
            canvas.setFont("Helvetica", 6.5)
            canvas.drawString(legend_x + 6 * mm, legend_y - 2, label)

    def draw(self):
        gap = 7 * mm
        panel_width = (self.width - gap) / 2.0
        panel_height = (self.height - gap) / 2.0
        panels = [
            (0, panel_height + gap, [("Triangle D", "tri_D"), ("Triangle D2", "tri_D2")], "Triangle density bias"),
            (panel_width + gap, panel_height + gap, [("B2 D", "b2_D"), ("B2 D2", "b2_D2")], "B2 density bias"),
            (0, 0, [("Inner target D2", "inner_D2")], "Inner-target D2 bias (r30-36)"),
            (panel_width + gap, 0, [("D", "tri_source_D"), ("D2", "tri_source_D2")], "Triangle source-inventory bias"),
        ]
        for x, y, series, title in panels:
            self.draw_panel(self.canv, x, y, panel_width, panel_height, series, title)


def summary_pdf(cases):
    buffer = io.BytesIO()
    document = SimpleDocTemplate(
        buffer,
        pagesize=landscape(A4),
        rightMargin=14 * mm,
        leftMargin=14 * mm,
        topMargin=12 * mm,
        bottomMargin=12 * mm,
        title="NT and EIRENE density scan",
    )
    styles = getSampleStyleSheet()
    story = [
        Paragraph("NT and EIRENE Density Scan", styles["Title"]),
        Paragraph(
            "All cases use EIRENE ILSIDE=1 back-side absorption. Values below are "
            "volume-weighted relative biases: (NT - EIRENE) / EIRENE.",
            styles["BodyText"],
        ),
        Spacer(1, 3 * mm),
        DensityTrendChart(cases),
    ]
    table_data = [["Density", "Tri D", "Tri D2", "B2 D", "B2 D2", "Inner D2", "Tri T(D)", "Tri T(D2)"]]
    for _, label, _, rows in cases:
        keys = ["tri_D", "tri_D2", "b2_D", "b2_D2", "inner_D2", "tri_T_D", "tri_T_D2"]
        table_data.append([label] + [percent(value(rows, FIELDS[key])) for key in keys])
    table = Table(table_data, colWidths=[28 * mm] + [30 * mm] * 7, repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#27364a")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#edf1f5")]),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#a9b2bd")),
        ("ALIGN", (1, 1), (-1, -1), "RIGHT"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.extend([PageBreak(), Paragraph("Per-case summary", styles["Heading2"]), table])
    document.build(story)
    buffer.seek(0)
    return buffer


def write_pdf(reader, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    writer = PdfWriter()
    for page in reader.pages:
        writer.add_page(page)
    with output.open("wb") as stream:
        writer.write(stream)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", required=True, help="DENSITY=benchmark_directory")
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--full-output")
    args = parser.parse_args()

    cases = sorted((load_case(specification) for specification in args.case), key=lambda case: case[0])
    summary_reader = PdfReader(summary_pdf(cases))
    write_pdf(summary_reader, Path(args.summary_output))
    print(f"wrote {args.summary_output} ({len(summary_reader.pages)} pages)")

    if args.full_output:
        writer = PdfWriter()
        for page in summary_reader.pages:
            writer.add_page(page)
        writer.add_outline_item("Density scan summary", 0)
        for _, label, benchmark_dir, _ in cases:
            report = benchmark_dir / "eirene_comparison_full.pdf"
            start_page = len(writer.pages)
            for page in PdfReader(report).pages:
                writer.add_page(page)
            writer.add_outline_item(label, start_page)
        output = Path(args.full_output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("wb") as stream:
            writer.write(stream)
        print(f"wrote {output} ({len(writer.pages)} pages)")


if __name__ == "__main__":
    main()
