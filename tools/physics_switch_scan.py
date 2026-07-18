#!/usr/bin/env python3
"""Prepare, submit, and summarize one-factor NT/EIRENE sensitivity scans."""

import argparse
import csv
import math
import re
import subprocess
from pathlib import Path


VARIANTS = [
    ("baseline_a", {}, "independent baseline replicate"),
    ("baseline_b", {}, "independent baseline replicate"),
    ("elastic_off", {"D2ElasticModel": "0"}, "disable D2-D+ elastic collisions"),
    ("elastic_h3", {"D2ElasticModel": "1"}, "H.3 first-moment elastic model"),
    ("elastic_rejection", {"D2ElasticModel": "3"}, "H.0/H.1/H.3 with sigma*v rejection"),
    ("nncs_off", {"K_NNCs": "0"}, "disable fixed-background neutral-neutral collisions"),
    ("boundary_off", {"K_DBoundarySource": "0"}, "remove FNIY interface strata 3-5"),
    (
        "boundary_direct",
        {"DBoundaryLaunchModel": "1"},
        "launch FNIY as direct outward D instead of surface recycling",
    ),
    (
        "target_legacy",
        {"DTargetIncidentModel": "0"},
        "legacy mean target energy and incidence angle",
    ),
    (
        "target_normal",
        {"DTargetIncidentModel": "2"},
        "normal monoenergetic target incidence sensitivity",
    ),
    (
        "target_flux_reconstructed",
        {"K_DTargetSourceMode": "2"},
        "reconstruct target flux from B2 ion moments",
    ),
    (
        "wall_transparent",
        {"K_EireneWallSide": "2"},
        "pass through additional-wall back-side crossings",
    ),
    (
        "trim_legacy",
        {"K_DWTrimReflection": "0"},
        "replace differential D-on-W model with legacy reflection",
    ),
    (
        "boundary_hist_3x",
        {"numPar_flight_DBoundary": "18000"},
        "three times the FNIY source histories",
    ),
    (
        "all_hist_3x",
        {
            "numPar_flight": "150000",
            "numPar_flight_Target": "150000",
            "numPar_flight_DBoundary": "18000",
        },
        "three times all active source histories",
    ),
]


METRICS = {
    "tri_D": "n_D_0",
    "tri_D2": "n_D2_0",
    "tri_TD": "T_D_0",
    "tri_TD2": "T_D2_0",
    "tri_inner_D2": "Tri_inner_target_r30_36_n_D2_0",
    "b2_D": "B2_n_D_0",
    "b2_D2": "B2_n_D2_0",
    "b2_TD": "B2_T_D_0",
    "b2_TD2": "B2_T_D2_0",
    "tri_mean_TD": "Tri_density_weighted_mean_T_D_0",
    "tri_mean_TD2": "Tri_density_weighted_mean_T_D2_0",
    "b2_mean_TD": "B2_density_weighted_mean_T_D_0",
    "b2_mean_TD2": "B2_density_weighted_mean_T_D2_0",
    "source_tri_D": "SourceStratum_tri_n_D_0_total",
    "source_tri_D2": "SourceStratum_tri_n_D2_0_total",
    "source_b2_D": "SourceStratum_b2_n_D_0_total",
    "source_b2_D2": "SourceStratum_b2_n_D2_0_total",
    "target_D2_shape": "Target_D2_source_shape_inner",
}


def replace_settings(text, changes):
    for key, value in changes.items():
        pattern = re.compile(rf"^(\s*{re.escape(key)}\s+)\S+", re.MULTILINE)
        text, count = pattern.subn(rf"\g<1>{value}", text, count=1)
        if count != 1:
            raise ValueError(f"setting {key!r} was found {count} times")
    return text


def write_manifest(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def prepare(args):
    root = Path(args.root).resolve()
    base = (root / args.base).resolve()
    outdir = (root / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    base_text = base.read_text()
    rows = []
    for label, changes, note in VARIANTS:
        setting = outdir / f"setting_{label}_tmp.log"
        setting.write_text(replace_settings(base_text, changes))
        rows.append(
            {
                "label": label,
                "setting_file": setting.relative_to(root).as_posix(),
                "changes": ";".join(f"{key}={value}" for key, value in changes.items()) or "baseline",
                "note": note,
                "simulation_job": "",
                "benchmark_job": "",
            }
        )
    manifest = outdir / "manifest.csv"
    write_manifest(manifest, rows)
    print(f"wrote {manifest}")
    return root, manifest, rows


def sbatch(root, arguments, dry_run):
    command = ["sbatch", "--parsable", "--chdir", str(root), *arguments]
    if dry_run:
        print(" ".join(command))
        return "dry-run"
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    return result.stdout.strip().split(";", 1)[0]


def submit(args):
    root, manifest, rows = prepare(args)
    for row in rows:
        simulation = sbatch(
            root,
            ["crslum.slurm", row["setting_file"]],
            args.dry_run,
        )
        benchmark = sbatch(
            root,
            [
                f"--dependency=afterok:{simulation}",
                "benchmark.slurm",
                args.case,
                simulation,
                f"5e19 switch scan: {row['label']}",
                row["setting_file"],
            ],
            args.dry_run,
        )
        row["simulation_job"] = simulation
        row["benchmark_job"] = benchmark
        write_manifest(manifest, rows)
        print(f"{row['label']}: simulation={simulation}, benchmark={benchmark}")


def read_csv(path, key):
    with path.open(newline="") as stream:
        return {row[key]: row for row in csv.DictReader(stream)}


def number(row, *keys):
    for key in keys:
        try:
            value = float(row.get(key, ""))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            return value
    return math.nan


def mean_finite(values):
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else math.nan


def wall_back_flux(path):
    if not path.exists():
        return math.nan
    total = 0.0
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row.get("side") == "back":
                total += float(row["incident_flux_s-1"])
    return total


def loss_weight(path, loss_type):
    if not path.exists():
        return math.nan
    rows = read_csv(path, "loss_type")
    return number(rows.get(loss_type, {}), "weight_s-1")


def summarize(args):
    root = Path(args.root).resolve()
    manifest = (root / args.manifest).resolve()
    with manifest.open(newline="") as stream:
        runs = list(csv.DictReader(stream))
    summaries = []
    for run in runs:
        job = run["simulation_job"]
        metrics_path = root / "Outputfile" / job / "fig/tri_benchmark/benchmark_tri_metrics.csv"
        if not job or not metrics_path.exists():
            print(f"skip incomplete {run['label']}: {metrics_path}")
            continue
        metrics = read_csv(metrics_path, "field")
        row = dict(run)
        for alias, field in METRICS.items():
            metric = metrics.get(field, {})
            row[f"{alias}_bias_pct"] = 100.0 * number(
                metric, "volume_weighted_rel_bias", "relative_bias"
            )
            row[f"{alias}_l1_pct"] = 100.0 * number(
                metric, "volume_weighted_rel_l1", "relative_l1"
            )
        data = root / "Outputfile" / job / "data"
        row["D_core_loss_s-1"] = loss_weight(data / "D_transport_loss_summary.csv", "core_loss")
        row["D_additional_no_hit_s-1"] = loss_weight(
            data / "D_transport_loss_summary.csv", "additional_cell_no_hit"
        )
        row["D2_additional_no_hit_s-1"] = loss_weight(
            data / "D2_transport_loss_summary.csv", "additional_cell_no_hit"
        )
        row["D_back_flux_s-1"] = wall_back_flux(data / "wall_side_impact_D.csv")
        row["D2_back_flux_s-1"] = wall_back_flux(data / "wall_side_impact_D2.csv")
        row["source_score_pct"] = mean_finite(
            abs(row[f"source_{mesh}_{species}_bias_pct"])
            for mesh in ("tri", "b2")
            for species in ("D", "D2")
        )
        row["density_spatial_score_pct"] = mean_finite(
            row[f"{mesh}_{species}_l1_pct"]
            for mesh in ("tri", "b2")
            for species in ("D", "D2")
        )
        row["temperature_bias_score_pct"] = mean_finite(
            abs(row[f"{mesh}_mean_{field}_bias_pct"])
            for mesh in ("tri", "b2")
            for field in ("TD", "TD2")
        )
        row["overall_diagnostic_score_pct"] = mean_finite(
            [
                row["source_score_pct"],
                row["density_spatial_score_pct"],
                abs(row["tri_inner_D2_bias_pct"]),
            ]
        )
        summaries.append(row)

    if not summaries:
        raise SystemExit("no completed benchmark results found")
    summaries.sort(key=lambda row: row["overall_diagnostic_score_pct"])
    output = (root / args.output).resolve()
    write_manifest(output, summaries)
    print(f"wrote {output}")
    print("rank label                 score source spatial innerD2  triD2   b2D2")
    for rank, row in enumerate(summaries, start=1):
        print(
            f"{rank:>4} {row['label']:<21} "
            f"{row['overall_diagnostic_score_pct']:>5.1f} "
            f"{row['source_score_pct']:>6.1f} "
            f"{row['density_spatial_score_pct']:>7.1f} "
            f"{row['tri_inner_D2_bias_pct']:>+7.1f} "
            f"{row['tri_D2_bias_pct']:>+7.1f} "
            f"{row['b2_D2_bias_pct']:>+7.1f}"
        )


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)
    for command in ("prepare", "submit"):
        sub = subparsers.add_parser(command)
        sub.add_argument("--root", default=".")
        sub.add_argument("--base", default="Inputfile/settingfile/setting_Trimesh_D_5.log")
        sub.add_argument("--outdir", default="tmp/physics_switch_scan")
        if command == "submit":
            sub.add_argument("--case", default="case_input/2MW-5e19")
            sub.add_argument("--dry-run", action="store_true")
            sub.set_defaults(func=submit)
        else:
            sub.set_defaults(func=prepare)
    sub = subparsers.add_parser("summarize")
    sub.add_argument("--root", default=".")
    sub.add_argument("--manifest", default="tmp/physics_switch_scan/manifest.csv")
    sub.add_argument("--output", default="tmp/physics_switch_scan/summary.csv")
    sub.set_defaults(func=summarize)
    return result


if __name__ == "__main__":
    arguments = parser().parse_args()
    arguments.func(arguments)
