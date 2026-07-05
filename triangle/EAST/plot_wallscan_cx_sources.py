#!/usr/bin/env python3
"""Plot density-scan CX source terms on the EAST triangle geometry."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

MPL_CACHE = Path(tempfile.gettempdir()) / "nt-matplotlib"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


CASE_DENSITIES = {
    1: 3.2e19,
    2: 3.6e19,
    3: 4.0e19,
    4: 4.6e19,
    5: 5.0e19,
    6: 5.4e19,
    7: 5.5e19,
    8: 5.6e19,
    9: 5.8e19,
}

CASE_PATHS = {
    1: "2MW-3.2e19",
    2: "2MW-3.6e19",
    3: "2MW-4e19",
    4: "2MW-4.6e19",
    5: "2MW-5e19",
    6: "2MW-5.4e19",
    7: "2MW-5.5e19",
    8: "2MW-5.6e19",
    9: "2MW-5.8e19",
}

DEFAULT_RUN_PREFIX = "wallscan_recyc500k_D"


def read_fort33(path: Path) -> tuple[np.ndarray, np.ndarray]:
    values = path.read_text().split()
    count = int(values[0])
    r = np.asarray(values[1 : 1 + count], dtype=float) / 100.0
    z = np.asarray(values[1 + count : 1 + 2 * count], dtype=float) / 100.0
    return r, z


def read_fort34(path: Path) -> np.ndarray:
    return np.loadtxt(path, dtype=int, skiprows=1)[:, 1:4] - 1


def triangle_volume(r: np.ndarray, z: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    rr = r[triangles]
    zz = z[triangles]
    area = 0.5 * np.abs(
        (rr[:, 1] - rr[:, 0]) * (zz[:, 2] - zz[:, 0])
        - (rr[:, 2] - rr[:, 0]) * (zz[:, 1] - zz[:, 0])
    )
    return 2.0 * np.pi * rr.mean(axis=1) * area


def load_geometry(root: Path, case: int) -> tuple[mtri.Triangulation, np.ndarray | None, np.ndarray]:
    case_dir = root / "case_input" / CASE_PATHS[case]
    r, z = read_fort33(case_dir / "solps_output" / "fort.33")
    triangles = read_fort34(case_dir / "solps_output" / "fort.34")
    tri = mtri.Triangulation(r, z, triangles)
    vol = triangle_volume(r, z, triangles)
    wall_path = case_dir / "shapedata" / "wall.txt"
    wall = np.loadtxt(wall_path, skiprows=1) if wall_path.exists() else None
    return tri, wall, vol


def positive_lognorm(values: list[np.ndarray], lower_pct: float = 2.0, upper_pct: float = 99.5) -> mcolors.Normalize:
    merged = np.concatenate([np.ravel(v) for v in values])
    positive = merged[np.isfinite(merged) & (merged > 0.0)]
    if positive.size == 0:
        return mcolors.Normalize(vmin=0.0, vmax=1.0)
    vmin = np.nanpercentile(positive, lower_pct)
    vmax = np.nanpercentile(positive, upper_pct)
    vmin = max(vmin, np.nanmin(positive))
    vmax = max(vmax, vmin * 10.0)
    return mcolors.LogNorm(vmin=vmin, vmax=vmax)


def signed_symlognorm(values: list[np.ndarray]) -> mcolors.SymLogNorm:
    merged = np.concatenate([np.ravel(v) for v in values])
    finite = merged[np.isfinite(merged)]
    if finite.size == 0:
        return mcolors.SymLogNorm(linthresh=1.0, vmin=-1.0, vmax=1.0)
    scale = np.nanpercentile(np.abs(finite), 99.5)
    scale = max(scale, np.nanmax(np.abs(finite)), 1.0)
    linthresh = max(np.nanpercentile(np.abs(finite[np.abs(finite) > 0.0]), 5.0), scale * 1e-4)
    return mcolors.SymLogNorm(linthresh=linthresh, vmin=-scale, vmax=scale)


def save_figure(fig: plt.Figure, png_path: Path) -> None:
    fig.savefig(png_path, dpi=240)
    fig.savefig(png_path.with_suffix(".pdf"))


def load_case_fields(root: Path, run_prefix: str, case: int, per_volume: bool) -> dict[str, np.ndarray | float]:
    data_dir = root / "Outputfile" / f"{run_prefix}_{case}" / "data"
    _, _, vol = load_geometry(root, case)

    sn = np.loadtxt(data_dir / "Sn_CX_D_0_Tri")
    se = np.loadtxt(data_dir / "SE_CX_D_0_Tri")
    smu = np.loadtxt(data_dir / "Smu_CX_D_0_Tri")

    if len(sn) != len(vol):
        raise ValueError(f"case {case}: source length {len(sn)} != triangle count {len(vol)}")

    divisor = vol if per_volume else np.ones_like(vol)
    divisor = np.where(divisor > 0.0, divisor, np.nan)
    return {
        "sn_raw": sn,
        "se_raw": se,
        "smu_raw": smu,
        "sn_plot": sn / divisor,
        "se_net_plot": se / divisor,
        "se_loss_plot": np.abs(np.minimum(se, 0.0)) / divisor,
        "smu_abs_plot": np.abs(smu) / divisor,
        "sn_sum": float(np.nansum(sn)),
        "se_net_sum": float(np.nansum(se)),
        "se_pos_sum": float(np.nansum(se[se > 0.0])),
        "se_loss_abs_sum": float(np.nansum(np.abs(se[se < 0.0]))),
        "smu_abs_sum": float(np.nansum(np.abs(smu))),
    }


def draw_case_grid(
    root: Path,
    results: dict[int, dict],
    field_name: str,
    title: str,
    cbar_label: str,
    out_path: Path,
    cmap: str,
    norm: mcolors.Normalize,
    selected: list[int] | None = None,
) -> None:
    cases = selected or sorted(results)
    ncols = 3 if len(cases) > 3 else len(cases)
    nrows = int(np.ceil(len(cases) / ncols))
    height_per_row = 5.2 if nrows == 1 else 4.2
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, height_per_row * nrows), squeeze=False)
    mappable = None
    for axis, case in zip(axes.flat, cases):
        tri, wall, _ = load_geometry(root, case)
        values = np.asarray(results[case][field_name], dtype=float)
        plot_values = np.ma.masked_invalid(values)
        if isinstance(norm, mcolors.LogNorm):
            plot_values = np.ma.masked_where(values <= 0.0, plot_values)
        mappable = axis.tripcolor(
            tri,
            facecolors=plot_values,
            cmap=cmap,
            norm=norm,
            edgecolors="none",
            rasterized=True,
        )
        if wall is not None:
            axis.plot(wall[:, 0], wall[:, 1], color="black", linewidth=0.35)
        axis.set_title(f"case {case}: n={CASE_DENSITIES[case]/1e19:.1f}")
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlim(1.25, 2.85)
        axis.set_ylim(-1.22, 1.18)
        axis.set_xlabel("R [m]")
        axis.set_ylabel("Z [m]")
        axis.grid(alpha=0.10)
    for axis in axes.flat[len(cases) :]:
        axis.set_axis_off()
    fig.suptitle(title, fontsize=15, y=0.985)
    top = 0.80 if nrows == 1 else 0.91
    bottom = 0.22 if nrows == 1 else 0.16
    cbar_bottom = 0.095 if nrows == 1 else 0.065
    fig.subplots_adjust(left=0.06, right=0.98, top=top, bottom=bottom, wspace=0.22, hspace=0.28)
    if mappable is not None:
        cax = fig.add_axes([0.17, cbar_bottom, 0.66, 0.025])
        cbar = fig.colorbar(mappable, cax=cax, orientation="horizontal")
        cbar.set_label(cbar_label)
    save_figure(fig, out_path)
    plt.close(fig)


def write_summary_csv(results: dict[int, dict], out_path: Path) -> None:
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "case",
                "density_m^-3",
                "Sn_CX_D0_sum",
                "SE_CX_D0_net_sum",
                "SE_CX_D0_positive_sum",
                "SE_CX_D0_loss_abs_sum",
                "abs_Smu_CX_D0_sum",
            ]
        )
        for case in sorted(results):
            row = results[case]
            writer.writerow(
                [
                    case,
                    f"{CASE_DENSITIES[case]:.6e}",
                    f"{row['sn_sum']:.8e}",
                    f"{row['se_net_sum']:.8e}",
                    f"{row['se_pos_sum']:.8e}",
                    f"{row['se_loss_abs_sum']:.8e}",
                    f"{row['smu_abs_sum']:.8e}",
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", help="Repository root")
    parser.add_argument("--run-prefix", default=DEFAULT_RUN_PREFIX)
    parser.add_argument("--out-dir", default=None, help="Default: ROOT/Outputfile/wallscan_summary")
    parser.add_argument("--integrated", action="store_true", help="Plot triangle-integrated source instead of source/volume")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    out_dir = Path(args.out_dir) if args.out_dir else root / "Outputfile" / "wallscan_summary"
    out_dir.mkdir(parents=True, exist_ok=True)
    per_volume = not args.integrated

    results = {case: load_case_fields(root, args.run_prefix, case, per_volume=per_volume) for case in CASE_DENSITIES}
    suffix = "density" if per_volume else "integrated"
    unit_suffix = " / volume" if per_volume else ""

    sn_norm = positive_lognorm([np.asarray(results[c]["sn_plot"]) for c in results])
    se_loss_norm = positive_lognorm([np.asarray(results[c]["se_loss_plot"]) for c in results])
    smu_norm = positive_lognorm([np.asarray(results[c]["smu_abs_plot"]) for c in results])
    se_net_norm = signed_symlognorm([np.asarray(results[c]["se_net_plot"]) for c in results])

    draw_case_grid(
        root,
        results,
        "sn_plot",
        "CX particle source term on triangle mesh: Sn_CX_D_0_Tri",
        f"CX particle source{unit_suffix} (shared log scale)",
        out_dir / f"physics_cx_particle_source_{suffix}_all_cases.png",
        "viridis",
        sn_norm,
    )
    draw_case_grid(
        root,
        results,
        "se_net_plot",
        "Net CX energy source term on triangle mesh: SE_CX_D_0_Tri",
        f"net CX energy source{unit_suffix} (shared symmetric log scale)",
        out_dir / f"physics_cx_net_energy_source_{suffix}_all_cases.png",
        "coolwarm",
        se_net_norm,
    )
    draw_case_grid(
        root,
        results,
        "se_loss_plot",
        "CX energy-loss component on triangle mesh: abs(min(SE_CX_D_0_Tri, 0))",
        f"CX energy-loss component{unit_suffix} (shared log scale)",
        out_dir / f"physics_cx_energy_loss_source_{suffix}_all_cases.png",
        "magma",
        se_loss_norm,
    )
    draw_case_grid(
        root,
        results,
        "smu_abs_plot",
        "CX momentum source magnitude on triangle mesh: abs(Smu_CX_D_0_Tri)",
        f"CX momentum source magnitude{unit_suffix} (shared log scale)",
        out_dir / f"physics_cx_momentum_source_{suffix}_all_cases.png",
        "cividis",
        smu_norm,
    )

    selected = [1, 5, 9]
    draw_case_grid(
        root,
        results,
        "sn_plot",
        "Selected density cases: CX particle source",
        f"CX particle source{unit_suffix} (shared log scale)",
        out_dir / f"physics_cx_particle_source_{suffix}_selected.png",
        "viridis",
        sn_norm,
        selected=selected,
    )
    draw_case_grid(
        root,
        results,
        "se_net_plot",
        "Selected density cases: net CX energy source",
        f"net CX energy source{unit_suffix} (shared symmetric log scale)",
        out_dir / f"physics_cx_net_energy_source_{suffix}_selected.png",
        "coolwarm",
        se_net_norm,
        selected=selected,
    )
    draw_case_grid(
        root,
        results,
        "se_loss_plot",
        "Selected density cases: CX energy-loss component",
        f"CX energy-loss component{unit_suffix} (shared log scale)",
        out_dir / f"physics_cx_energy_loss_source_{suffix}_selected.png",
        "magma",
        se_loss_norm,
        selected=selected,
    )

    write_summary_csv(results, out_dir / f"physics_cx_source_summary_{suffix}.csv")
    print(f"Wrote CX source maps to {out_dir}")


if __name__ == "__main__":
    main()
