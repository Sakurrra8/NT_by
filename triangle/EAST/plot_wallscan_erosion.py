#!/usr/bin/env python3
"""Summarize and plot neutral-induced first-wall erosion density scans.

Expected input layout:
  Outputfile/wallscan_D_1/data/D_W_stats.h5
  ...
  Outputfile/wallscan_D_9/data/D_W_stats.h5

The script writes:
  Outputfile/wallscan_summary/wallscan_erosion_summary.csv
  Outputfile/wallscan_summary/*.png
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("Outputfile") / "wallscan_summary" / ".mplconfig"))

try:
    import h5py
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError as exc:  # pragma: no cover - user-facing dependency gate
    raise SystemExit(
        "Missing Python dependency. In WSL, run:\n"
        "  python3 -m pip install h5py matplotlib numpy\n"
        f"Original error: {exc}"
    )


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

DEFAULT_RUN_PREFIX = "wallscan_recyc500k_D"


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path)
    fig.savefig(path.with_suffix(".pdf"))


def read_case(path: Path) -> dict[str, np.ndarray | float | int]:
    with h5py.File(path, "r") as h5:
        g = h5["D_W"]
        total_flux = np.asarray(g["totalFlux"], dtype=float)
        total_sputter = np.asarray(g["totalSputter"], dtype=float)
        energy_hist = np.asarray(g["energyHist"], dtype=float)
        sputter_energy_hist = np.asarray(g["sputterEnergyHist"], dtype=float)
        angle_hist = np.asarray(g["angleHist"], dtype=float)
        sputter_angle_hist = np.asarray(g["sputterAngleHist"], dtype=float)

        e_min = float(g.attrs["E_min"])
        e_max = float(g.attrs["E_max"])
        mu_min = float(g.attrs["mu_min"])
        mu_max = float(g.attrs["mu_max"])
        physical_wall_count = int(g.attrs.get("physical_wall_count", len(total_flux)))
        target_count = int(g.attrs.get("target_count", 0))

    e_centers = np.linspace(e_min, e_max, energy_hist.shape[1], endpoint=False)
    e_centers += 0.5 * (e_max - e_min) / energy_hist.shape[1]
    mu_centers = np.linspace(mu_min, mu_max, angle_hist.shape[1], endpoint=False)
    mu_centers += 0.5 * (mu_max - mu_min) / angle_hist.shape[1]

    flux_sum = float(total_flux.sum())
    sputter_sum = float(total_sputter.sum())
    peak_wall = int(total_sputter.argmax()) if total_sputter.size else -1
    peak_sputter = float(total_sputter[peak_wall]) if peak_wall >= 0 else 0.0
    y_eff = sputter_sum / flux_sum if flux_sum > 0 else 0.0

    energy_sum = energy_hist.sum(axis=0)
    e_mean = float((energy_sum * e_centers).sum() / energy_sum.sum()) if energy_sum.sum() > 0 else 0.0
    high_e_flux = float(energy_sum[e_centers >= 50.0].sum())
    high_e_fraction = high_e_flux / flux_sum if flux_sum > 0 else 0.0
    high_e_fraction_178 = float(energy_sum[e_centers >= 178.0].sum() / flux_sum) if flux_sum > 0 else 0.0
    high_e_fraction_200 = float(energy_sum[e_centers >= 200.0].sum() / flux_sum) if flux_sum > 0 else 0.0
    high_e_fraction_400 = float(energy_sum[e_centers >= 400.0].sum() / flux_sum) if flux_sum > 0 else 0.0

    angle_sum = angle_hist.sum(axis=0)
    mu_mean = float((angle_sum * mu_centers).sum() / angle_sum.sum()) if angle_sum.sum() > 0 else 0.0

    return {
        "totalFlux": total_flux,
        "totalSputter": total_sputter,
        "energyHist": energy_hist,
        "sputterEnergyHist": sputter_energy_hist,
        "angleHist": angle_hist,
        "sputterAngleHist": sputter_angle_hist,
        "eCenters": e_centers,
        "muCenters": mu_centers,
        "flux_sum": flux_sum,
        "sputter_sum": sputter_sum,
        "peak_wall": peak_wall,
        "peak_sputter": peak_sputter,
        "y_eff": y_eff,
        "e_mean": e_mean,
        "high_e_fraction": high_e_fraction,
        "high_e_fraction_178": high_e_fraction_178,
        "high_e_fraction_200": high_e_fraction_200,
        "high_e_fraction_400": high_e_fraction_400,
        "mu_mean": mu_mean,
        "physical_wall_count": physical_wall_count,
        "target_count": target_count,
    }


def plot_wall_profiles(results: dict[int, dict], out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.8))
    for case, result in results.items():
        density = CASE_DENSITIES[case] / 1e19
        ax.plot(result["totalSputter"], label=f"{density:.1f}")
    ax.set_xlabel("Surface segment (wall first, target last)")
    ax.set_ylabel("Sputtered W source [weighted particles]")
    ax.set_title("First-wall neutral-induced sputtering profile")
    ax.legend(title=r"$n_e$ scan")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    save_figure(fig, out_dir / "wall_sputter_profiles.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4.8))
    for case, result in results.items():
        density = CASE_DENSITIES[case] / 1e19
        ax.plot(result["totalFlux"], label=f"{density:.1f}")
    ax.set_xlabel("Surface segment (wall first, target last)")
    ax.set_ylabel("Neutral wall flux [weighted particles]")
    ax.set_title("First-wall neutral impact profile")
    ax.legend(title=r"$n_e$ scan")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    save_figure(fig, out_dir / "wall_flux_profiles.png")
    plt.close(fig)


def plot_scan_scalars(results: dict[int, dict], out_dir: Path) -> None:
    cases = sorted(results)
    density = np.array([CASE_DENSITIES[c] for c in cases]) / 1e19
    flux = np.array([results[c]["flux_sum"] for c in cases])
    sputter = np.array([results[c]["sputter_sum"] for c in cases])
    y_eff = np.array([results[c]["y_eff"] for c in cases])
    high_e = np.array([results[c]["high_e_fraction"] for c in cases])

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
    panels = [
        (axes[0, 0], flux, "Total wall flux"),
        (axes[0, 1], sputter, "Total sputtering source"),
        (axes[1, 0], y_eff, "Effective yield"),
        (axes[1, 1], high_e, "Flux fraction E >= 50 eV"),
    ]
    for ax, y, label in panels:
        ax.plot(density, y, marker="o")
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
    axes[1, 0].set_xlabel(r"Case density [$10^{19}$ m$^{-3}$]")
    axes[1, 1].set_xlabel(r"Case density [$10^{19}$ m$^{-3}$]")
    fig.suptitle("Density scan summary")
    fig.tight_layout()
    save_figure(fig, out_dir / "density_scan_summary.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    series = [
        (flux / flux[0], "wall flux / case 1"),
        (sputter / sputter[0], "sputtering / case 1"),
        (y_eff / y_eff[0], "effective yield / case 1"),
        (high_e / high_e[0], "E >= 50 eV fraction / case 1"),
    ]
    for y, label in series:
        ax.plot(density, y, marker="o", label=label)
    ax.set_yscale("log")
    ax.set_xlabel(r"Case density [$10^{19}$ m$^{-3}$]")
    ax.set_ylabel("Normalized value")
    ax.set_title("Flux-erosion decoupling across density scan")
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    fig.tight_layout()
    save_figure(fig, out_dir / "density_scan_normalized_trends.png")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
    e_mean = np.array([results[c]["e_mean"] for c in cases])
    high_178 = np.array([results[c]["high_e_fraction_178"] for c in cases])
    high_200 = np.array([results[c]["high_e_fraction_200"] for c in cases])
    high_400 = np.array([results[c]["high_e_fraction_400"] for c in cases])
    mu_mean = np.array([results[c]["mu_mean"] for c in cases])
    panels = [
        (axes[0, 0], e_mean, "Mean impact energy [eV]"),
        (axes[0, 1], high_178, "Flux fraction E >= 178 eV"),
        (axes[1, 0], high_200, "Flux fraction E >= 200 eV"),
        (axes[1, 1], high_400, "Flux fraction E >= 400 eV"),
    ]
    for ax, y, label in panels:
        ax.plot(density, y, marker="o")
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
    axes[1, 0].set_xlabel(r"Case density [$10^{19}$ m$^{-3}$]")
    axes[1, 1].set_xlabel(r"Case density [$10^{19}$ m$^{-3}$]")
    fig.suptitle("Density-driven spectral softening")
    fig.tight_layout()
    save_figure(fig, out_dir / "density_scan_spectral_softening.png")
    plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(8, 4.8))
    ax1.plot(density, flux / flux[0], marker="o", label="wall flux")
    ax1.plot(density, y_eff / y_eff[0], marker="o", label="effective yield")
    ax1.plot(density, sputter / sputter[0], marker="o", label="sputtering")
    ax1.plot(density, high_200 / high_200[0], marker="o", label="E >= 200 eV fraction")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"Case density [$10^{19}$ m$^{-3}$]")
    ax1.set_ylabel("Normalized value")
    ax1.set_title("Mechanism decomposition: flux times spectral yield")
    ax1.grid(alpha=0.25, which="both")
    ax1.legend()
    fig.tight_layout()
    save_figure(fig, out_dir / "density_scan_mechanism_decomposition.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(density, mu_mean, marker="o")
    ax.set_xlabel(r"Case density [$10^{19}$ m$^{-3}$]")
    ax.set_ylabel(r"Mean incident $\mu$")
    ax.set_title("Wall-impact angle trend")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    save_figure(fig, out_dir / "density_scan_angle_trend.png")
    plt.close(fig)


def plot_peak_spectra(results: dict[int, dict], out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for case, result in results.items():
        peak = result["peak_wall"]
        spec = result["energyHist"][peak]
        norm = spec.sum()
        if norm > 0:
            spec = spec / norm
        ax.plot(result["eCenters"], spec, label=f"{CASE_DENSITIES[case] / 1e19:.1f}")
    ax.set_xlim(0, 300)
    ax.set_xlabel("Impact energy [eV]")
    ax.set_ylabel("Normalized flux spectrum at peak erosion wall")
    ax.set_title("Peak-wall neutral impact energy spectra")
    ax.legend(title=r"$n_e$ scan")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    save_figure(fig, out_dir / "peak_wall_energy_spectra.png")
    plt.close(fig)


def write_summary(results: dict[int, dict], out_dir: Path) -> None:
    fields = [
        "case",
        "density_m-3",
        "total_wall_flux",
        "total_sputter_source",
        "peak_wall_segment",
        "peak_sputter_source",
        "effective_yield",
        "mean_impact_energy_eV",
        "fraction_E_ge_50eV",
        "fraction_E_ge_178eV",
        "fraction_E_ge_200eV",
        "fraction_E_ge_400eV",
        "mean_incident_mu",
        "physical_wall_count",
        "target_count",
    ]
    with (out_dir / "wallscan_erosion_summary.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for case in sorted(results):
            r = results[case]
            writer.writerow(
                {
                    "case": case,
                    "density_m-3": CASE_DENSITIES[case],
                    "total_wall_flux": r["flux_sum"],
                    "total_sputter_source": r["sputter_sum"],
                    "peak_wall_segment": r["peak_wall"],
                    "peak_sputter_source": r["peak_sputter"],
                    "effective_yield": r["y_eff"],
                    "mean_impact_energy_eV": r["e_mean"],
                    "fraction_E_ge_50eV": r["high_e_fraction"],
                    "fraction_E_ge_178eV": r["high_e_fraction_178"],
                    "fraction_E_ge_200eV": r["high_e_fraction_200"],
                    "fraction_E_ge_400eV": r["high_e_fraction_400"],
                    "mean_incident_mu": r["mu_mean"],
                    "physical_wall_count": r["physical_wall_count"],
                    "target_count": r["target_count"],
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", help="Repository root")
    parser.add_argument("--cases", nargs="*", type=int, default=list(range(1, 10)))
    parser.add_argument("--run-prefix", default=DEFAULT_RUN_PREFIX, help="Outputfile case directory prefix")
    parser.add_argument("--summary-name", default="wallscan_summary", help="Outputfile summary directory name")
    args = parser.parse_args()

    root = Path(args.root)
    out_dir = root / "Outputfile" / args.summary_name
    out_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    missing = []
    for case in args.cases:
        h5 = root / "Outputfile" / f"{args.run_prefix}_{case}" / "data" / "D_W_stats.h5"
        if not h5.exists():
            missing.append(str(h5))
            continue
        results[case] = read_case(h5)

    if not results:
        raise SystemExit("No D_W_stats.h5 files found. Missing examples:\n" + "\n".join(missing[:3]))

    write_summary(results, out_dir)
    plot_wall_profiles(results, out_dir)
    plot_scan_scalars(results, out_dir)
    plot_peak_spectra(results, out_dir)

    print(f"Wrote summary and plots to {out_dir}")
    if missing:
        print("Missing cases:")
        for path in missing:
            print(f"  {path}")


if __name__ == "__main__":
    main()
