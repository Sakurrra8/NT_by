#!/usr/bin/env python3
"""Make compact physics-story plots for the EAST wall density scan."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("Outputfile") / "wallscan_summary" / ".mplconfig"))

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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

RUN_PREFIX = "wallscan_recyc500k_D"


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path)
    fig.savefig(path.with_suffix(".pdf"))


def segment_areas(segments: np.ndarray) -> np.ndarray:
    r1 = segments[:, 0, 0]
    z1 = segments[:, 0, 1]
    r2 = segments[:, 1, 0]
    z2 = segments[:, 1, 1]
    length = np.hypot(r2 - r1, z2 - z1)
    return np.pi * (r1 + r2) * length


def load_case(root: Path, case: int) -> dict[str, np.ndarray | float]:
    h5_path = root / "Outputfile" / f"{RUN_PREFIX}_{case}" / "data" / "D_W_stats.h5"
    with h5py.File(h5_path, "r") as h5:
        g = h5["D_W"]
        total_flux = np.asarray(g["totalFlux"], dtype=float)
        total_sputter = np.asarray(g["totalSputter"], dtype=float)
        energy_hist = np.asarray(g["energyHist"], dtype=float)
        angle_hist = np.asarray(g["angleHist"], dtype=float)
        sputter_energy_hist = np.asarray(g["sputterEnergyHist"], dtype=float)
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

    e_sum = energy_hist.sum(axis=0)
    a_sum = angle_hist.sum(axis=0)
    flux_sum = float(total_flux.sum())
    sputter_sum = float(total_sputter.sum())
    mean_e = float((e_sum * e_centers).sum() / e_sum.sum())
    mean_mu = float((a_sum * mu_centers).sum() / a_sum.sum())

    return {
        "density": CASE_DENSITIES[case],
        "total_flux": total_flux,
        "total_sputter": total_sputter,
        "energy_hist": energy_hist,
        "sputter_energy_hist": sputter_energy_hist,
        "angle_hist": angle_hist,
        "sputter_angle_hist": sputter_angle_hist,
        "e_centers": e_centers,
        "mu_centers": mu_centers,
        "flux_sum": flux_sum,
        "sputter_sum": sputter_sum,
        "effective_yield": sputter_sum / flux_sum,
        "mean_e": mean_e,
        "mean_mu": mean_mu,
        "frac_mu_lt_03": float(a_sum[mu_centers < 0.3].sum() / a_sum.sum()),
        "frac_mu_gt_07": float(a_sum[mu_centers > 0.7].sum() / a_sum.sum()),
        "frac_e50": float(e_sum[e_centers >= 50.0].sum() / e_sum.sum()),
        "frac_e178": float(e_sum[e_centers >= 178.0].sum() / e_sum.sum()),
        "frac_e200": float(e_sum[e_centers >= 200.0].sum() / e_sum.sum()),
        "frac_e400": float(e_sum[e_centers >= 400.0].sum() / e_sum.sum()),
        "physical_wall_count": physical_wall_count,
        "target_count": target_count,
    }


def load_wall_points(root: Path, case: int = 1) -> np.ndarray:
    wall_path = root / "case_input" / CASE_PATHS[case] / "shapedata" / "wall_segment.txt"
    rows = wall_path.read_text().splitlines()
    n_points = int(rows[0].strip())
    points = np.array([[float(x) for x in line.split()[:2]] for line in rows[1 : 1 + n_points]])
    return points


def load_target_segments(root: Path, case: int = 1) -> np.ndarray:
    shape_dir = root / "case_input" / CASE_PATHS[case] / "shapedata"

    def arr(name: str) -> np.ndarray:
        return np.loadtxt(shape_dir / name)

    x0, y0 = arr("com_nx1"), arr("com_ny1")
    x1, y1 = arr("com_nx2"), arr("com_ny2")
    x2, y2 = arr("com_nx3"), arr("com_ny3")
    x3, y3 = arr("com_nx4"), arr("com_ny4")
    n_poloidal, n_radial = x0.shape

    segments = []
    p_inner = 1
    p_outer = n_poloidal - 2
    for i in range(n_radial):
        segments.append([[x3[p_inner, i], y3[p_inner, i]], [x0[p_inner, i], y0[p_inner, i]]])
    for i in range(n_radial):
        segments.append([[x1[p_outer, i], y1[p_outer, i]], [x2[p_outer, i], y2[p_outer, i]]])
    return np.asarray(segments, dtype=float)


def load_surface_geometry(root: Path, case: int = 1) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    points = load_wall_points(root, case)
    wall_segments = np.stack([points[:-1], points[1:]], axis=1)
    target_segments = load_target_segments(root, case)
    segments = np.concatenate([wall_segments, target_segments], axis=0)
    areas = segment_areas(segments)
    kinds = np.array(["wall"] * len(wall_segments) + ["target"] * len(target_segments))
    return segments, areas, kinds


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8.5,
            "figure.dpi": 130,
            "savefig.dpi": 220,
        }
    )


def positive_or_nan(values: np.ndarray) -> np.ndarray:
    out = np.asarray(values, dtype=float).copy()
    out[out <= 0.0] = np.nan
    return out


def plot_density_story(results: dict[int, dict], out_dir: Path) -> None:
    cases = sorted(results)
    density = np.array([results[c]["density"] for c in cases]) / 1e19
    flux = np.array([results[c]["flux_sum"] for c in cases])
    sputter = np.array([results[c]["sputter_sum"] for c in cases])
    y_eff = np.array([results[c]["effective_yield"] for c in cases])
    mean_e = np.array([results[c]["mean_e"] for c in cases])
    frac200 = np.array([results[c]["frac_e200"] for c in cases])

    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.5), sharex=True)
    panels = [
        (axes[0, 0], flux, "Total neutral wall flux", "weighted s$^{-1}$", "C0"),
        (axes[0, 1], sputter, "Total W sputtering source", "weighted s$^{-1}$", "C3"),
        (axes[1, 0], y_eff, "Effective sputtering yield", "$S_W/\\Gamma_\\mathrm{wall}$", "C2"),
        (axes[1, 1], frac200, "High-energy fraction", "$E \\geq 200$ eV", "C4"),
    ]
    for ax, y, title, ylabel, color in panels:
        ax.plot(density, y, marker="o", color=color, lw=1.8)
        ax.set_yscale("log")
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25, which="both")
    axes[1, 0].set_xlabel("Background density label [$10^{19}$ m$^{-3}$]")
    axes[1, 1].set_xlabel("Background density label [$10^{19}$ m$^{-3}$]")
    fig.suptitle("Density scan: neutral loading increases while W sputtering weakens")
    fig.tight_layout()
    save_figure(fig, out_dir / "physics_density_scan_main.png")
    plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(7.6, 4.6))
    ax1.plot(density, mean_e, marker="o", color="C1", label="mean impact energy")
    ax1.set_xlabel("Background density label [$10^{19}$ m$^{-3}$]")
    ax1.set_ylabel("Mean impact energy [eV]", color="C1")
    ax1.tick_params(axis="y", labelcolor="C1")
    ax1.grid(alpha=0.25)
    ax2 = ax1.twinx()
    ax2.plot(density, frac200 * 100.0, marker="s", color="C4", label="$E \\geq 200$ eV")
    ax2.set_yscale("log")
    ax2.set_ylabel("Flux fraction with $E \\geq 200$ eV [%]", color="C4")
    ax2.tick_params(axis="y", labelcolor="C4")
    fig.suptitle("Density-driven spectral softening of wall-reaching neutrals")
    fig.tight_layout()
    save_figure(fig, out_dir / "physics_spectral_softening.png")
    plt.close(fig)


def plot_profiles_and_spectra(results: dict[int, dict], out_dir: Path) -> None:
    selected = [1, 3, 5, 7, 9]
    colors = plt.cm.viridis(np.linspace(0.08, 0.9, len(selected)))

    fig, axes = plt.subplots(2, 1, figsize=(8.6, 6.4), sharex=True)
    wall_count = int(results[selected[0]]["physical_wall_count"])
    target_count = int(results[selected[0]]["target_count"])
    inner_target_count = target_count // 2
    inner_target_end = wall_count + inner_target_count - 0.5
    for color, case in zip(colors, selected):
        label = f"{results[case]['density'] / 1e19:.1f}"
        axes[0].plot(positive_or_nan(results[case]["flux_intensity"]), lw=1.25, color=color, label=label)
        axes[1].plot(positive_or_nan(results[case]["sputter_intensity"]), lw=1.25, color=color, label=label)
    axes[0].set_yscale("log")
    axes[1].set_yscale("log")
    axes[0].set_ylabel("Neutral flux intensity [m$^{-2}$ s$^{-1}$]")
    axes[1].set_ylabel("W sputtering intensity [m$^{-2}$ s$^{-1}$]")
    axes[1].set_xlabel("Surface segment (wall first, targets last)")
    axes[0].set_title("Surface-resolved neutral loading intensity")
    axes[1].set_title("Surface-resolved W sputtering intensity")
    for ax in axes:
        ax.axvspan(wall_count - 0.5, inner_target_end, color="C0", alpha=0.055)
        ax.axvspan(inner_target_end, wall_count + target_count - 0.5, color="C1", alpha=0.055)
        ax.axvline(wall_count - 0.5, color="0.35", ls="--", lw=1.0)
        ax.axvline(inner_target_end, color="0.45", ls=":", lw=1.0)
        ax.text(
            wall_count + 2,
            0.05,
            "inner target",
            transform=ax.get_xaxis_transform(),
            color="0.25",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.5},
        )
        ax.text(
            inner_target_end + 2,
            0.05,
            "outer target",
            transform=ax.get_xaxis_transform(),
            color="0.25",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.5},
        )
        ax.grid(alpha=0.22, which="both")
        ax.legend(title="$n$ [$10^{19}$ m$^{-3}$]", ncol=5, loc="upper right")
    fig.tight_layout()
    save_figure(fig, out_dir / "physics_surface_intensity_profiles_selected.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    for color, case in zip(colors, selected):
        e = results[case]["e_centers"]
        spec = results[case]["energy_hist"].sum(axis=0)
        spec = spec / spec.sum()
        ax.semilogy(e, spec, lw=1.7, color=color, label=f"{results[case]['density'] / 1e19:.1f}")
    ax.axvline(200.0, ls="--", lw=1.0, color="0.35", label="200 eV marker")
    ax.set_xlim(0, 500)
    ax.set_ylim(1e-6, None)
    ax.set_xlabel("Wall-impact energy [eV]")
    ax.set_ylabel("Normalized flux spectrum")
    ax.set_title("Integrated wall-impact energy spectra")
    ax.grid(alpha=0.25, which="both")
    ax.legend(title="$n$ [$10^{19}$ m$^{-3}$]")
    fig.tight_layout()
    save_figure(fig, out_dir / "physics_energy_spectra_selected.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    for color, case in zip(colors, selected):
        e = results[case]["e_centers"]
        spec = results[case]["energy_hist"].sum(axis=0)
        tail = np.cumsum(spec[::-1])[::-1] / spec.sum()
        ax.semilogy(e, tail, lw=1.9, color=color, label=f"{results[case]['density'] / 1e19:.1f}")
    for marker in [50.0, 200.0]:
        ax.axvline(marker, ls="--", lw=0.9, color="0.45")
    ax.set_xlim(0, 500)
    ax.set_ylim(1e-4, 1.1)
    ax.set_xlabel("Energy threshold $E_\\mathrm{th}$ [eV]")
    ax.set_ylabel("Flux fraction with $E \\geq E_\\mathrm{th}$")
    ax.set_title("Cumulative high-energy tail of wall-impact neutrals")
    ax.grid(alpha=0.25, which="both")
    ax.legend(title="$n$ [$10^{19}$ m$^{-3}$]")
    fig.tight_layout()
    save_figure(fig, out_dir / "physics_cumulative_high_energy_tail.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    for color, case in zip(colors, selected):
        e = results[case]["e_centers"]
        spec = results[case]["sputter_energy_hist"].sum(axis=0)
        if spec.sum() > 0:
            spec = spec / spec.sum()
        ax.semilogy(e, spec, lw=1.7, color=color, label=f"{results[case]['density'] / 1e19:.1f}")
    ax.set_xlim(0, 1500)
    ax.set_ylim(1e-6, None)
    ax.set_xlabel("Incident D energy [eV]")
    ax.set_ylabel("Normalized sputtering contribution")
    ax.set_title("Energy channels contributing to W sputtering")
    ax.grid(alpha=0.25, which="both")
    ax.legend(title="$n$ [$10^{19}$ m$^{-3}$]")
    fig.tight_layout()
    save_figure(fig, out_dir / "physics_sputter_energy_channels.png")
    plt.close(fig)


def plot_direction_story(results: dict[int, dict], out_dir: Path) -> None:
    selected = [1, 3, 5, 7, 9]
    colors = plt.cm.viridis(np.linspace(0.08, 0.9, len(selected)))

    fig, axes = plt.subplots(2, 1, figsize=(7.8, 6.2), sharex=True)
    for color, case in zip(colors, selected):
        mu = results[case]["mu_centers"]
        flux_angle = results[case]["angle_hist"].sum(axis=0)
        sput_angle = results[case]["sputter_angle_hist"].sum(axis=0)
        if flux_angle.sum() > 0:
            flux_angle = flux_angle / flux_angle.sum()
        if sput_angle.sum() > 0:
            sput_angle = sput_angle / sput_angle.sum()
        label = f"{results[case]['density'] / 1e19:.1f}"
        axes[0].plot(mu, flux_angle, lw=1.7, color=color, label=label)
        axes[1].plot(mu, sput_angle, lw=1.7, color=color, label=label)
    axes[0].set_ylabel("Normalized neutral flux")
    axes[1].set_ylabel("Normalized sputtering source")
    axes[1].set_xlabel(r"Incident $\mu = \cos(\theta)$")
    axes[0].set_title("Wall-impact direction distribution")
    axes[1].set_title("Direction channels contributing to W sputtering")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(title="$n$ [$10^{19}$ m$^{-3}$]", ncol=5, loc="upper center")
    fig.tight_layout()
    save_figure(fig, out_dir / "physics_incident_mu_distributions.png")
    plt.close(fig)

    cases = sorted(results)
    density = np.array([results[c]["density"] for c in cases]) / 1e19
    mean_mu = np.array([results[c]["mean_mu"] for c in cases])
    frac_grazing = np.array([results[c]["frac_mu_lt_03"] for c in cases])
    frac_normal = np.array([results[c]["frac_mu_gt_07"] for c in cases])

    fig, ax1 = plt.subplots(figsize=(7.6, 4.6))
    ax1.plot(density, mean_mu, marker="o", color="C0", label=r"mean $\mu$")
    ax1.set_xlabel("Background density label [$10^{19}$ m$^{-3}$]")
    ax1.set_ylabel(r"Mean incident $\mu$", color="C0")
    ax1.tick_params(axis="y", labelcolor="C0")
    ax1.grid(alpha=0.25)
    ax2 = ax1.twinx()
    ax2.plot(density, frac_grazing * 100.0, marker="s", color="C2", label=r"$\mu < 0.3$")
    ax2.plot(density, frac_normal * 100.0, marker="^", color="C3", label=r"$\mu > 0.7$")
    ax2.set_ylabel("Flux fraction [%]", color="0.2")
    ax2.tick_params(axis="y", labelcolor="0.2")
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines], loc="best")
    fig.suptitle("Density dependence of wall-impact direction")
    fig.tight_layout()
    save_figure(fig, out_dir / "physics_incident_mu_trend.png")
    plt.close(fig)


def plot_wall_location_maps(root: Path, results: dict[int, dict], out_dir: Path) -> None:
    segments, _, kinds = load_surface_geometry(root, 1)
    selected = [1, 5, 9]
    target_mask = kinds == "target"
    target_segments = segments[target_mask]
    inner_target_segments = target_segments[: len(target_segments) // 2]
    outer_target_segments = target_segments[len(target_segments) // 2 :]
    inner_target_center = inner_target_segments.mean(axis=(0, 1))
    outer_target_center = outer_target_segments.mean(axis=(0, 1))

    def draw_segment_map(values_name: str, title: str, out_name: str, log_scale: bool = True) -> None:
        fig, axes = plt.subplots(1, len(selected), figsize=(12, 4.8), sharex=True, sharey=True)
        all_vals = np.concatenate([np.asarray(results[c][values_name], dtype=float) for c in selected])
        positive = all_vals[all_vals > 0]
        vmin = np.percentile(positive, 5) if positive.size else 1.0
        vmax = np.percentile(positive, 99) if positive.size else 1.0
        norm = matplotlib.colors.LogNorm(vmin=max(vmin, 1e-300), vmax=max(vmax, vmin * 1.01)) if log_scale else None
        cmap = plt.get_cmap("inferno")
        for ax, case in zip(axes, selected):
            vals = np.asarray(results[case][values_name], dtype=float)
            for idx, seg in enumerate(segments):
                color = cmap(norm(vals[idx])) if norm is not None and vals[idx] > 0 else "0.85"
                lw = 2.0 if kinds[idx] == "wall" else 4.2
                ax.plot(seg[:, 0], seg[:, 1], color=color, lw=lw, solid_capstyle="round")
            mids = segments[target_mask].mean(axis=1)
            target_vals = vals[target_mask]
            target_positive = target_vals > 0
            if np.any(target_positive):
                marker_colors = [cmap(norm(v)) for v in target_vals[target_positive]]
                ax.scatter(
                    mids[target_positive, 0],
                    mids[target_positive, 1],
                    c=marker_colors,
                    s=18,
                    edgecolors="black",
                    linewidths=0.25,
                    zorder=5,
                )
            ax.annotate(
                "inner target",
                xy=inner_target_center,
                xytext=(inner_target_center[0] - 0.32, inner_target_center[1] + 0.18),
                arrowprops={"arrowstyle": "->", "lw": 0.8, "color": "0.2"},
                fontsize=8,
                color="0.2",
            )
            ax.annotate(
                "outer target",
                xy=outer_target_center,
                xytext=(outer_target_center[0] + 0.08, outer_target_center[1] + 0.18),
                arrowprops={"arrowstyle": "->", "lw": 0.8, "color": "0.2"},
                fontsize=8,
                color="0.2",
            )
            ax.set_aspect("equal", adjustable="box")
            ax.set_title(f"n = {results[case]['density'] / 1e19:.1f}")
            ax.set_xlabel("R [m]")
            ax.grid(alpha=0.15)
        axes[0].set_ylabel("Z [m]")
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        fig.colorbar(sm, ax=axes, orientation="horizontal", fraction=0.055, pad=0.12, aspect=45)
        fig.suptitle(title)
        fig.subplots_adjust(left=0.06, right=0.98, top=0.86, bottom=0.2, wspace=0.2)
        save_figure(fig, out_dir / out_name)
        plt.close(fig)

    high_flux_by_wall = {}
    for case, r in results.items():
        e = r["e_centers"]
        high_flux_by_wall[case] = r["energy_hist"][:, e >= 200.0].sum(axis=1)
        r["high_flux_intensity_e200"] = high_flux_by_wall[case] / r["surface_area"]

    draw_segment_map(
        "high_flux_intensity_e200",
        "Surface locations reached by high-energy neutrals (E >= 200 eV)",
        "physics_wall_location_high_energy_flux.png",
    )
    draw_segment_map(
        "sputter_intensity",
        "Surface locations contributing to W sputtering intensity",
        "physics_wall_location_sputter_source.png",
    )


def write_physics_note(results: dict[int, dict], out_dir: Path) -> None:
    r1 = results[1]
    r9 = results[9]
    flux_ratio = r9["flux_sum"] / r1["flux_sum"]
    sputter_ratio = r1["sputter_sum"] / r9["sputter_sum"]
    yield_ratio = r1["effective_yield"] / r9["effective_yield"]
    frac200_ratio = r1["frac_e200"] / r9["frac_e200"]

    note = f"""# EAST 高能中性粒子第一壁加载与潜在侵蚀：密度扫描小结

## 核心观察

这组密度扫描里，总的中性粒子壁面加载随密度升高而增加，但用 W 作为高阈值材料诊断时，物理溅射源反而下降。这里最关键的控制量不是总中性粒子通量，而是到达第一壁的能量谱，尤其是高能尾。

从 case 1 到 case 9：

- 背景密度标签：`3.2e19 -> 5.8e19 m^-3`
- 总中性粒子壁面通量：`{r1['flux_sum']:.3e} -> {r9['flux_sum']:.3e}`，增加 `{flux_ratio:.2f}` 倍
- 总 W 溅射源：`{r1['sputter_sum']:.3e} -> {r9['sputter_sum']:.3e}`，降低 `{sputter_ratio:.1f}` 倍
- 有效溅射产额：`{r1['effective_yield']:.3e} -> {r9['effective_yield']:.3e}`，降低 `{yield_ratio:.1f}` 倍
- 平均入射能量：`{r1['mean_e']:.2f} -> {r9['mean_e']:.2f} eV`
- `E >= 200 eV` 的通量比例：`{100*r1['frac_e200']:.2f}% -> {100*r9['frac_e200']:.3f}%`，降低 `{frac200_ratio:.1f}` 倍

## 物理解释

密度升高后，中性粒子总加载增强，但到达第一壁的粒子能谱明显软化。对 W 这类高阈值材料而言，高能尾的削弱比总通量增加更重要，因此会出现“总中性流增加、W 溅射反而减小”的趋势。

可采用的机制链条是：

```text
边缘/背景密度升高
  -> CX 中性粒子和背景等离子体的碰撞、再电离、能量退化增强
  -> 高能中性粒子的有效平均自由程变短
  -> 能够保持高能量到达第一壁的粒子减少
  -> 到壁能谱软化，高能尾减弱
  -> W 的有效物理溅射产额下降
  -> 中性粒子诱导的 W 溅射源降低
```

## 本文建议的边界

本文重点建议定为：非靶板第一壁区域的高能中性粒子加载与潜在侵蚀。靶板区域仍然保留在图中，用于说明几何边界和粒子去向；但靶板处通常有强离子流和热流主导侵蚀，因此不把靶板中性粒子侵蚀作为主要结论。

W 在这里更适合作为高阈值壁材料的溅射诊断，不建议把文章写成完整 W 杂质输运文章。后续如果需要，可以只加入轻量的“新产生 W 的第一次电离位置”诊断，作为补充而不是主线。

## 主图

### 1. 密度扫描总趋势

![密度扫描总趋势](./physics_density_scan_main.png)

### 2. 能谱软化：平均能量与高能比例

![能谱软化](./physics_spectral_softening.png)

### 3. 面积归一后的表面加载与潜在溅射强度

图中右侧浅蓝区为 inner target，浅橙区为 outer target。靶板处保留用于说明边界和粒子去向，但不作为本文侵蚀结论重点。

![表面强度分布](./physics_surface_intensity_profiles_selected.png)

### 4. 到壁中性粒子能谱

![能谱](./physics_energy_spectra_selected.png)

### 5. 累计高能尾

![累计高能尾](./physics_cumulative_high_energy_tail.png)

### 6. 贡献 W 溅射的能量通道

![溅射能量通道](./physics_sputter_energy_channels.png)

### 7. 到壁入射方向分布

这里的 `mu` 来自 wall-impact tracker，代码中用 `theta = acos(mu)` 计算角度相关的溅射产额。因此 `mu` 越接近 1，越接近法向入射；`mu` 越接近 0，越接近掠入射。

![到壁入射方向分布](./physics_incident_mu_distributions.png)

### 8. 入射方向随密度的变化

![入射方向趋势](./physics_incident_mu_trend.png)

### 9. 高能中性粒子到达位置

![高能中性粒子空间分布](./physics_wall_location_high_energy_flux.png)

### 10. W 溅射源空间分布

![W 溅射源空间分布](./physics_wall_location_sputter_source.png)

## 还需要核对

1. 量化所有 case 的 lost/stalled particle weights。
2. 解释高密度 case 中 `Ep_2D.dat` 和 `Er_2D.dat` 的读取 warning，尤其是 ExB drift 关闭时是否影响结论。
3. 如果要强调角度效应，后续最好导出真正的联合 `E-mu` 直方图。
4. 论文级表述前，建议直接从程序输出 wall/target segment 的几何面积或靶板面积，避免 Python 重构靶板面积带来的不确定性。
"""
    (out_dir / "physics_summary_note.md").write_text(note, encoding="utf-8")


def main() -> None:
    global RUN_PREFIX

    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", help="Repository root")
    parser.add_argument("--run-prefix", default=RUN_PREFIX, help="Outputfile run directory prefix")
    args = parser.parse_args()

    RUN_PREFIX = args.run_prefix

    root = Path(args.root)
    out_dir = root / "Outputfile" / "wallscan_summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    set_style()
    results = {case: load_case(root, case) for case in CASE_DENSITIES}
    _, surface_area, _ = load_surface_geometry(root, 1)
    for result in results.values():
        if len(surface_area) != len(result["total_flux"]):
            raise SystemExit(
                f"Surface geometry count {len(surface_area)} does not match H5 count {len(result['total_flux'])}"
            )
        result["surface_area"] = surface_area
        result["flux_intensity"] = result["total_flux"] / surface_area
        result["sputter_intensity"] = result["total_sputter"] / surface_area
    plot_density_story(results, out_dir)
    plot_profiles_and_spectra(results, out_dir)
    plot_direction_story(results, out_dir)
    plot_wall_location_maps(root, results, out_dir)
    write_physics_note(results, out_dir)
    print(f"Wrote physics-story plots and note to {out_dir}")


if __name__ == "__main__":
    main()
