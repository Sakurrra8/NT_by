#!/usr/bin/env python3
"""Create a shared-colorbar NT vs SOLPS/EIRENE triangle comparison PDF."""

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
from matplotlib.backends.backend_pdf import PdfPages


def read_fort33(path):
    values = Path(path).read_text().split()
    count = int(values[0])
    r = np.asarray(values[1 : 1 + count], dtype=float) / 100.0
    z = np.asarray(values[1 + count : 1 + 2 * count], dtype=float) / 100.0
    return r, z


def read_fort34(path):
    return np.loadtxt(path, dtype=int, skiprows=1)[:, 1:4] - 1


def read_fort35_b2(path):
    rows = np.loadtxt(path, dtype=int, skiprows=1)
    b2 = np.full((rows.shape[0], 2), -1, dtype=int)
    for index, row in enumerate(rows):
        poloidal, radial = row[-2], row[-1]
        if poloidal < 0 or radial < 0:
            continue
        if poloidal > 74:
            poloidal -= 2
        elif poloidal > 25:
            poloidal -= 1
        b2[index] = (poloidal, radial)
    return b2


def triangle_volume(r, z, triangles):
    rr = r[triangles]
    zz = z[triangles]
    area = 0.5 * np.abs(
        (rr[:, 1] - rr[:, 0]) * (zz[:, 2] - zz[:, 0])
        - (rr[:, 2] - rr[:, 0]) * (zz[:, 1] - zz[:, 0])
    )
    return 2.0 * np.pi * rr.mean(axis=1) * area


def read_b2_matrix(path):
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 2 and data.shape == (98, 38):
        return data
    if data.ndim == 2 and data.shape == (38, 98):
        return data.T
    if data.ndim == 2 and data.shape[1] >= 4:
        return data[:, 3].reshape(38, 98).T
    raise ValueError(f"Unrecognized B2 field shape for {path}: {data.shape}")


def map_b2_value_to_tri(field, b2):
    result = np.full(len(b2), np.nan)
    valid = (
        (b2[:, 0] >= 0)
        & (b2[:, 0] < field.shape[0])
        & (b2[:, 1] >= 0)
        & (b2[:, 1] < field.shape[1])
    )
    result[valid] = field[b2[valid, 0], b2[valid, 1]]
    return result


def distribute_b2_integral_to_tri(field, b2, tri_volume):
    """Distribute a B2 cell-integrated quantity using triangle volume fractions."""
    result = np.full(len(b2), np.nan)
    valid = (
        (b2[:, 0] >= 0)
        & (b2[:, 0] < field.shape[0])
        & (b2[:, 1] >= 0)
        & (b2[:, 1] < field.shape[1])
    )
    volume_sum = np.zeros(field.shape)
    np.add.at(volume_sum, (b2[valid, 0], b2[valid, 1]), tri_volume[valid])
    denominator = volume_sum[b2[valid, 0], b2[valid, 1]]
    ok = denominator > 0.0
    valid_indices = np.flatnonzero(valid)
    result[valid_indices[ok]] = (
        field[b2[valid_indices[ok], 0], b2[valid_indices[ok], 1]]
        * tri_volume[valid_indices[ok]]
        / denominator[ok]
    )
    return result


def combined_log_norm(first, second):
    values = np.concatenate((np.asarray(first), np.asarray(second)))
    positive = values[np.isfinite(values) & (values > 0.0)]
    if positive.size == 0:
        return mcolors.Normalize(vmin=0.0, vmax=1.0)
    vmin = np.nanpercentile(positive, 1.0)
    vmax = np.nanpercentile(positive, 99.5)
    vmin = max(vmin, np.nanmin(positive))
    vmax = max(vmax, vmin * 10.0)
    return mcolors.LogNorm(vmin=vmin, vmax=vmax)


def metrics(name, nt, reference):
    valid = np.isfinite(nt) & np.isfinite(reference)
    active = valid & ((nt != 0.0) | (reference != 0.0))
    positive = valid & (nt > 0.0) & (reference > 0.0)
    row = {"field": name, "active_triangles": int(active.sum())}
    if active.any():
        denominator = np.sum(np.abs(reference[active]))
        row["relative_l1"] = (
            np.sum(np.abs(nt[active] - reference[active])) / denominator
            if denominator > 0.0
            else np.nan
        )
        row["relative_bias"] = (
            np.sum(nt[active] - reference[active]) / denominator
            if denominator > 0.0
            else np.nan
        )
        row["nt_sum"] = float(np.sum(nt[active]))
        row["solps_sum"] = float(np.sum(reference[active]))
    if positive.any():
        ratio = nt[positive] / reference[positive]
        row["median_positive_ratio"] = float(np.median(ratio))
        row["log10_rmse"] = float(
            np.sqrt(np.mean(np.square(np.log10(ratio))))
        )
    return row


def add_comparison_page(
    pdf,
    triangulation,
    wall,
    name,
    unit,
    nt,
    reference,
    xlim,
    ylim,
):
    shared_norm = combined_log_norm(nt, reference)
    ratio = np.full_like(nt, np.nan, dtype=float)
    valid = np.isfinite(nt) & np.isfinite(reference) & (nt > 0.0) & (reference > 0.0)
    ratio[valid] = nt[valid] / reference[valid]

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(14.0, 5.4),
        gridspec_kw={"width_ratios": [1.0, 1.0, 1.0]},
        constrained_layout=True,
    )
    shared_mappable = None
    for axis, title, values in zip(
        axes[:2], ("NT triangle", "SOLPS/EIRENE"), (nt, reference)
    ):
        shared_mappable = axis.tripcolor(
            triangulation,
            facecolors=np.ma.masked_invalid(values),
            cmap="turbo",
            norm=shared_norm,
            edgecolors="none",
            rasterized=True,
        )
        axis.set_title(title)

    ratio_mappable = axes[2].tripcolor(
        triangulation,
        facecolors=np.ma.masked_invalid(ratio),
        cmap="coolwarm",
        norm=mcolors.LogNorm(vmin=0.1, vmax=10.0),
        edgecolors="none",
        rasterized=True,
    )
    axes[2].set_title("NT / SOLPS")

    for axis in axes:
        if wall is not None:
            axis.plot(wall[:, 0], wall[:, 1], color="black", linewidth=0.45)
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlim(xlim)
        axis.set_ylim(ylim)
        axis.set_xlabel("R [m]")
        axis.set_ylabel("Z [m]")

    shared_colorbar = fig.colorbar(
        shared_mappable, ax=axes[:2], location="bottom", shrink=0.82, pad=0.08
    )
    shared_colorbar.set_label(unit + " (shared scale)")
    ratio_colorbar = fig.colorbar(
        ratio_mappable, ax=axes[2], location="bottom", shrink=0.82, pad=0.08
    )
    ratio_colorbar.set_label("dimensionless")
    fig.suptitle(name, fontsize=15)
    pdf.savefig(fig, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare NT and SOLPS/EIRENE on the same triangle mesh. "
            "The two physical panels always share one color scale."
        )
    )
    parser.add_argument("--root", default=".")
    parser.add_argument("--case", default="case_input/2MW-5e19")
    parser.add_argument("--output", default="Outputfile/d2p_mesh3_prototype")
    parser.add_argument(
        "--pdf",
        default=None,
        help="Default: CASE/fig/tri_nt_solps_comparison.pdf",
    )
    parser.add_argument(
        "--metrics",
        default=None,
        help="Default: CASE/fig/tri_nt_solps_comparison_metrics.csv",
    )
    parser.add_argument("--xlim", nargs=2, type=float, default=(1.30, 1.80))
    parser.add_argument("--ylim", nargs=2, type=float, default=(-1.12, -0.60))
    args = parser.parse_args()

    root = Path(args.root).resolve()
    case = root / args.case
    output_data = root / args.output / "data"
    figure_dir = case / "fig"
    figure_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = Path(args.pdf) if args.pdf else figure_dir / "tri_nt_solps_comparison.pdf"
    if not pdf_path.is_absolute():
        pdf_path = root / pdf_path
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    r, z = read_fort33(case / "solps_output/fort.33")
    triangles = read_fort34(case / "solps_output/fort.34")
    b2 = read_fort35_b2(case / "solps_output/fort.35")
    tri_volumes = triangle_volume(r, z, triangles)
    triangulation = mtri.Triangulation(r, z, triangles)
    wall_path = case / "shapedata/wall.txt"
    wall = np.loadtxt(wall_path, skiprows=1) if wall_path.exists() else None

    fields = []

    def add_direct(name, unit, nt_name, reference_name):
        nt_path = output_data / nt_name
        reference_path = case / reference_name
        if nt_path.exists() and reference_path.exists():
            fields.append(
                (name, unit, np.loadtxt(nt_path), np.loadtxt(reference_path))
            )

    add_direct("Atomic deuterium density", r"$m^{-3}$", "n_D_0_Tri", "D_0_n_Tri")
    add_direct("Molecular deuterium density", r"$m^{-3}$", "n_D2_0_Tri", "D2_0_n_Tri")

    d2_ion_nt = output_data / "n_D2_1_Tri"
    d2_ion_reference = case / "nDmoleculeion_2D.data"
    if d2_ion_nt.exists() and d2_ion_reference.exists():
        fields.append(
            (
                "Molecular deuterium ion density",
                r"$m^{-3}$",
                np.loadtxt(d2_ion_nt),
                map_b2_value_to_tri(read_b2_matrix(d2_ion_reference), b2),
            )
        )

    ei_terms = {
        "Sn_Ion_D_0_Tri": 1.0,
        "Sn_Diss2_D2_0_Tri": 1.0,
        "Sn_DS1_D2_1_Tri": 2.0,
        "Sn_DS2_D2_1_Tri": 1.0,
    }
    ei_arrays = [
        coefficient * np.loadtxt(output_data / filename)
        for filename, coefficient in ei_terms.items()
        if (output_data / filename).exists()
    ]
    ei_reference = case / "sum_strata_ei_sna_Dion.data"
    if ei_arrays and ei_reference.exists():
        fields.append(
            (
                "Electron-impact D+ production source",
                r"$s^{-1}$ per triangle",
                np.sum(ei_arrays, axis=0),
                distribute_b2_integral_to_tri(
                    read_b2_matrix(ei_reference), b2, tri_volumes
                ),
            )
        )

    dissociation_terms = [
        output_data / "Sn_Diss1_D2_0_Tri",
        output_data / "Sn_Diss2_D2_0_Tri",
    ]
    dissociation_reference = case / "mol_ds_2D.data"
    if all(path.exists() for path in dissociation_terms) and dissociation_reference.exists():
        fields.append(
            (
                "D2 dissociation collision source",
                r"$s^{-1}$ per triangle",
                sum(np.loadtxt(path) for path in dissociation_terms),
                distribute_b2_integral_to_tri(
                    read_b2_matrix(dissociation_reference), b2, tri_volumes
                ),
            )
        )

    if not fields:
        raise RuntimeError("No comparable fields were found.")
    for name, _, nt, reference in fields:
        if len(nt) != len(triangles) or len(reference) != len(triangles):
            raise ValueError(
                f"{name}: expected {len(triangles)} triangle values, "
                f"got NT={len(nt)}, SOLPS={len(reference)}"
            )

    metric_rows = []
    with PdfPages(pdf_path) as pdf:
        metadata = pdf.infodict()
        metadata["Title"] = "NT and SOLPS/EIRENE triangle-grid comparison"
        metadata["Subject"] = "Densities and collision source terms with shared colorbars"
        for name, unit, nt, reference in fields:
            add_comparison_page(
                pdf,
                triangulation,
                wall,
                name,
                unit,
                nt,
                reference,
                args.xlim,
                args.ylim,
            )
            metric_rows.append(metrics(name, nt, reference))

    metric_path = (
        (root / args.metrics).resolve()
        if args.metrics
        else figure_dir / "tri_nt_solps_comparison_metrics.csv"
    )
    metric_path.parent.mkdir(parents=True, exist_ok=True)
    keys = []
    for row in metric_rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with metric_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(metric_rows)

    print(f"Wrote {pdf_path}")
    print(f"Wrote {metric_path}")


if __name__ == "__main__":
    main()
