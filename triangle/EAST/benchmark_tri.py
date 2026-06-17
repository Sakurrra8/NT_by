#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

matplotlib.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'stix',
    'font.size': 12,
    'figure.autolayout': True,
})


def read_fort33(path):
    values = Path(path).read_text().split()
    n = int(values[0])
    r = np.array([float(x) for x in values[1:1 + n]]) / 100.0
    z = np.array([float(x) for x in values[1 + n:1 + 2 * n]]) / 100.0
    return r, z


def read_fort34(path):
    raw = np.loadtxt(path, dtype=int, skiprows=1)
    return raw[:, 1:4] - 1


def read_fort35_b2(path):
    rows = np.loadtxt(path, dtype=int, skiprows=1)
    b2 = np.full((rows.shape[0], 2), -1, dtype=int)
    for i, row in enumerate(rows):
        m = row[-2]
        n = row[-1]
        if m < 0 or n < 0:
            continue
        if m > 74:
            b2[i, 0] = m - 2
        elif m > 25:
            b2[i, 0] = m - 1
        else:
            b2[i, 0] = m
        b2[i, 1] = n
    return b2


def tri_volume(r, z, triangles):
    rr = r[triangles]
    zz = z[triangles]
    area = 0.5 * np.abs(
        (rr[:, 1] - rr[:, 0]) * (zz[:, 2] - zz[:, 0])
        - (rr[:, 2] - rr[:, 0]) * (zz[:, 1] - zz[:, 0])
    )
    rc = rr.mean(axis=1)
    return 2.0 * np.pi * rc * area


def read_2d_data(path):
    data = np.loadtxt(path, skiprows=1)
    return np.transpose(data[:, 3].reshape(38, 98))


def map_b2_to_tri(field, b2):
    out = np.full(b2.shape[0], np.nan)
    ok = (
        (b2[:, 0] >= 0)
        & (b2[:, 0] < field.shape[0])
        & (b2[:, 1] >= 0)
        & (b2[:, 1] < field.shape[1])
    )
    out[ok] = field[b2[ok, 0], b2[ok, 1]]
    return out


def metric_row(name, code, ref, vol):
    finite = np.isfinite(code) & np.isfinite(ref)
    positive = finite & (code > 0.0) & (ref > 0.0)
    active = finite & ((code != 0.0) | (ref != 0.0))
    if not np.any(active):
        return {'field': name, 'cells': 0}

    w = vol[active]
    c = code[active]
    r = ref[active]
    denom = np.sum(np.abs(r) * w)
    rel_l1 = np.nan if denom == 0.0 else np.sum(np.abs(c - r) * w) / denom
    rel_bias = np.nan if denom == 0.0 else np.sum((c - r) * w) / denom
    row = {
        'field': name,
        'cells': int(np.sum(active)),
        'positive_cells': int(np.sum(positive)),
        'code_volume_integral': float(np.sum(c * w)),
        'ref_volume_integral': float(np.sum(r * w)),
        'volume_weighted_rel_l1': float(rel_l1),
        'volume_weighted_rel_bias': float(rel_bias),
        'code_min': float(np.nanmin(c)),
        'code_max': float(np.nanmax(c)),
        'ref_min': float(np.nanmin(r)),
        'ref_max': float(np.nanmax(r)),
    }
    if np.any(positive):
        ratio = code[positive] / ref[positive]
        logdiff = np.log10(code[positive]) - np.log10(ref[positive])
        row.update({
            'median_ratio_positive': float(np.nanmedian(ratio)),
            'p10_ratio_positive': float(np.nanpercentile(ratio, 10)),
            'p90_ratio_positive': float(np.nanpercentile(ratio, 90)),
            'log10_rmse_positive': float(np.sqrt(np.nanmean(logdiff * logdiff))),
        })
    return row


def safe_lognorm(values):
    pos = values[np.isfinite(values) & (values > 0.0)]
    if pos.size == 0:
        return None
    vmin = max(np.nanpercentile(pos, 1), np.nanmin(pos[pos > 0.0]))
    vmax = np.nanpercentile(pos, 99)
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
        return None
    return colors.LogNorm(vmin=vmin, vmax=vmax)


def plot_field(outdir, triangulation, wall, name, code, ref, xlim, ylim):
    ratio = np.full_like(code, np.nan, dtype=float)
    ok = np.isfinite(code) & np.isfinite(ref) & (code > 0.0) & (ref > 0.0)
    ratio[ok] = code[ok] / ref[ok]

    fig, ax = plt.subplots(1, 3, figsize=(15, 6), dpi=180)
    panels = [
        ('code', code, safe_lognorm(code), 'jet'),
        ('reference', ref, safe_lognorm(ref), 'jet'),
        ('code/reference', ratio, colors.LogNorm(vmin=1e-2, vmax=1e2), 'coolwarm'),
    ]
    for axis, (title, data, norm, cmap) in zip(ax, panels):
        tc = axis.tripcolor(
            triangulation,
            facecolors=data,
            cmap=cmap,
            norm=norm,
            edgecolors='none',
        )
        axis.set_aspect('equal', adjustable='box')
        axis.set_xlim(xlim)
        axis.set_ylim(ylim)
        axis.set_title(name + ' ' + title)
        axis.set_xlabel('R m')
        axis.set_ylabel('Z m')
        if wall is not None:
            axis.plot(wall[:, 0], wall[:, 1], color='k', lw=0.3)
        fig.colorbar(tc, ax=axis, fraction=0.046, pad=0.04)
    fig.savefig(outdir / (name + '_tri_compare.pdf'))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description='Compare current triangle mesh output with SOLPS/EIRENE triangle references.'
    )
    parser.add_argument('--root', default='.', help='Repository root. Default: current directory')
    parser.add_argument('--case', default='case_input/2MW-5e19')
    parser.add_argument('--output', default='Outputfile/test')
    parser.add_argument('--outdir', default=None, help='Default: output/fig/tri_benchmark')
    parser.add_argument('--xlim', nargs=2, type=float, default=[1.30, 1.80])
    parser.add_argument('--ylim', nargs=2, type=float, default=[-1.12, -0.60])
    args = parser.parse_args()

    root = Path(args.root)
    case = root / args.case
    output_root = root / args.output
    out = output_root / 'data'
    if args.outdir is None:
        outdir = output_root / 'fig' / 'tri_benchmark'
    else:
        outdir = root / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    r, z = read_fort33(case / 'solps_output/fort.33')
    triangles = read_fort34(case / 'solps_output/fort.34')
    b2 = read_fort35_b2(case / 'solps_output/fort.35')
    vol = tri_volume(r, z, triangles)
    triangulation = mtri.Triangulation(r, z, triangles)

    wall_path = case / 'shapedata/wall.txt'
    wall = np.loadtxt(wall_path, skiprows=1) if wall_path.exists() else None

    fields = []
    direct = [
        ('n_D_0', out / 'n_D_0_Tri', case / 'D_0_n_Tri'),
        ('n_D2_0', out / 'n_D2_0_Tri', case / 'D2_0_n_Tri'),
        ('T_D_0', out / 'T_D_0_Tri', case / 'D_0_T_Tri'),
        ('T_D2_0', out / 'T_D2_0_Tri', case / 'D2_0_T_Tri'),
    ]
    for name, code_path, ref_path in direct:
        if code_path.exists() and ref_path.exists():
            fields.append((name, np.loadtxt(code_path), np.loadtxt(ref_path)))

    d2ion_path = out / 'n_D2_1_Tri'
    solps_d2ion_path = case / '2D_data/nDmoleculeion_2D.data'
    if d2ion_path.exists() and solps_d2ion_path.exists():
        fields.append((
            'n_D2_1_mapped',
            np.loadtxt(d2ion_path),
            map_b2_to_tri(read_2d_data(solps_d2ion_path), b2),
        ))

    rows = []
    for name, code, ref in fields:
        if len(code) != len(ref) or len(code) != len(vol):
            raise ValueError(name + ' length mismatch')
        rows.append(metric_row(name, code, ref, vol))
        plot_field(outdir, triangulation, wall, name, code, ref, args.xlim, args.ylim)

    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with open(outdir / 'benchmark_tri_metrics.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(row)
    print('wrote', outdir / 'benchmark_tri_metrics.csv')
    print('figures in', outdir)


if __name__ == '__main__':
    main()
