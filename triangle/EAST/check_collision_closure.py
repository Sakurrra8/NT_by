#!/usr/bin/env python3
"""Compare sampled electron collisions with the track-length estimator."""

import argparse
import csv
from pathlib import Path

import numpy as np

from benchmark_tri import read_fort33, read_fort34, read_fort35_b2, tri_volume


CHANNELS = (
    (
        "D",
        "Ion",
        "n_D_0_Tri",
        "Sn_Ion_D_0_Tri",
        "R2_1_5_H4_H.4_2.1.5.csv",
    ),
    (
        "D2",
        "Ion",
        "n_D2_0_Tri",
        "Sn_Ion_D2_0_Tri",
        "R2_2_9_H4_H.4_2.2.9.csv",
    ),
    (
        "D2",
        "Diss1",
        "n_D2_0_Tri",
        "Sn_Diss1_D2_0_Tri",
        "R2_2_5g_H4_H.4_2.2.5g.csv",
    ),
    (
        "D2",
        "Diss2",
        "n_D2_0_Tri",
        "Sn_Diss2_D2_0_Tri",
        "R2_2_10_H4_H.4_2.2.10.csv",
    ),
)


def read_frequency(path, shape=(98, 38)):
    frequency = np.zeros(shape)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            i = int(row["i"])
            j = int(row["j"])
            frequency[i, j] = (
                float(row["ne_m-3"]) * float(row["value_m3_s"])
            )
    return frequency


def main():
    parser = argparse.ArgumentParser(
        description="Check sampled D/D2 electron collisions against track length."
    )
    parser.add_argument("--root", default=".")
    parser.add_argument(
        "--case",
        required=True,
        help="Case directory used by the corresponding NT run.",
    )
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--rates",
        required=True,
        help="Case-specific eirene_reaction_b2 directory.",
    )
    parser.add_argument(
        "--csv",
        default=None,
        help="Default: OUTPUT/fig/collision_closure.csv",
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    case = root / args.case
    output = root / args.output
    data_dir = output / "data"
    rate_dir = root / args.rates
    csv_path = Path(args.csv) if args.csv else output / "fig/collision_closure.csv"
    if not csv_path.is_absolute():
        csv_path = root / csv_path

    r, z = read_fort33(case / "solps_output/fort.33")
    triangles = read_fort34(case / "solps_output/fort.34")
    b2 = read_fort35_b2(case / "solps_output/fort.35")
    volumes = tri_volume(r, z, triangles)
    mapped = (
        (b2[:, 0] >= 0)
        & (b2[:, 0] < 98)
        & (b2[:, 1] >= 0)
        & (b2[:, 1] < 38)
    )

    rows = []
    for species, reaction, density_name, source_name, rate_name in CHANNELS:
        density = np.loadtxt(data_dir / density_name)
        source = np.loadtxt(data_dir / source_name)
        frequency = read_frequency(rate_dir / rate_name)
        expected = np.sum(
            density[mapped]
            * volumes[mapped]
            * frequency[b2[mapped, 0], b2[mapped, 1]]
        )
        actual = np.sum(source)
        ratio = actual / expected if expected > 0.0 else np.nan
        rows.append(
            {
                "species": species,
                "reaction": reaction,
                "actual_collision_s-1": actual,
                "track_length_estimate_s-1": expected,
                "actual_over_track": ratio,
                "relative_residual": ratio - 1.0,
            }
        )

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {csv_path}")
    for row in rows:
        print(
            f"{row['species']:>2} {row['reaction']:<5} "
            f"actual/track={row['actual_over_track']:.6f}"
        )


if __name__ == "__main__":
    main()
