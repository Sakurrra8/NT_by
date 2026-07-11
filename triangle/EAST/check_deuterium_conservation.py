#!/usr/bin/env python3
"""Check D/D2 reaction-source conservation on NT output files.

The script verifies the internal reaction ledger, i.e. source terms that are
written in Outputfile/<job>/data/Sn_*_Tri. Boundary recycling, wall, target, and
pump balance needs species-resolved fate diagnostics and is reported separately
when such diagnostics are available.
"""

import argparse
import csv
import math
from pathlib import Path


def read_array(data_dir, name, template=None):
    path = data_dir / name
    if not path.exists():
        if template is None:
            return None
        return [0.0] * len(template)
    return [float(value) for value in path.read_text().split()]


def require_array(data_dir, name):
    path = data_dir / name
    if not path.exists():
        raise FileNotFoundError(f"missing required file: {path}")
    return [float(value) for value in path.read_text().split()]


def finite_sum(array):
    return sum(value for value in array if math.isfinite(value))


def relerr(value, scale):
    if scale <= 0.0:
        return 0.0 if abs(value) == 0.0 else np.inf
    return value / scale


def add_row(rows, quantity, terms):
    positive = sum(value for _, value in terms if value > 0.0)
    negative = sum(value for _, value in terms if value < 0.0)
    net = positive + negative
    scale = max(positive, -negative, 0.0)
    rows.append(
        {
            "quantity": quantity,
            "positive_s-1": positive,
            "negative_s-1": negative,
            "net_s-1": net,
            "scale_s-1": scale,
            "relative_residual": relerr(net, scale),
            "terms": "; ".join(f"{name}={value:.8e}" for name, value in terms),
        }
    )


def read_fate_summary(path):
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {
            row["fate"]: float(row["weighted_value_s-1"])
            for row in csv.DictReader(handle)
        }


def read_source_launch_summary(path):
    launched = {"D0": 0.0, "D2": 0.0, "D_nuclei": 0.0}
    if not path.exists():
        return launched
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            species = row["species"]
            charge = int(row["charge"])
            weight = float(row["launched_weight_s-1"])
            if species == "D" and charge == 0:
                launched["D0"] += weight
                launched["D_nuclei"] += weight
            elif species == "D2" and charge == 0:
                launched["D2"] += weight
                launched["D_nuclei"] += 2.0 * weight
    return launched


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "quantity",
        "positive_s-1",
        "negative_s-1",
        "net_s-1",
        "scale_s-1",
        "relative_residual",
        "terms",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Check deuterium reaction-source conservation for NT output."
    )
    parser.add_argument("--root", default=".")
    parser.add_argument("--output", default="Outputfile/10283830")
    parser.add_argument(
        "--csv",
        default=None,
        help="Default: OUTPUT/fig/deuterium_conservation_check.csv",
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    output = root / args.output
    data_dir = output / "data"
    if not data_dir.exists():
        raise FileNotFoundError(f"missing output data directory: {data_dir}")

    template = require_array(data_dir, "Sn_Ion_D_0_Tri")
    names = [
        "Sn_Ion_D_0_Tri",
        "Sn_Rec_D_1_Tri",
        "Sn_Diss1_D2_0_Tri",
        "Sn_Diss2_D2_0_Tri",
        "Sn_Ion_D2_0_Tri",
        "Sn_CX_D2_0_Tri",
        "Sn_DS1_D2_1_Tri",
        "Sn_DS2_D2_1_Tri",
        "Sn_DS3_D2_1_Tri",
        "Sn_DS4_D2_1_Tri",
    ]
    arrays = {name: read_array(data_dir, name, template) for name in names}
    sums = {name: finite_sum(array) for name, array in arrays.items()}

    rows = []
    add_row(
        rows,
        "D0_internal_reaction_source",
        [
            ("-Ion_D0", -sums["Sn_Ion_D_0_Tri"]),
            ("+Rec_D1", sums["Sn_Rec_D_1_Tri"]),
            ("+2*Diss1_D2", 2.0 * sums["Sn_Diss1_D2_0_Tri"]),
            ("+Diss2_D2", sums["Sn_Diss2_D2_0_Tri"]),
            ("+CX_D2", sums["Sn_CX_D2_0_Tri"]),
            ("+DS2_D2p", sums["Sn_DS2_D2_1_Tri"]),
            ("+2*DS3_D2p", 2.0 * sums["Sn_DS3_D2_1_Tri"]),
        ],
    )
    add_row(
        rows,
        "Dplus_internal_reaction_source",
        [
            ("+Ion_D0", sums["Sn_Ion_D_0_Tri"]),
            ("-Rec_D1", -sums["Sn_Rec_D_1_Tri"]),
            ("+Diss2_D2", sums["Sn_Diss2_D2_0_Tri"]),
            ("-CX_D2", -sums["Sn_CX_D2_0_Tri"]),
            ("+2*DS1_D2p", 2.0 * sums["Sn_DS1_D2_1_Tri"]),
            ("+DS2_D2p", sums["Sn_DS2_D2_1_Tri"]),
        ],
    )
    add_row(
        rows,
        "D2_internal_reaction_source",
        [
            ("-Diss1_D2", -sums["Sn_Diss1_D2_0_Tri"]),
            ("-Diss2_D2", -sums["Sn_Diss2_D2_0_Tri"]),
            ("-Ion_D2", -sums["Sn_Ion_D2_0_Tri"]),
            ("-CX_D2", -sums["Sn_CX_D2_0_Tri"]),
        ],
    )
    add_row(
        rows,
        "D2plus_internal_reaction_source",
        [
            ("+Ion_D2", sums["Sn_Ion_D2_0_Tri"]),
            ("+CX_D2", sums["Sn_CX_D2_0_Tri"]),
            ("-DS1_D2p", -sums["Sn_DS1_D2_1_Tri"]),
            ("-DS2_D2p", -sums["Sn_DS2_D2_1_Tri"]),
            ("-DS3_D2p", -sums["Sn_DS3_D2_1_Tri"]),
            ("-DS4_D2p", -sums["Sn_DS4_D2_1_Tri"]),
        ],
    )

    d_nuclei_terms = [
        ("D0", rows[0]["net_s-1"]),
        ("Dplus", rows[1]["net_s-1"]),
        ("2*D2", 2.0 * rows[2]["net_s-1"]),
        ("2*D2plus", 2.0 * rows[3]["net_s-1"]),
    ]
    add_row(rows, "D_nuclei_internal_reaction_residual", d_nuclei_terms)

    fate = read_fate_summary(data_dir / "D2p_mesh3_tracking_fate_summary.csv")
    if fate:
        created = fate.get("D2p_created_by_ion", 0.0) + fate.get(
            "D2p_created_by_cx", 0.0
        )
        lost = (
            fate.get("D2p_DS1_events", 0.0)
            + fate.get("D2p_DS2_events", 0.0)
            + fate.get("D2p_DS3_events", 0.0)
            + fate.get("D2p_boundary_loss", 0.0)
            + fate.get("D2p_max_steps_loss", 0.0)
        )
        add_row(
            rows,
            "D2plus_transport_fate_residual",
            [("+created", created), ("-tracked_fates", -lost)],
        )

    launched = read_source_launch_summary(data_dir / "source_stratum_summary.csv")
    add_row(
        rows,
        "reported_neutral_launch_D0_not_closure",
        [("+launched_D0", launched["D0"])],
    )
    add_row(
        rows,
        "reported_neutral_launch_D2_not_closure",
        [("+launched_D2", launched["D2"])],
    )
    add_row(
        rows,
        "reported_neutral_launch_D_nuclei_not_closure",
        [("+launched_D_nuclei", launched["D_nuclei"])],
    )

    csv_path = Path(args.csv) if args.csv else output / "fig" / "deuterium_conservation_check.csv"
    if not csv_path.is_absolute():
        csv_path = root / csv_path
    write_csv(csv_path, rows)

    print(f"Wrote {csv_path}")
    print("quantity,net_s-1,relative_residual")
    for row in rows:
        print(
            f"{row['quantity']},{row['net_s-1']:.8e},"
            f"{row['relative_residual']:.8e}"
        )
    print(
        "Note: neutral wall/target/pump closure is not included unless "
        "species-resolved neutral fate diagnostics are written. Rows ending "
        "in not_closure are launch magnitudes only."
    )


if __name__ == "__main__":
    main()
