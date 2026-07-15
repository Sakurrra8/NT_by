#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


SUMMARY_FIELDS = {
    'tri_D': 'n_D_0',
    'tri_D2': 'n_D2_0',
    'b2_D': 'B2_n_D_0',
    'b2_D_inner': 'B2_n_D_0_inner_target',
    'b2_D_inner_r30_36': 'B2_n_D_0_inner_target_r30_36',
    'b2_D_outer': 'B2_n_D_0_outer_target',
    'b2_D2': 'B2_n_D2_0',
    'b2_D2_inner': 'B2_n_D2_0_inner_target',
    'b2_D2_inner_r30_36': 'B2_n_D2_0_inner_target_r30_36',
    'b2_D2_outer': 'B2_n_D2_0_outer_target',
    'b2_T_D_inner': 'B2_T_D_0_inner_target',
    'b2_T_D2_inner': 'B2_T_D2_0_inner_target',
    'b2_D2_sink': 'B2_D2_sink_vs_srcml_all',
    'b2_D2_sink_inner': 'B2_D2_sink_vs_srcml_inner_target',
    'b2_D2_sink_inner_r30_36': 'B2_D2_sink_vs_srcml_inner_target_r30_36',
    'target_D2_shape': 'Target_D2_source_shape_inner',
}


def parse_run(value):
    try:
        density, output = value.split('=', 1)
        return float(density), Path(output)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            'run must use DENSITY_E19=OUTPUT_DIRECTORY'
        ) from error


def read_rows(path, key):
    with path.open(newline='') as stream:
        return {row[key]: row for row in csv.DictReader(stream)}


def first_number(row, keys):
    for key in keys:
        value = row.get(key, '')
        if value != '':
            return float(value)
    return None


def write_csv(path, rows):
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description='Combine NT/EIRENE density-scan benchmark outputs.'
    )
    parser.add_argument('--root', default='.')
    parser.add_argument('--run', action='append', type=parse_run, required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    root = Path(args.root)
    summary_rows = []
    long_rows = []
    for density, relative_output in sorted(args.run):
        output = root / relative_output
        metrics_path = output / 'fig/tri_benchmark/benchmark_tri_metrics.csv'
        metrics = read_rows(metrics_path, 'field')
        job_id = output.name

        summary = {'density_e19_m-3': density, 'job_id': job_id}
        for alias, field in SUMMARY_FIELDS.items():
            row = metrics[field]
            bias = first_number(
                row, ('volume_weighted_rel_bias', 'relative_bias')
            )
            l1 = first_number(row, ('volume_weighted_rel_l1', 'relative_l1'))
            code = first_number(
                row, ('code_volume_integral', 'code_cell_integral')
            )
            reference = first_number(
                row, ('ref_volume_integral', 'ref_cell_integral')
            )
            if bias is not None:
                summary[f'{alias}_bias_pct'] = 100.0 * bias
            if l1 is not None:
                summary[f'{alias}_l1_pct'] = 100.0 * l1
            if code is not None and reference not in (None, 0.0):
                summary[f'{alias}_ratio'] = code / reference
            lifetime = first_number(row, ('inventory_over_sink_ratio',))
            if lifetime is not None:
                summary[f'{alias}_lifetime_ratio'] = lifetime

        conservation = read_rows(
            output / 'fig/deuterium_conservation_check.csv', 'quantity'
        )
        summary['D_nuclei_reaction_residual'] = float(
            conservation['D_nuclei_internal_reaction_residual'][
                'relative_residual'
            ]
        )
        with (output / 'fig/collision_closure.csv').open(newline='') as stream:
            for row in csv.DictReader(stream):
                summary[
                    f"{row['species']}_{row['reaction']}_actual_over_track"
                ] = float(row['actual_over_track'])
        summary_rows.append(summary)

        for field, row in metrics.items():
            combined = {
                'density_e19_m-3': density,
                'job_id': job_id,
                **row,
            }
            long_rows.append(combined)

    outdir = root / args.outdir
    write_csv(outdir / 'density_scan_summary.csv', summary_rows)
    write_csv(outdir / 'density_scan_metrics_long.csv', long_rows)
    print('wrote', outdir / 'density_scan_summary.csv')
    print('wrote', outdir / 'density_scan_metrics_long.csv')


if __name__ == '__main__':
    main()
