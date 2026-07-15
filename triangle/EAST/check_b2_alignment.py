#!/usr/bin/env python3
"""Rank possible B2 axis flips and integer offsets against an EIRENE field."""

import argparse
import io
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_tri import read_b2_matrix, read_eirene_field


def read_code_matrix(path):
    if path != '-':
        return read_b2_matrix(path)
    data = np.loadtxt(io.StringIO(sys.stdin.read().lstrip('\ufeff')))
    if data.shape == (98, 38):
        return data
    if data.shape == (38, 98):
        return data.T
    if data.ndim == 2 and data.shape[1] == 4:
        return np.transpose(data[:, 3].reshape(38, 98))
    raise ValueError('stdin is not a recognized B2 98x38 field')


def overlapping(a, b, shift_i, shift_j):
    ni, nj = a.shape
    ai = slice(max(0, shift_i), min(ni, ni + shift_i))
    aj = slice(max(0, shift_j), min(nj, nj + shift_j))
    bi = slice(max(0, -shift_i), min(ni, ni - shift_i))
    bj = slice(max(0, -shift_j), min(nj, nj - shift_j))
    return a[ai, aj], b[bi, bj]


def score(code, reference, flip_i, flip_j, shift_i, shift_j):
    candidate = reference
    if flip_i:
        candidate = np.flip(candidate, axis=0)
    if flip_j:
        candidate = np.flip(candidate, axis=1)
    code_overlap, reference_overlap = overlapping(
        code, candidate, shift_i, shift_j
    )
    positive = (code_overlap > 0.0) & (reference_overlap > 0.0)
    if np.sum(positive) < 2:
        return None
    log_code = np.log10(code_overlap[positive])
    log_reference = np.log10(reference_overlap[positive])
    return {
        'log10_rmse': float(np.sqrt(np.mean((log_code - log_reference) ** 2))),
        'log10_correlation': float(np.corrcoef(log_code, log_reference)[0, 1]),
        'median_ratio': float(np.median(code_overlap[positive] / reference_overlap[positive])),
        'positive_cells': int(np.sum(positive)),
        'flip_i': flip_i,
        'flip_j': flip_j,
        'shift_i': shift_i,
        'shift_j': shift_j,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--code', required=True, help='Code 98x38 B2 field, or - for stdin'
    )
    parser.add_argument('--fort44', required=True)
    parser.add_argument('--field', default='dmb2')
    parser.add_argument('--max-shift', type=int, default=3)
    parser.add_argument('--top', type=int, default=12)
    args = parser.parse_args()

    code = read_code_matrix(args.code)[1:97, 1:37]
    reference = read_eirene_field(args.fort44, args.field)
    if reference.size != 96 * 36:
        raise ValueError(f'{args.field} has unexpected size {reference.size}')
    reference = reference.reshape(36, 96).T

    rows = []
    for flip_i in (False, True):
        for flip_j in (False, True):
            for shift_i in range(-args.max_shift, args.max_shift + 1):
                for shift_j in range(-args.max_shift, args.max_shift + 1):
                    row = score(
                        code, reference, flip_i, flip_j, shift_i, shift_j
                    )
                    if row is not None:
                        rows.append(row)
    rows.sort(key=lambda row: (row['log10_rmse'], -row['log10_correlation']))
    print(
        'rank,log10_rmse,log10_correlation,median_ratio,positive_cells,'
        'flip_i,flip_j,shift_i,shift_j'
    )
    for rank, row in enumerate(rows[:args.top], start=1):
        print(
            f"{rank},{row['log10_rmse']},{row['log10_correlation']},"
            f"{row['median_ratio']},{row['positive_cells']},"
            f"{int(row['flip_i'])},{int(row['flip_j'])},"
            f"{row['shift_i']},{row['shift_j']}"
        )


if __name__ == '__main__':
    main()
