#!/usr/bin/env python3
"""Export NT neutral-neutral background moments from an EIRENE fort.46."""

import argparse
from pathlib import Path

import numpy as np


ELECTRON_CHARGE = 1.602176634e-19
ATOMIC_MASS_UNIT = 1.66053906660e-27


def read_field(path, name):
    lines = path.read_text().splitlines()
    marker = f'*eirene data field {name} '
    start = next(
        (index for index, line in enumerate(lines) if line.startswith(marker)),
        None,
    )
    if start is None:
        raise ValueError(f'{name} is missing from {path}')
    tokens = lines[start].split()
    count = int(tokens[tokens.index('size') + 1])
    values = []
    for line in lines[start + 1:]:
        values.extend(float(value) for value in line.split())
        if len(values) >= count:
            break
    if len(values) < count:
        raise ValueError(f'{name} is truncated in {path}')
    return np.asarray(values[:count])


def neutral_moments(path, density_name, energy_name, suffix, mass_number):
    density_cm3 = read_field(path, density_name)
    energy_density = read_field(path, energy_name)
    momentum = np.column_stack([
        read_field(path, f'{component}{suffix}')
        for component in ('vx', 'vy', 'vz')
    ])
    mass_kg = mass_number * ATOMIC_MASS_UNIT
    velocity = np.divide(
        momentum,
        density_cm3[:, None] * (mass_kg * 1.0e3),
        out=np.zeros_like(momentum),
        where=density_cm3[:, None] > 0.0,
    ) * 1.0e-2
    total_energy = np.divide(
        energy_density,
        density_cm3,
        out=np.zeros_like(energy_density),
        where=density_cm3 > 0.0,
    )
    flow_energy = (
        0.5 * mass_kg * np.sum(velocity * velocity, axis=1)
        / ELECTRON_CHARGE
    )
    thermal_temperature = (
        (2.0 / 3.0) * np.maximum(total_energy - flow_energy, 0.0)
    )
    return density_cm3 * 1.0e6, thermal_temperature, velocity


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Write the fixed D/D2 density, thermal-temperature, and drift '
            'backgrounds used by NT neutral-neutral collisions.'
        )
    )
    parser.add_argument('--fort46', required=True, type=Path)
    parser.add_argument('--outdir', required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    for species, density_name, energy_name, suffix, mass_number in (
        ('D', 'pdena', 'edena', 'dena', 2.0),
        ('D2', 'pdenm', 'edenm', 'denm', 4.0),
    ):
        density, thermal_temperature, velocity = neutral_moments(
            args.fort46, density_name, energy_name, suffix, mass_number
        )
        np.savetxt(args.outdir / f'{species}_0_n_Tri', density)
        np.savetxt(args.outdir / f'{species}_0_T_Tri', thermal_temperature)
        np.savetxt(args.outdir / f'{species}_0_V_Tri', velocity)
        print(
            f'{species}: cells={density.size} '
            f'density_sum={np.sum(density):.8e} '
            f'density_weighted_thermal_T_eV='
            f'{np.sum(density * thermal_temperature) / np.sum(density):.8g}'
        )


if __name__ == '__main__':
    main()
