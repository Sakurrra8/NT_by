#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path

import numpy as np


ELECTRON_CHARGE = 1.602176634e-19
ATOMIC_MASS_UNIT = 1.66053906660e-27
NT_MASSES = {
    'D': 3.34e-27,
    'D2': 6.695e-27,
}
MASS_NUMBERS = {
    'D': 2.0,
    'D2': 4.0,
}

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    from matplotlib.backends.backend_pdf import PdfPages

    matplotlib.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['DejaVu Serif'],
        'mathtext.fontset': 'stix',
        'font.size': 12,
        'figure.autolayout': True,
    })
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    colors = None
    plt = None
    mtri = None
    PdfPages = None


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


def read_eirene_field(path, field_name):
    lines = Path(path).read_text().splitlines()
    marker = f'*eirene data field {field_name} '
    marker_index = next(
        (index for index, line in enumerate(lines) if line.startswith(marker)),
        None,
    )
    if marker_index is None:
        raise ValueError(f'{field_name} is missing from {path}')
    marker_values = lines[marker_index].split()
    count = int(marker_values[marker_values.index('size') + 1])
    values = []
    line_index = marker_index + 1
    while len(values) < count:
        values.extend(float(value) for value in lines[line_index].split())
        line_index += 1
    return np.asarray(values[:count])


def read_eirene_neutral_triangles(path):
    fields = {}
    for species, density_name, energy_name, momentum_suffix in (
        ('D', 'pdena', 'edena', 'dena'),
        ('D2', 'pdenm', 'edenm', 'denm'),
    ):
        density_cm3 = read_eirene_field(path, density_name)
        energy_density = read_eirene_field(path, energy_name)
        fields[f'n_{species}_0'] = density_cm3 * 1.0e6
        total_energy = np.divide(
            energy_density,
            density_cm3,
            out=np.zeros_like(energy_density),
            where=density_cm3 > 0.0,
        )
        fields[f'T_{species}_0'] = (2.0 / 3.0) * total_energy

        momentum = np.column_stack([
            read_eirene_field(path, f'{component}{momentum_suffix}')
            for component in ('vx', 'vy', 'vz')
        ])
        mass_kg = MASS_NUMBERS[species] * ATOMIC_MASS_UNIT
        velocity = np.divide(
            momentum,
            density_cm3[:, None] * (mass_kg * 1.0e3),
            out=np.zeros_like(momentum),
            where=density_cm3[:, None] > 0.0,
        ) * 1.0e-2
        flow_energy = (
            0.5 * mass_kg * np.sum(velocity * velocity, axis=1)
            / ELECTRON_CHARGE
        )
        fields[f'V_{species}_0'] = velocity
        fields[f'T_{species}_0_thermal'] = (
            (2.0 / 3.0) * np.maximum(total_energy - flow_energy, 0.0)
        )
    return fields


def read_eirene_b2_field(path, field_name, scale=1.0):
    active = read_eirene_field(path, field_name)
    if active.size != 96 * 36:
        raise ValueError(f'{field_name} has unexpected B2 size {active.size}')
    field = np.zeros((98, 38))
    field[1:97, 1:37] = active.reshape(36, 96).T * scale
    return field


def read_named_setting(path, name):
    for raw_line in Path(path).read_text().splitlines():
        tokens = raw_line.split('#', 1)[0].split()
        if len(tokens) >= 2 and tokens[0] == name:
            return tokens[1]
    return None


def read_eirene_neutral_reactions(path):
    lines = Path(path).read_text().splitlines()
    reactions = {}
    in_neutral_species = False
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('** 4a.'):
            in_neutral_species = True
            continue
        if stripped.startswith('** 4c.'):
            break
        if not in_neutral_species:
            continue

        tokens = stripped.split()
        if len(tokens) < 7:
            continue
        try:
            integer_values = [int(value) for value in tokens]
        except ValueError:
            continue
        reaction = integer_values[0]
        if reaction not in (1, 2, 6, 7) or index + 1 >= len(lines):
            continue
        try:
            coefficients = [float(value) for value in lines[index + 1].split()]
        except ValueError:
            continue
        if len(coefficients) < 5:
            continue
        reactions[reaction] = {
            'freac': coefficients[4],
            'ibgk': integer_values[6],
        }
    return reactions


def reference_physics_audit_rows(input_path, log_path, setting_path):
    rows = []
    reactions = (
        read_eirene_neutral_reactions(input_path)
        if input_path is not None and input_path.exists()
        else {}
    )
    for reaction in (1, 2, 6, 7):
        if reaction not in reactions:
            continue
        rows.append({
            'field': f'EIRENE_reaction_{reaction}_FREAC',
            'mesh': 'reference_audit',
            'reference_value': reactions[reaction]['freac'],
            'eirene_ibgk': reactions[reaction]['ibgk'],
            'note': 'neutral-neutral AMMONX species-deck multiplier',
        })

    if len(reactions) == 4:
        enabled = {
            reaction: values['freac'] > 1.0e-20
            for reaction, values in reactions.items()
        }
        if all(enabled.values()):
            eirene_state = 'enabled'
            expected_nt = 1
        elif not any(enabled.values()):
            eirene_state = 'disabled'
            expected_nt = 0
        else:
            eirene_state = 'mixed'
            expected_nt = None
        nt_value = None
        if setting_path is not None and setting_path.exists():
            raw_value = read_named_setting(setting_path, 'K_NNCs')
            nt_value = int(raw_value) if raw_value is not None else None
        status = 'not_checked'
        if expected_nt is None:
            status = 'cannot_match_mixed_eirene_deck'
        elif nt_value is not None:
            status = 'match' if nt_value == expected_nt else 'mismatch'
        rows.append({
            'field': 'Neutral_neutral_collision_switch',
            'mesh': 'reference_audit',
            'eirene_state': eirene_state,
            'nt_K_NNCs': nt_value,
            'expected_nt_K_NNCs': expected_nt,
            'status': status,
            'note': (
                'Switch parity only. EIRENE IBGK=111/112 updates virtual '
                'neutral backgrounds, while NT samples binary events against '
                'a fixed background.'
            ),
        })

    if log_path is None or not log_path.exists():
        return rows
    log_text = log_path.read_text(errors='replace')
    bgk_matches = list(re.finditer(
        r'MODBGK CALLED AFTER ITERATION IITER=\s*(\d+)', log_text
    ))
    iterations = [int(match.group(1)) for match in bgk_matches]
    if not iterations:
        rows.append({
            'field': 'EIRENE_BGK_export_order',
            'mesh': 'reference_audit',
            'bgk_calls': 0,
            'status': 'bgk_not_called',
            'reference_log': str(log_path),
        })
        return rows

    fort46_write = log_text.rfind('Writing fort.46')
    final_bgk = bgk_matches[-1].start()
    written_before_bgk = (
        fort46_write >= 0 and fort46_write < final_bgk
    )
    density_residuals = [
        float(value)
        for value in re.findall(
            r'RESDEN:\s+TOT,\s+DEL,\s+REL\s+\S+\s+\S+\s+([0-9.Ee+-]+)',
            log_text,
        )
    ]
    energy_residuals = [
        float(value)
        for value in re.findall(
            r'RESENE:\s+TOT,\s+DEL,\s+REL\s+\S+\s+\S+\s+([0-9.Ee+-]+)',
            log_text,
        )
    ]
    final_bgk_text = log_text[final_bgk:]
    final_density_residuals = [
        float(value)
        for value in re.findall(
            r'RESDEN:\s+TOT,\s+DEL,\s+REL\s+\S+\s+\S+\s+([0-9.Ee+-]+)',
            final_bgk_text,
        )
    ]
    final_energy_residuals = [
        float(value)
        for value in re.findall(
            r'RESENE:\s+TOT,\s+DEL,\s+REL\s+\S+\s+\S+\s+([0-9.Ee+-]+)',
            final_bgk_text,
        )
    ]
    rows.append({
        'field': 'EIRENE_BGK_export_order',
        'mesh': 'reference_audit',
        'bgk_calls': len(iterations),
        'last_bgk_iteration': iterations[-1],
        'fort46_written_before_bgk': int(written_before_bgk),
        'max_density_relative_change': (
            max(density_residuals) if density_residuals else np.nan
        ),
        'max_energy_relative_change': (
            max(energy_residuals) if energy_residuals else np.nan
        ),
        'final_bgk_max_density_relative_change': (
            max(final_density_residuals)
            if final_density_residuals else np.nan
        ),
        'final_bgk_max_energy_relative_change': (
            max(final_energy_residuals)
            if final_energy_residuals else np.nan
        ),
        'status': 'pre_bgk_export' if written_before_bgk else 'post_bgk_export',
        'reference_log': str(log_path),
        'note': (
            'The final fort.46 write is compared with the final MODBGK call. '
            'A pre-BGK file omits that last reported background update.'
        ),
    })
    return rows


def select_reference_log(eirene_output):
    for name in ('b2.log', 'run.log', 'eirene_output'):
        candidate = eirene_output / name
        if candidate.exists():
            return candidate
    return None


def neutral_background_rows(case, reference, volume):
    rows = []
    for species in ('D', 'D2'):
        paths = {
            'density': case / f'{species}_0_n_Tri',
            'temperature': case / f'{species}_0_T_Tri',
            'velocity': case / f'{species}_0_V_Tri',
        }
        if not all(path.exists() for path in paths.values()):
            continue
        density = np.loadtxt(paths['density'])
        thermal_temperature = np.maximum(
            np.loadtxt(paths['temperature']), 0.1
        )
        velocity = np.loadtxt(paths['velocity'])
        if (
            len(density) != len(volume)
            or len(thermal_temperature) != len(volume)
            or velocity.shape != (len(volume), 3)
        ):
            raise ValueError(f'{species} runtime neutral background length mismatch')

        comparisons = (
            (
                f'InputBackground_n_{species}_0',
                density,
                reference[f'n_{species}_0'],
                'runtime fixed density compared with the selected fort.46',
            ),
            (
                f'InputBackground_T_{species}_0_thermal',
                thermal_temperature,
                reference[f'T_{species}_0_thermal'],
                'runtime thermal temperature compared after subtracting EIRENE flow energy',
            ),
            (
                f'InputBackground_speed_{species}_0',
                np.linalg.norm(velocity, axis=1),
                np.linalg.norm(reference[f'V_{species}_0'], axis=1),
                'speed magnitude; independent of an in-plane coordinate rotation',
            ),
        )
        for name, code, ref, note in comparisons:
            row = metric_row(name, code, ref, volume)
            row['mesh'] = 'input_background'
            row['note'] = note
            rows.append(row)

        effective_total_temperature = (
            thermal_temperature
            + NT_MASSES[species]
            * np.sum(velocity * velocity, axis=1)
            / (3.0 * ELECTRON_CHARGE)
        )
        row = metric_row(
            f'InputBackground_effective_total_T_{species}_0',
            effective_total_temperature,
            reference[f'T_{species}_0'],
            volume,
        )
        row['mesh'] = 'input_background'
        row['note'] = (
            'temperature equivalent of the total energy sampled by NT: '
            'thermal Maxwellian plus fixed drift'
        )
        rows.append(row)

        reference_density = reference[f'n_{species}_0']
        reference_velocity = reference[f'V_{species}_0']
        reference_total_temperature = reference[f'T_{species}_0']
        energy_fields = (
            (
                f'InputBackground_flow_energy_density_{species}_0',
                density
                * 0.5
                * NT_MASSES[species]
                * np.sum(velocity * velocity, axis=1)
                / ELECTRON_CHARGE,
                reference_density
                * 0.5
                * MASS_NUMBERS[species]
                * ATOMIC_MASS_UNIT
                * np.sum(reference_velocity * reference_velocity, axis=1)
                / ELECTRON_CHARGE,
                'drift kinetic-energy density in eV m^-3',
            ),
            (
                f'InputBackground_total_energy_density_{species}_0',
                density * 1.5 * effective_total_temperature,
                reference_density * 1.5 * reference_total_temperature,
                'total sampled energy density in eV m^-3',
            ),
        )
        for name, code, ref, note in energy_fields:
            row = metric_row(name, code, ref, volume)
            row['mesh'] = 'input_background'
            row['note'] = note
            rows.append(row)
    return rows


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


def read_b2_matrix(path):
    data = np.loadtxt(path)
    if data.ndim == 2 and data.shape == (98, 38):
        return data
    if data.ndim == 2 and data.shape == (38, 98):
        return data.T
    if data.ndim == 2 and data.shape[1] == 4:
        return np.transpose(data[:, 3].reshape(38, 98))
    raise ValueError(f'{path} is not a recognized B2 98x38 field')


def source_stratum_rows(summary_path, fort44_path):
    eirene_fields = {
        ('D', 'tri'): read_eirene_field(fort44_path, 'pdena_int'),
        ('D2', 'tri'): read_eirene_field(fort44_path, 'pdenm_int'),
        ('D', 'b2'): read_eirene_field(fort44_path, 'pdena_int_b2'),
        ('D2', 'b2'): read_eirene_field(fort44_path, 'pdenm_int_b2'),
    }
    stratum_index = {
        'IT': 1,
        'OT': 2,
        'PFRSide1': 3,
        'PFRSide2': 4,
        'OuterSide': 5,
        'Recombination': 6,
    }
    rows = []
    code_totals = {
        (species, mesh): 0.0
        for species in ('D', 'D2')
        for mesh in ('tri', 'b2')
    }
    with Path(summary_path).open(newline='') as stream:
        reader = csv.DictReader(stream)
        inventory_columns = {
            'b2_track_length_inventory',
            'tri_track_length_inventory',
        }
        if not inventory_columns.issubset(reader.fieldnames or []):
            return rows
        for record in reader:
            species = record['species']
            stratum = record['source_stratum']
            if record['charge'] != '0' or species not in ('D', 'D2'):
                continue
            if stratum not in stratum_index:
                continue
            for mesh, code_column in (
                ('tri', 'tri_track_length_inventory'),
                ('b2', 'b2_track_length_inventory'),
            ):
                code = float(record[code_column])
                code_totals[(species, mesh)] += code
                ref = eirene_fields[(species, mesh)][stratum_index[stratum]]
                row = metric_row(
                    f'SourceStratum_{mesh}_n_{species}_0_{stratum}',
                    np.asarray([code]),
                    np.asarray([ref]),
                    np.ones(1),
                )
                row['mesh'] = mesh
                row['source_stratum'] = stratum
                row['note'] = 'track-length inventory compared with fort.44 source-stratum integral'
                rows.append(row)
    for species in ('D', 'D2'):
        for mesh in ('tri', 'b2'):
            code = code_totals[(species, mesh)]
            ref = eirene_fields[(species, mesh)][0]
            row = metric_row(
                f'SourceStratum_{mesh}_n_{species}_0_total',
                np.asarray([code]),
                np.asarray([ref]),
                np.ones(1),
            )
            row['mesh'] = mesh
            row['source_stratum'] = 'Total'
            row['note'] = (
                'sum of the six matched primary-source track-length inventories '
                'compared with fort.44 pden*_int index 0'
            )
            rows.append(row)
    return rows


def b2_source_spatial_rows(spatial_path, summary_path, density, volume):
    shape = density.shape
    by_source = {}
    with Path(spatial_path).open(newline='') as stream:
        reader = csv.DictReader(stream)
        required = {
            'i', 'j', 'source_stratum', 'track_length_inventory',
        }
        if not required.issubset(reader.fieldnames or []):
            return []
        for record in reader:
            i = int(record['i'])
            j = int(record['j'])
            if not (0 <= i < shape[0] and 0 <= j < shape[1]):
                continue
            source = record['source_stratum']
            field = by_source.setdefault(source, np.zeros(shape))
            field[i, j] += float(record['track_length_inventory'])

    summary = {}
    with Path(summary_path).open(newline='') as stream:
        for record in csv.DictReader(stream):
            if record['species'] == 'D2' and record['charge'] == '0':
                summary[record['source_stratum']] = float(
                    record['b2_track_length_inventory']
                )

    masks = {
        'all': np.ones(shape, dtype=bool),
        'inner_target': np.zeros(shape, dtype=bool),
        'inner_target_r30_36': np.zeros(shape, dtype=bool),
        'outer_target': np.zeros(shape, dtype=bool),
        'outer_target_r30_36': np.zeros(shape, dtype=bool),
    }
    masks['inner_target'][1, 1:37] = True
    masks['inner_target_r30_36'][1, 30:37] = True
    masks['outer_target'][96, 1:37] = True
    masks['outer_target_r30_36'][96, 30:37] = True

    rows = []
    for region, mask in masks.items():
        source_inventory = {
            source: float(np.sum(field[mask]))
            for source, field in by_source.items()
        }
        spatial_total = sum(source_inventory.values())
        density_total = float(np.sum(density[mask] * volume[mask]))
        rows.append({
            'field': f'B2_D2_source_spatial_closure_{region}',
            'mesh': 'b2',
            'region': region,
            'source_spatial_inventory': spatial_total,
            'density_volume_integral': density_total,
            'spatial_to_density_ratio': (
                spatial_total / density_total if density_total > 0.0 else np.nan
            ),
            'note': 'sum of source-resolved track length divided by n_D2_0 times B2 volume',
        })
        for source, inventory in source_inventory.items():
            if inventory == 0.0:
                continue
            row = {
                'field': f'B2_D2_inventory_by_source_{region}_{source}',
                'mesh': 'b2',
                'region': region,
                'source_stratum': source,
                'code_volume_integral': inventory,
                'fraction_of_region': (
                    inventory / spatial_total if spatial_total > 0.0 else np.nan
                ),
                'note': 'D2 track-length inventory attributed to the primary source stratum',
            }
            if region == 'all' and source in summary:
                row['summary_inventory'] = summary[source]
                row['spatial_to_summary_ratio'] = (
                    inventory / summary[source]
                    if summary[source] > 0.0 else np.nan
                )
            rows.append(row)
    return rows


def inner_target_d2_source_shape_row(
    output_data, case, fort44_path, target_profiles
):
    target_path = output_data / 'target_incident_D.csv'
    area_path = case / '2D_data/sarea_l.data'
    if not target_path.exists() or not fort44_path.exists() or not area_path.exists():
        return None

    with target_path.open(newline='') as stream:
        target_rows = {
            int(record['target']): record
            for record in csv.DictReader(stream)
        }
    if not all(index in target_rows for index in range(1, 37)):
        return None

    code = np.asarray([
        float(target_rows[index]['D2_source_s-1'])
        for index in range(1, 37)
    ])
    resolved_area = read_eirene_field(fort44_path, 'sarea_res')
    resolved_molecule_recycling = read_eirene_field(
        fort44_path, 'ewldmr_res'
    )
    inner_segment_area = np.abs(np.loadtxt(area_path)[:, 1])[1:37]
    segment_offset = next(
        (
            offset
            for offset in range(len(resolved_area) - len(inner_segment_area) + 1)
            if np.allclose(
                resolved_area[offset:offset + len(inner_segment_area)],
                inner_segment_area,
                rtol=1.0e-10,
                atol=1.0e-12,
            )
        ),
        None,
    )
    if segment_offset is None:
        return None

    reference = (
        resolved_molecule_recycling[
            segment_offset:segment_offset + len(inner_segment_area)
        ]
        * inner_segment_area
    )
    if np.sum(code) <= 0.0 or np.sum(reference) <= 0.0:
        return None

    code_share = code / np.sum(code)
    reference_share = reference / np.sum(reference)
    for radial_index, (code_value, ref_value) in enumerate(
        zip(code_share, reference_share), start=1
    ):
        target_profiles.append({
            'field': 'Target_D2_source_shape',
            'side': 'inner',
            'radial_index': radial_index,
            'code': code_value,
            'reference': ref_value,
            'ratio': code_value / ref_value if ref_value > 0.0 else np.nan,
        })

    row = source_metric_row(
        'Target_D2_source_shape_inner', code_share, reference_share
    )
    row['mesh'] = 'target'
    row['region'] = 'inner_target_radial_1_36'
    row['note'] = (
        'normalized radial shape only: code D2_source_s-1 versus '
        'EIRENE ewldmr_res times resolved segment area; absolute units differ'
    )
    return row


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


def map_tri_density_to_b2(field, tri_vol, b2, shape):
    if len(field) != len(tri_vol) or len(field) != len(b2):
        raise ValueError('triangle field, volume, and B2 map length mismatch')
    inventory = np.zeros(shape)
    mapped_volume = np.zeros(shape)
    valid = (
        (b2[:, 0] >= 0)
        & (b2[:, 0] < shape[0])
        & (b2[:, 1] >= 0)
        & (b2[:, 1] < shape[1])
    )
    np.add.at(
        inventory,
        (b2[valid, 0], b2[valid, 1]),
        field[valid] * tri_vol[valid],
    )
    np.add.at(
        mapped_volume,
        (b2[valid, 0], b2[valid, 1]),
        tri_vol[valid],
    )
    density = np.divide(
        inventory,
        mapped_volume,
        out=np.zeros_like(inventory),
        where=mapped_volume > 0.0,
    )
    return density, mapped_volume


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


def density_weighted_temperature_rows(
    prefix, species, code_density, ref_density, code_temperature,
    ref_temperature, volume,
):
    finite = (
        np.isfinite(code_density)
        & np.isfinite(ref_density)
        & np.isfinite(code_temperature)
        & np.isfinite(ref_temperature)
        & np.isfinite(volume)
    )
    active = finite & ((code_density > 0.0) | (ref_density > 0.0))
    if not np.any(active):
        return [{'field': f'{prefix}_density_weighted_mean_T_{species}_0', 'cells': 0}]

    code_inventory = np.sum(code_density[active] * volume[active])
    ref_inventory = np.sum(ref_density[active] * volume[active])
    code_n_t = np.sum(
        code_density[active] * code_temperature[active] * volume[active]
    )
    ref_n_t = np.sum(
        ref_density[active] * ref_temperature[active] * volume[active]
    )
    code_mean = code_n_t / code_inventory if code_inventory > 0.0 else np.nan
    ref_mean = ref_n_t / ref_inventory if ref_inventory > 0.0 else np.nan
    denominator = abs(ref_mean)
    bias = (code_mean - ref_mean) / denominator if denominator > 0.0 else np.nan
    mean_row = {
        'field': f'{prefix}_density_weighted_mean_T_{species}_0',
        'cells': int(np.sum(active)),
        'code_volume_integral': float(code_mean),
        'ref_volume_integral': float(ref_mean),
        'volume_weighted_rel_l1': float(abs(bias)),
        'volume_weighted_rel_bias': float(bias),
        'code_density_integral': float(code_inventory),
        'ref_density_integral': float(ref_inventory),
        'code_nT_integral': float(code_n_t),
        'ref_nT_integral': float(ref_n_t),
        'note': 'density-weighted mean temperature: integral(n*T*dV) / integral(n*dV)',
    }
    energy_row = metric_row(
        f'{prefix}_nT_{species}_0',
        code_density * code_temperature,
        ref_density * ref_temperature,
        volume,
    )
    energy_row['note'] = 'neutral thermal-energy proxy n*T in eV m^-3'
    return [mean_row, energy_row]


def source_metric_row(name, code, ref):
    finite = np.isfinite(code) & np.isfinite(ref)
    positive = finite & (code > 0.0) & (ref > 0.0)
    active = finite & ((code != 0.0) | (ref != 0.0))
    if not np.any(active):
        return {'field': name, 'cells': 0}

    c = code[active]
    r = ref[active]
    ref_sum = np.sum(r)
    row = {
        'field': name,
        'cells': int(np.sum(active)),
        'positive_cells': int(np.sum(positive)),
        'code_cell_integral': float(np.sum(c)),
        'ref_cell_integral': float(ref_sum),
        'relative_l1': float(np.sum(np.abs(c - r)) / np.sum(np.abs(r))),
        'relative_bias': float((np.sum(c) - ref_sum) / ref_sum),
        'pearson_active': float(np.corrcoef(c, r)[0, 1]),
        'code_min': float(np.nanmin(c)),
        'code_max': float(np.nanmax(c)),
        'ref_min': float(np.nanmin(r)),
        'ref_max': float(np.nanmax(r)),
    }
    if np.any(positive):
        ratio = code[positive] / ref[positive]
        log_code = np.log10(code[positive])
        log_ref = np.log10(ref[positive])
        row.update({
            'median_ratio_positive': float(np.nanmedian(ratio)),
            'p10_ratio_positive': float(np.nanpercentile(ratio, 10)),
            'p90_ratio_positive': float(np.nanpercentile(ratio, 90)),
            'log10_rmse_positive': float(np.sqrt(np.nanmean((log_code - log_ref) ** 2))),
            'log10_pearson_positive': float(np.corrcoef(log_code, log_ref)[0, 1]),
        })
    return row


def b2_tri_volume_row(b2_vol, tri_vol, b2):
    tri_sum = np.zeros_like(b2_vol)
    tri_count = np.zeros_like(b2_vol, dtype=int)
    valid = (
        (b2[:, 0] >= 0)
        & (b2[:, 0] < b2_vol.shape[0])
        & (b2[:, 1] >= 0)
        & (b2[:, 1] < b2_vol.shape[1])
    )
    np.add.at(tri_sum, (b2[valid, 0], b2[valid, 1]), tri_vol[valid])
    np.add.at(tri_count, (b2[valid, 0], b2[valid, 1]), 1)
    covered = (tri_count > 0) & (b2_vol > 0.0)
    ratio = tri_sum[covered] / b2_vol[covered]
    return {
        'field': 'B2_vs_mapped_triangle_volume',
        'mesh': 'geometry',
        'cells': int(np.sum(covered)),
        'mapped_triangles': int(np.sum(valid)),
        'tri_volume_sum': float(np.sum(tri_vol[valid])),
        'b2_volume_sum': float(np.sum(b2_vol[covered])),
        'volume_sum_ratio': float(np.sum(tri_vol[valid]) / np.sum(b2_vol[covered])),
        'median_ratio': float(np.median(ratio)),
        'p10_ratio': float(np.percentile(ratio, 10)),
        'p90_ratio': float(np.percentile(ratio, 90)),
        'max_abs_relative_error': float(np.max(np.abs(ratio - 1.0))),
        'note': 'sum of mapped triangle volumes divided by the corresponding B2 cell volume',
    }


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
    if not HAS_MATPLOTLIB:
        return
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


def plot_signed_field(outdir, triangulation, wall, name, code, ref, xlim, ylim):
    if not HAS_MATPLOTLIB:
        return
    difference = code - ref
    combined = np.concatenate((np.abs(code[np.isfinite(code)]), np.abs(ref[np.isfinite(ref)])))
    vmax = np.nanpercentile(combined, 99) if combined.size else 1.0
    vmax = max(float(vmax), 1.0e-30)
    field_norm = colors.SymLogNorm(
        linthresh=max(vmax * 1.0e-3, 1.0e-30), vmin=-vmax, vmax=vmax
    )
    finite_difference = np.abs(difference[np.isfinite(difference)])
    difference_max = np.nanpercentile(finite_difference, 99) if finite_difference.size else 1.0
    difference_max = max(float(difference_max), 1.0e-30)
    difference_norm = colors.SymLogNorm(
        linthresh=max(difference_max * 1.0e-3, 1.0e-30),
        vmin=-difference_max,
        vmax=difference_max,
    )
    fig, ax = plt.subplots(1, 3, figsize=(15, 6), dpi=180)
    panels = [
        ('code', code, field_norm),
        ('reference', ref, field_norm),
        ('code-reference', difference, difference_norm),
    ]
    for axis, (title, data, norm) in zip(ax, panels):
        tc = axis.tripcolor(
            triangulation,
            facecolors=data,
            cmap='coolwarm',
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


def plot_b2_field(outdir, name, code, ref):
    if not HAS_MATPLOTLIB:
        return
    ratio = np.full_like(code, np.nan, dtype=float)
    ok = np.isfinite(code) & np.isfinite(ref) & (code > 0.0) & (ref > 0.0)
    ratio[ok] = code[ok] / ref[ok]

    fig, ax = plt.subplots(1, 3, figsize=(15, 5), dpi=180)
    panels = [
        ('code B2', code, safe_lognorm(code), 'jet'),
        ('reference B2', ref, safe_lognorm(ref), 'jet'),
        ('code/reference', ratio, colors.LogNorm(vmin=1e-2, vmax=1e2), 'coolwarm'),
    ]
    for axis, (title, data, norm, cmap) in zip(ax, panels):
        im = axis.imshow(data, origin='lower', aspect='auto', cmap=cmap, norm=norm)
        axis.set_title(name + ' ' + title)
        axis.set_xlabel('radial index')
        axis.set_ylabel('poloidal index')
        fig.colorbar(im, ax=axis, fraction=0.046, pad=0.04)
    fig.savefig(outdir / (name + '_b2_compare.pdf'))
    plt.close(fig)


def plot_target_profiles(outdir, target_profiles):
    if not HAS_MATPLOTLIB or not target_profiles:
        return
    by_field = {}
    for record in target_profiles:
        by_field.setdefault(record['field'], []).append(record)
    output = outdir / 'target_profiles_compare.pdf'
    with PdfPages(output) as pdf:
        for field in sorted(by_field):
            records = by_field[field]
            sides = [
                side for side in ('inner', 'outer')
                if any(record['side'] == side for record in records)
            ]
            if not sides:
                continue
            fig, axes = plt.subplots(
                2, len(sides), figsize=(7.2 * len(sides), 7.5),
                dpi=180, squeeze=False,
            )
            for column, side in enumerate(sides):
                side_records = sorted(
                    (record for record in records if record['side'] == side),
                    key=lambda record: record['radial_index'],
                )
                radial = np.asarray([record['radial_index'] for record in side_records])
                code = np.asarray([record['code'] for record in side_records])
                ref = np.asarray([record['reference'] for record in side_records])
                ratio = np.asarray([record['ratio'] for record in side_records])
                upper = axes[0, column]
                if np.all(code >= 0.0) and np.all(ref >= 0.0) and np.any(code > 0.0) and np.any(ref > 0.0):
                    upper.semilogy(radial, code, marker='o', ms=3, label='NT')
                    upper.semilogy(radial, ref, marker='s', ms=3, label='EIRENE')
                else:
                    upper.plot(radial, code, marker='o', ms=3, label='NT')
                    upper.plot(radial, ref, marker='s', ms=3, label='EIRENE')
                upper.set_title(f'{field} - {side} target')
                upper.set_ylabel('value')
                upper.grid(alpha=0.25)
                upper.legend(frameon=False)
                lower = axes[1, column]
                lower.plot(radial, ratio, color='#d1495b', marker='o', ms=3)
                lower.axhline(1.0, color='k', lw=0.8)
                lower.set_xlabel('radial index')
                lower.set_ylabel('NT / EIRENE')
                lower.set_ylim(bottom=0.0)
                lower.grid(alpha=0.25)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print('wrote', output)


def main():
    parser = argparse.ArgumentParser(
        description='Compare current triangle mesh output with SOLPS/EIRENE triangle references.'
    )
    parser.add_argument('--root', default='.', help='Repository root. Default: current directory')
    parser.add_argument('--case', default='case_input/2MW-5e19')
    parser.add_argument('--output', default='Outputfile/10087101')
    parser.add_argument(
        '--eirene-output',
        default=None,
        help='Directory containing fort.44 and fort.46. Default: case/solps_output',
    )
    parser.add_argument(
        '--eirene-input',
        default=None,
        help='EIRENE input.dat used for the reference. Default: eirene-output/input.dat',
    )
    parser.add_argument(
        '--setting',
        default='Inputfile/settingfile/setting_Trimesh_D_5.log',
        help='NT setting file used for switch-parity checks.',
    )
    parser.add_argument('--outdir', default=None, help='Default: output/fig/tri_benchmark')
    parser.add_argument('--xlim', nargs=2, type=float, default=[1.30, 1.80])
    parser.add_argument('--ylim', nargs=2, type=float, default=[-1.12, -0.60])
    parser.add_argument(
        '--skip-b2-sources',
        action='store_true',
        help='Do not compare ionization source terms on the B2 plasma grid.',
    )
    args = parser.parse_args()

    root = Path(args.root)
    case = root / args.case
    eirene_output = (
        case / 'solps_output'
        if args.eirene_output is None
        else root / args.eirene_output
    )
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
    triangulation = mtri.Triangulation(r, z, triangles) if HAS_MATPLOTLIB else None

    wall_path = case / 'shapedata/wall.txt'
    wall = np.loadtxt(wall_path, skiprows=1) if wall_path.exists() else None

    fields = []
    fort44 = eirene_output / 'fort.44'
    fort46 = eirene_output / 'fort.46'
    eirene_triangles = read_eirene_neutral_triangles(fort46)
    eirene_input = (
        eirene_output / 'input.dat'
        if args.eirene_input is None
        else root / args.eirene_input
    )
    setting_path = root / args.setting if args.setting else None
    eirene_log = select_reference_log(eirene_output)
    direct = [
        ('n_D_0', out / 'n_D_0_Tri'),
        ('n_D2_0', out / 'n_D2_0_Tri'),
        ('T_D_0', out / 'T_D_0_Tri'),
        ('T_D2_0', out / 'T_D2_0_Tri'),
    ]
    for name, code_path in direct:
        if code_path.exists():
            fields.append((name, np.loadtxt(code_path), eirene_triangles[name]))

    velocity_fields = []
    for species in ('D', 'D2'):
        code_path = out / f'{species}_0_V_Tri'
        if not code_path.exists():
            continue
        code = np.loadtxt(code_path)
        ref = eirene_triangles[f'V_{species}_0']
        if code.shape != ref.shape or code.shape != (len(vol), 3):
            raise ValueError(f'V_{species}_0 shape mismatch')
        velocity_fields.append((species, code, ref))

    d2ion_path = out / 'n_D2_1_Tri'
    solps_d2ion_path = case / '2D_data/nDmoleculeion_2D.data'
    if d2ion_path.exists() and solps_d2ion_path.exists():
        fields.append((
            'n_D2_1_mapped',
            np.loadtxt(d2ion_path),
            map_b2_to_tri(read_2d_data(solps_d2ion_path), b2),
        ))

    rows = reference_physics_audit_rows(
        eirene_input, eirene_log, setting_path
    )
    rows.extend(neutral_background_rows(case, eirene_triangles, vol))
    target_profiles = []
    for name, code, ref in fields:
        if len(code) != len(ref) or len(code) != len(vol):
            raise ValueError(name + ' length mismatch')
        row = metric_row(name, code, ref, vol)
        row['mesh'] = 'tri'
        rows.append(row)
        plot_field(outdir, triangulation, wall, name, code, ref, args.xlim, args.ylim)

    field_values = {name: (code, ref) for name, code, ref in fields}
    for species in ('D', 'D2'):
        density_name = f'n_{species}_0'
        temperature_name = f'T_{species}_0'
        if density_name not in field_values or temperature_name not in field_values:
            continue
        code_density, ref_density = field_values[density_name]
        code_temperature, ref_temperature = field_values[temperature_name]
        temperature_rows = density_weighted_temperature_rows(
            'Tri', species, code_density, ref_density,
            code_temperature, ref_temperature, vol,
        )
        for row in temperature_rows:
            row['mesh'] = 'tri'
        rows.extend(temperature_rows)

    inner_tail_triangles = (
        (b2[:, 0] == 1) & (b2[:, 1] >= 30) & (b2[:, 1] <= 36)
    )
    for name in ('n_D2_0', 'T_D2_0'):
        if name not in field_values or not np.any(inner_tail_triangles):
            continue
        code, ref = field_values[name]
        row = metric_row(
            f'Tri_inner_target_r30_36_{name}',
            code[inner_tail_triangles],
            ref[inner_tail_triangles],
            vol[inner_tail_triangles],
        )
        row['mesh'] = 'tri'
        row['region'] = 'inner_target_b2_i1_radial_30_36'
        rows.append(row)

    component_names = ('R', 'Z', 'toroidal')
    for species, code, ref in velocity_fields:
        for component, component_name in enumerate(component_names):
            row = metric_row(
                f'V_{species}_0_{component_name}',
                code[:, component], ref[:, component], vol,
            )
            row['mesh'] = 'tri'
            row['note'] = 'track-length-weighted neutral mean velocity in m/s'
            rows.append(row)
            plot_signed_field(
                outdir, triangulation, wall,
                f'V_{species}_0_{component_name}',
                code[:, component], ref[:, component], args.xlim, args.ylim,
            )

            if np.any(inner_tail_triangles):
                row = metric_row(
                    f'Tri_inner_target_r30_36_V_{species}_0_{component_name}',
                    code[inner_tail_triangles, component],
                    ref[inner_tail_triangles, component],
                    vol[inner_tail_triangles],
                )
                row['mesh'] = 'tri'
                row['region'] = 'inner_target_b2_i1_radial_30_36'
                row['note'] = 'track-length-weighted neutral mean velocity in m/s'
                rows.append(row)

        density_name = f'n_{species}_0'
        if density_name not in field_values:
            continue
        code_density, ref_density = field_values[density_name]
        for component, component_name in enumerate(component_names):
            code_flux = code_density * code[:, component]
            ref_flux = ref_density * ref[:, component]
            row = metric_row(
                f'Particle_flux_{species}_0_{component_name}',
                code_flux, ref_flux, vol,
            )
            row['mesh'] = 'tri'
            row['note'] = 'density times mean velocity in m^-2 s^-1'
            rows.append(row)
            plot_signed_field(
                outdir, triangulation, wall,
                f'Particle_flux_{species}_0_{component_name}',
                code_flux, ref_flux, args.xlim, args.ylim,
            )

            if np.any(inner_tail_triangles):
                row = metric_row(
                    f'Tri_inner_target_r30_36_Particle_flux_{species}_0_{component_name}',
                    code_flux[inner_tail_triangles],
                    ref_flux[inner_tail_triangles],
                    vol[inner_tail_triangles],
                )
                row['mesh'] = 'tri'
                row['region'] = 'inner_target_b2_i1_radial_30_36'
                row['note'] = 'density times mean velocity in m^-2 s^-1'
                rows.append(row)

    if not args.skip_b2_sources:
        vol_path = case / 'vol_2D.data'
        if not vol_path.exists():
            vol_path = case / '2D_data/vol_2D.dat'
        if vol_path.exists():
            b2_vol = read_b2_matrix(vol_path)
            rows.append(b2_tri_volume_row(b2_vol, vol, b2))
            mapped_eirene_b2 = {
                name: map_tri_density_to_b2(
                    eirene_triangles[name], vol, b2, b2_vol.shape
                )[0]
                for name in ('n_D_0', 'n_D2_0')
            }
            b2_neutral_fields = [
                ('B2_n_D_0', 'n_D_0', 'dab2', 1.0),
                ('B2_n_D2_0', 'n_D2_0', 'dmb2', 1.0),
                ('B2_T_D_0', 'T_D_0', 'tab2', 1.0 / 1.602176634e-19),
                ('B2_T_D2_0', 'T_D2_0', 'tmb2', 1.0 / 1.602176634e-19),
            ]
            b2_field_values = {}
            for name, code_name, ref_name, ref_scale in b2_neutral_fields:
                code_path = out / code_name
                if not code_path.exists():
                    continue
                code = read_b2_matrix(code_path)
                ref = read_eirene_b2_field(fort44, ref_name, ref_scale)
                b2_field_values[code_name] = (code, ref)
                row = metric_row(name, code.ravel(), ref.ravel(), b2_vol.ravel())
                row['mesh'] = 'b2'
                row['note'] = 'direct comparison with the exported EIRENE B2 field'
                rows.append(row)
                mapped_name = {
                    'dab2': 'n_D_0',
                    'dmb2': 'n_D2_0',
                }.get(ref_name)
                if mapped_name is not None:
                    mapping_row = metric_row(
                        f'EIRENE_{ref_name}_vs_fort46_mapped',
                        ref.ravel(),
                        mapped_eirene_b2[mapped_name].ravel(),
                        b2_vol.ravel(),
                    )
                    mapping_row['mesh'] = 'eirene_reference'
                    mapping_row['note'] = (
                        f'fort.44 {ref_name} compared with a volume-weighted '
                        f'map of fort.46 {mapped_name}; this validates the '
                        'fine-triangle-to-B2 density mapping'
                    )
                    rows.append(mapping_row)
                integral_name = {
                    'dab2': 'pdena_int_b2',
                    'dmb2': 'pdenm_int_b2',
                }.get(ref_name)
                if integral_name is not None:
                    density_integral = float(np.sum(ref * b2_vol))
                    reported_integral = read_eirene_field(
                        fort44, integral_name
                    )[0]
                    consistency = metric_row(
                        f'EIRENE_{ref_name}_volume_integral_consistency',
                        np.asarray([density_integral]),
                        np.asarray([reported_integral]),
                        np.ones(1),
                    )
                    consistency['mesh'] = 'b2'
                    consistency['note'] = (
                        f'{ref_name} times local vol_2D.dat compared with '
                        f'fort.44 {integral_name} index 0; these are separate '
                        'EIRENE density and source-stratum integral estimators, '
                        'so this is a tally-normalization diagnostic rather than '
                        'a code or mesh-mapping comparison'
                    )
                    rows.append(consistency)
                plot_b2_field(outdir, name, code, ref)

                for side, poloidal in (('inner', 1), ('outer', 96)):
                    for radial_index in range(1, 37):
                        ref_value = ref[poloidal, radial_index]
                        code_value = code[poloidal, radial_index]
                        target_profiles.append({
                            'field': name,
                            'side': side,
                            'radial_index': radial_index,
                            'code': code_value,
                            'reference': ref_value,
                            'ratio': code_value / ref_value if ref_value > 0.0 else np.nan,
                        })

                for region_name, poloidal, first_radial, last_radial in (
                    ('inner_target', 1, 1, 36),
                    ('inner_target_r30_36', 1, 30, 36),
                    ('outer_target', 96, 1, 36),
                ):
                    radial = slice(first_radial, last_radial + 1)
                    region_row = metric_row(
                        f'{name}_{region_name}',
                        code[poloidal, radial],
                        ref[poloidal, radial],
                        b2_vol[poloidal, radial],
                    )
                    region_row['mesh'] = 'b2'
                    region_row['region'] = region_name
                    region_row['note'] = (
                        f'B2 poloidal index {poloidal}, radial indices '
                        f'{first_radial}-{last_radial}'
                    )
                    rows.append(region_row)
            for species in ('D', 'D2'):
                density_name = f'n_{species}_0'
                temperature_name = f'T_{species}_0'
                if (
                    density_name not in b2_field_values
                    or temperature_name not in b2_field_values
                ):
                    continue
                code_density, ref_density = b2_field_values[density_name]
                code_temperature, ref_temperature = b2_field_values[temperature_name]
                temperature_rows = density_weighted_temperature_rows(
                    'B2', species, code_density.ravel(), ref_density.ravel(),
                    code_temperature.ravel(), ref_temperature.ravel(),
                    b2_vol.ravel(),
                )
                for row in temperature_rows:
                    row['mesh'] = 'b2'
                rows.extend(temperature_rows)
            d2_spatial_path = out / 'D2_b2_density_by_source.csv'
            d2_summary_path = out / 'source_stratum_summary.csv'
            d2_density_path = out / 'n_D2_0'
            if (
                d2_spatial_path.exists()
                and d2_summary_path.exists()
                and d2_density_path.exists()
            ):
                rows.extend(b2_source_spatial_rows(
                    d2_spatial_path,
                    d2_summary_path,
                    read_b2_matrix(d2_density_path),
                    b2_vol,
                ))
            code_d2ion_path = out / 'n_D2_1'
            if code_d2ion_path.exists():
                code_d2ion = read_b2_matrix(code_d2ion_path)
                ref_d2ion = read_eirene_b2_field(fort44, 'dib2')
                row = metric_row(
                    'B2_n_D2_1',
                    code_d2ion.ravel(),
                    ref_d2ion.ravel(),
                    b2_vol.ravel(),
                )
                row['mesh'] = 'b2'
                row['note'] = 'molecular-ion density compared directly on the B2 plasma grid'
                rows.append(row)
                plot_b2_field(outdir, 'B2_n_D2_1', code_d2ion, ref_d2ion)
            source_fields = [
                (
                    'B2_Dplus_EI_sources_vs_sum_strata_EI',
                    {
                        'Sn_Ion_D_0': 1.0,
                        'Sn_Diss2_D2_0': 1.0,
                        'Sn_DS1_D2_1': 2.0,
                        'Sn_DS2_D2_1': 1.0,
                    },
                    case / 'sum_strata_ei_sna_Dion.data',
                ),
                ('B2_Sn_Ion_D_0_vs_sum_strata_EI', {'Sn_Ion_D_0': 1.0}, case / 'sum_strata_ei_sna_Dion.data'),
                ('B2_Sn_Ion_D2_0_vs_sum_strata_EI', {'Sn_Ion_D2_0': 1.0}, case / 'sum_strata_ei_sna_Dion.data'),
                (
                    'B2_Dplus_EI_sources_vs_IT_EI',
                    {
                        'Sn_Ion_D_0': 1.0,
                        'Sn_Diss2_D2_0': 1.0,
                        'Sn_DS1_D2_1': 2.0,
                        'Sn_DS2_D2_1': 1.0,
                    },
                    case / 'it_ei_sna_Dion.data',
                ),
                (
                    'B2_Dplus_EI_sources_vs_OT_EI',
                    {
                        'Sn_Ion_D_0': 1.0,
                        'Sn_Diss2_D2_0': 1.0,
                        'Sn_DS1_D2_1': 2.0,
                        'Sn_DS2_D2_1': 1.0,
                    },
                    case / 'ot_ei_sna_Dion.data',
                ),
                (
                    'B2_Dplus_EI_sources_vs_MCW_EI',
                    {
                        'Sn_Ion_D_0': 1.0,
                        'Sn_Diss2_D2_0': 1.0,
                        'Sn_DS1_D2_1': 2.0,
                        'Sn_DS2_D2_1': 1.0,
                    },
                    case / 'mcw_ei_sna_Dion.data',
                ),
            ]
            for name, code_terms, ref_path in source_fields:
                if not ref_path.exists():
                    continue
                code_arrays = [
                    coefficient * read_b2_matrix(out / code_name)
                    for code_name, coefficient in code_terms.items()
                    if (out / code_name).exists()
                ]
                if not code_arrays:
                    continue
                code = np.sum(code_arrays, axis=0)
                ref = read_b2_matrix(ref_path)
                row = source_metric_row(name, code.ravel(), ref.ravel())
                row['mesh'] = 'b2'
                if ref_path.name == 'sum_strata_ei_sna_Dion.data':
                    row['note'] = (
                        'cell-integrated ionization sources in s^-1; summed directly '
                        'without multiplying by cell volume'
                    )
                else:
                    row['note'] = (
                        'diagnostic only: total code D+ EI source is compared with one '
                        'EIRENE source stratum, so the fields are not physically equivalent'
                    )
                rows.append(row)
                plot_b2_field(outdir, name, code, ref)

            d2_sink_paths = [
                out / 'Sn_Ion_D2_0',
                out / 'Sn_Diss1_D2_0',
                out / 'Sn_Diss2_D2_0',
                out / 'Sn_CX_D2_0',
            ]
            if all(path.exists() for path in d2_sink_paths):
                code_sink = np.sum(
                    [read_b2_matrix(path) for path in d2_sink_paths], axis=0
                )
                ref_sink = np.maximum(
                    0.,
                    -read_eirene_b2_field(fort44, 'srcml')
                    * b2_vol * 1.0e6 / 1.602176634e-19,
                )
                code_density_path = out / 'n_D2_0'
                code_density = (
                    read_b2_matrix(code_density_path)
                    if code_density_path.exists() else None
                )
                ref_density = read_eirene_b2_field(fort44, 'dmb2')
                for region_name, region in (
                    ('all', (slice(None), slice(None))),
                    ('inner_target', (1, slice(1, 37))),
                    ('inner_target_r30_36', (1, slice(30, 37))),
                    ('outer_target', (96, slice(1, 37))),
                ):
                    row = source_metric_row(
                        f'B2_D2_sink_vs_srcml_{region_name}',
                        code_sink[region].ravel(),
                        ref_sink[region].ravel(),
                    )
                    row['mesh'] = 'b2'
                    row['region'] = region_name
                    row['note'] = (
                        'code 2.2.9 + 2.2.5g + 2.2.10 + 3.2.3 event rates in s^-1; '
                        'EIRENE -srcml converted from A/cm^3 using B2 cell volume'
                    )
                    if code_density is not None:
                        code_inventory = np.sum(
                            code_density[region] * b2_vol[region]
                        )
                        ref_inventory = np.sum(
                            ref_density[region] * b2_vol[region]
                        )
                        code_loss = np.sum(code_sink[region])
                        ref_loss = np.sum(ref_sink[region])
                        if code_loss > 0. and ref_loss > 0.:
                            code_lifetime = code_inventory / code_loss
                            ref_lifetime = ref_inventory / ref_loss
                            row['code_inventory_over_sink_s'] = code_lifetime
                            row['ref_inventory_over_sink_s'] = ref_lifetime
                            row['inventory_over_sink_ratio'] = (
                                code_lifetime / ref_lifetime
                            )
                    rows.append(row)
                plot_b2_field(
                    outdir, 'B2_D2_sink_vs_srcml', code_sink, ref_sink
                )
        else:
            rows.append({
                'field': 'B2_source_comparison_skipped',
                'mesh': 'b2',
                'note': 'missing vol_2D.data for B2 volume weighting',
            })

    source_summary = out / 'source_stratum_summary.csv'
    if source_summary.exists() and fort44.exists():
        rows.extend(source_stratum_rows(source_summary, fort44))

    target_source_row = inner_target_d2_source_shape_row(
        out, case, fort44, target_profiles
    )
    if target_source_row is not None:
        rows.append(target_source_row)

    if target_profiles:
        with open(outdir / 'b2_target_profiles.csv', 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=target_profiles[0].keys())
            writer.writeheader()
            writer.writerows(target_profiles)
        plot_target_profiles(outdir, target_profiles)

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
    if target_profiles:
        print('wrote', outdir / 'b2_target_profiles.csv')
    if HAS_MATPLOTLIB:
        print('figures in', outdir)
    else:
        print('matplotlib is not available; wrote metrics only')


if __name__ == '__main__':
    main()
