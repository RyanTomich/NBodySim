#!/usr/bin/env python3
"""Collect validated copies of the inputs used by our completed experiments."""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]


def validate(path, expected):
    with path.open() as source:
        if int(source.readline()) != expected:
            raise ValueError(f'Unexpected body count in {path}')
        positions = set()
        count = 0
        for line in source:
            if not line.strip():
                continue
            row = [float(x) for x in line.split()]
            if len(row) != 10 or not all(math.isfinite(x) for x in row) or row[0] <= 0:
                raise ValueError(f'Invalid body in {path}')
            positions.add(tuple(row[1:4])); count += 1
        if count != expected or len(positions) != count:
            raise ValueError(f'Wrong count or coincident positions in {path}')


def completed(path):
    status = json.loads((path / 'status.json').read_text())
    if status['status'] != 'ok':
        raise ValueError(f'Source experiment is incomplete: {path}')
    return status


def source_location(path):
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('inputs'))
    parser.add_argument('--repository-study', type=Path, default=ROOT / 'results/repository-distributions/20261123T194941133900Z')
    parser.add_argument('--synthetic-studies', type=Path, nargs='+', default=[ROOT / 'results/initial-setups/20261112T210148160514Z', ROOT / 'results/initial-setups/20261114T000327343262Z'])
    parser.add_argument('--spiral-study', type=Path, default=ROOT / 'results/galaxy-barnes-hut/20261123T203250768870Z')
    args = parser.parse_args()
    plan = []
    completed(args.repository_study)
    with (args.repository_study / 'repository-selected.csv').open() as source:
        cases = sorted({(r['distribution'], int(r['bodies'])) for r in csv.DictReader(source)})
    for family, n in cases:
        case = args.repository_study / f'{family}-n{n}'
        metadata = json.loads((case / 'source.json').read_text())
        plan.append((case / 'input.txt', args.output / family / f'n{n}.txt', n, metadata))
    for study in args.synthetic_studies:
        status = completed(study)
        for family in status['settings']['setups']:
            for n in status['settings']['bodies']:
                for seed in status['settings']['seeds']:
                    source = study / f'{family}-n{n}-seed{seed}/input.txt'
                    metadata = dict(source='experiments/initial_setups/generate.py', setup=family, bodies=n, seed=seed,
                        total_mass=1, initial_rms_radius=1, velocity='zero/cold start',
                        physics='unsoftened G=1', source_result=source_location(study), reference_settings=status['settings'])
                    plan.append((source, args.output / ('synthetic_' + family) / f'n{n}_seed{seed}.txt', n, metadata))
    status = completed(args.spiral_study)
    source = args.spiral_study / 'galaxy_movie/input.txt'
    n = status['settings']['bodies']
    metadata = {k: status[k] for k in ('model', 'provenance', 'masses', 'scales', 'model_reference')}
    metadata.update(bodies=n, original_movie_physics='Plummer-softened G=1, epsilon=0.06',
                    suite_physics='same initial state adapted to unsoftened G=1 for all three solvers')
    plan.append((source, args.output / 'spiral_merger' / f'n{n}.txt', n, metadata))
    targets = [target for _, target, _, _ in plan]
    if len(targets) != len(set(targets)):
        raise ValueError('Duplicate target names')
    # Validate every source before copying. Preserve locally edited existing inputs.
    for source, target, n, _ in plan:
        validate(source, n)
        if target.exists() and target.read_bytes() != source.read_bytes():
            raise ValueError(f'Existing input differs; choose a new --output folder: {target}')
    for source, target, _, metadata in plan:
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            shutil.copy2(source, target)
        metadata_path = target.with_suffix('.source.json')
        if not metadata_path.exists():
            metadata_path.write_text(json.dumps(metadata, indent=2) + '\n')
    print(f'Prepared {len(plan)} inputs in {len(set(p.parent.name for p in targets))} families: {args.output.resolve()}')


if __name__ == '__main__':
    main()
