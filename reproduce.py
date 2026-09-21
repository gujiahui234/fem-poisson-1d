"""Reproduce the two 1D experiments, run tests, and compare reference errors."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'results')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    subprocess.run([sys.executable, '-B', '-X', 'utf8', str(ROOT / 'code/run_experiments.py'),
                    '--output-dir', str(output)], check=True)
    subprocess.run([sys.executable, '-B', '-X', 'utf8', str(ROOT / 'code/run_tests.py')], check=True)
    with (ROOT / 'examples/expected_errors.csv').open(encoding='utf-8') as f:
        expected = list(csv.DictReader(f))
    with (output / 'errors.csv').open(encoding='utf-8') as f:
        actual = list(csv.DictReader(f))
    if len(expected) != len(actual):
        raise ValueError('Unexpected number of result rows')
    largest_difference = 0.0
    for reference, row in zip(expected, actual):
        if (row['case'], row['N']) != (reference['case'], reference['N']):
            raise ValueError('Reference cases or meshes do not match')
        for key in ('h', 'E_inf', 'E_L2', 'E_H1sem', 'p_inf', 'p_L2', 'p_H1sem'):
            if not reference[key] and not row[key]:
                continue
            a, b = float(row[key]), float(reference[key])
            if not math.isclose(a, b, rel_tol=2e-6, abs_tol=1e-12):
                raise ValueError(f'Reference mismatch: {row["case"]} N={row["N"]} {key}')
            largest_difference = max(largest_difference, abs(a-b)/abs(b) if b else abs(a-b))
    summary = {'passed': True, 'rows_compared': len(actual),
               'maximum_relative_difference': largest_difference,
               'comparison_rtol': 2e-6, 'comparison_atol': 1e-12,
               'python_version': sys.version.split()[0]}
    (output / 'reproduction_summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print('REPRODUCTION PASSED: 10 rows match the reference within tolerance.')
    print(f'Results: {output}')


if __name__ == '__main__':
    main()
