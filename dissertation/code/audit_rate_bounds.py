"""Check the new rate intervals from saved vectors; no assembly or solve."""
from pathlib import Path
import csv, hashlib, json, math
import numpy as np
from fem_experiments import sample2

ROOT = Path(__file__).resolve().parents[1]
paths = sorted((ROOT / 'results').glob('*'))
before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()}
rows = [r for r in csv.DictReader((ROOT / 'results/two_dimensional.csv').open())
        if r['case'] == 'power' and int(r['degree']) == 5]
rows.sort(key=lambda r: int(r['n']))
samples = []
checks = []
for r in rows:
    n = int(r['n'])
    z = np.load(ROOT / f'results/solution2d_power_{n}_d5.npz')
    m = {'tri': z['tri'], 'pts': z['xy'][z['tri']]}
    s = sample2(m, z['Uref'], 'power', 16)
    q = float(r['quadrature_solution_gap'])
    s5 = float(r['E_sample'])
    checks.append({'name': f'sampled_difference_bound_n{n}',
                   'passed': abs(s - s5) <= q + 2e-15})
    samples.append({'n': n, 'S_reference': s, 'S5': s5, 'Q5': q})
intervals = []
for a, b in zip(samples, samples[1:]):
    lo = math.log2((a['S5'] - a['Q5']) / (b['S5'] + b['Q5']))
    hi = math.log2((a['S5'] + a['Q5']) / (b['S5'] - b['Q5']))
    rate = math.log2(a['S_reference'] / b['S_reference'])
    checks.append({'name': f'reference_rate_in_interval_{a["n"]}_{b["n"]}',
                   'passed': lo <= rate <= hi})
    intervals.append({'coarse_n': a['n'], 'fine_n': b['n'],
                      'lower': lo, 'reference_rate': rate, 'upper': hi})
checks.append({'name': 'all_result_files_unchanged',
               'passed': all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == h for p, h in before.items())})
report = {'passed': all(c['passed'] for c in checks), 'checks': checks,
          'samples': samples, 'intervals': intervals, 'fem_solves_performed': 0}
(ROOT / 'qa/rate_bound_audit.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
print(json.dumps(report, indent=2))
assert report['passed']
