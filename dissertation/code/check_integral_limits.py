"""Validate normalized-error table against saved data and independent identities.

No finite-element system is assembled or solved. The mathematical limit is
proved in the manuscript, not inferred from these finite-mesh checks.
"""
from pathlib import Path
import csv, hashlib, json
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 60
manifest = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (ROOT/'results').iterdir() if p.is_file()}
checks = []
def check(name, ok):
    checks.append({'name': name, 'passed': bool(ok)})
rows = [r for r in csv.DictReader((ROOT/'results/one_dimensional.csv').open()) if r['case']=='power']
c2, ce = mp.sqrt(10)/18, mp.mpf(5)/9
check('quadratic_L2_coefficient', abs(mp.quad(lambda t: t*t*(t-1)**2/4, [0,1])-mp.mpf(1)/120)<mp.mpf('1e-55'))
check('quadratic_energy_coefficient', abs(mp.quad(lambda t: (t-mp.mpf('0.5'))**2, [0,1])-mp.mpf(1)/12)<mp.mpf('1e-55'))
check('power_L2_limit_simplification', abs(c2-mp.sqrt((mp.mpf(100)/27)/120))<mp.mpf('1e-55'))
check('power_energy_limit_simplification', abs(ce-mp.sqrt((mp.mpf(100)/27)/12))<mp.mpf('1e-55'))
normalized = []
for r in rows:
    n=int(r['n']); h=mp.mpf(1)/n; alpha=mp.mpf(5)/3
    norm2=mp.mpf(0); norme=mp.mpf(0)
    for j in range(n):
        a=mp.mpf(j)/n; b=a+h; slope=(b**alpha-a**alpha)/h; c=a**alpha-slope*a
        D=lambda k: b**k-a**k
        norm2 += D(2*alpha+1)/(2*alpha+1)-2*slope*D(alpha+2)/(alpha+2)-2*c*D(alpha+1)/(alpha+1)+slope*slope*D(3)/3+slope*c*D(2)+c*c*h
        norme += alpha*alpha*D(2*alpha-1)/(2*alpha-1)-slope*slope*h
    observed2=mp.mpf(r['E_L2'])/h**2; observede=mp.mpf(r['E_energy'])/h
    reference2=mp.sqrt(norm2)/h**2; referencee=mp.sqrt(norme)/h
    check(f'normalized_L2_independent_identity_n{n}', abs(observed2-reference2)<mp.mpf('1e-9'))
    check(f'normalized_energy_independent_identity_n{n}', abs(observede-referencee)<mp.mpf('1e-10'))
    normalized.append({'n':n,'L2_over_h2':float(observed2),'energy_over_h':float(observede),
                       'L2_below_limit_percent':float(100*(1-observed2/c2)),
                       'energy_below_limit_percent':float(100*(1-observede/ce))})
check('saved_result_files_unchanged', all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in manifest.items()))
lines=[r'\begin{table}[htbp]',r'\centering\small',
       r'\caption{Power-case integral errors normalized by their proved asymptotic powers. Finite rows use the saved assembled-system errors; the last row contains the analytic limits.}\label{tab:power-normalized}',
       r'\begin{tabular}{rrr}',r'\toprule',r'$n$ & $E_2/h^2$ & $E_E/h$ \\',r'\midrule']
for r in normalized:
    lines.append(f"{r['n']} & {r['L2_over_h2']:.6f} & {r['energy_over_h']:.6f} "+r'\\')
lines += [r'\midrule',r'$\infty$ (limit) & $\sqrt{10}/18\approx0.175682$ & $5/9\approx0.555556$ \\',r'\bottomrule',r'\end{tabular}',r'\end{table}']
(ROOT/'tables/power_normalized.tex').write_text('\n'.join(lines)+'\n',encoding='utf-8')
report={'passed':all(c['passed'] for c in checks),'checks':checks,'rows':normalized,
        'limits':{'L2_over_h2':float(c2),'energy_over_h':float(ce)},'fem_solves_performed':0}
(ROOT/'qa/integral_limits_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2))
assert report['passed']
