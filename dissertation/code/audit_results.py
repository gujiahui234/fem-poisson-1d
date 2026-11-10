"""Read-only audit of saved FEM results; no assembly or solver rerun."""
from pathlib import Path
import csv
import hashlib
import json
import math
import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
checks = []
def test(name, valid):
    checks.append({'name': name, 'passed': bool(valid)})

def rows(name):
    with (ROOT/'results'/name).open(encoding='utf-8') as f:
        return [{k: (v if k == 'case' or not v else float(v)) for k, v in row.items()} for row in csv.DictReader(f)]

one = rows('one_dimensional.csv')
two = rows('two_dimensional.csv')
sampling = rows('sampling_sensitivity.csv')
refs = rows('reference_convergence.csv')
test('row_counts_and_case_sets', len(one)==10 and len(two)==32 and {r['case'] for r in one+two}=={'sine','power'})
for data, group_keys, measures in [(one,['case'],['E_inf','E_D_inf','E_L2','E_energy']), (two,['case','degree'],['E_L2','E_energy','E_sample'])]:
    previous={}
    for row in data:
        tag=tuple(row[k] for k in group_keys)
        test(f'positive_values_{tag}_{row["n"]}', all(math.isfinite(row[k]) and row[k]>0 for k in measures) and row['h']==1/row['n'])
        if tag in previous:
            old=previous[tag]
            test(f'independent_rates_{tag}_{row["n"]}', all(abs(math.log(old[k]/row[k],2)-row['p_'+k])<1e-12 for k in measures))
        previous[tag]=row

mp.mp.dps=60
for row in one:
    n=int(row['n']); h=mp.mpf(1)/n
    if row['case']=='sine':
        theta=mp.pi/n; B=2*(1-mp.cos(theta))/theta**2
        E2=mp.sqrt((5+mp.cos(theta))/6-B); EE=mp.sqrt(mp.pi**2*(1-B)/2)
    else:
        alpha=mp.mpf(5)/3; l2=mp.mpf(0); energy=mp.mpf(0)
        for i in range(n):
            a=mp.mpf(i)/n;b=mp.mpf(i+1)/n
            slope=(b**alpha-a**alpha)/h; c=a**alpha-slope*a
            D=lambda k:b**k-a**k
            l2+=D(2*alpha+1)/(2*alpha+1)-2*slope*D(alpha+2)/(alpha+2)-2*c*D(alpha+1)/(alpha+1)+slope**2*D(3)/3+slope*c*D(2)+c*c*h
            energy+=alpha**2*D(2*alpha-1)/(2*alpha-1)-slope**2*h
        E2=mp.sqrt(l2);EE=mp.sqrt(energy)
        first=(mp.mpf(2)/5)*(mp.mpf(3)/5)**(mp.mpf(3)/2)*h**alpha
        test(f'power_maximum_identities_{n}', abs(row['E_inf']/float(first)-1)<2e-8 and abs(row['E_D_inf']/float(h**(alpha-1))-1)<2e-8)
    test(f'closed_form_integral_norms_{row["case"]}_{n}', abs(float(E2)-row['E_L2'])<5e-13 and abs(float(EE)-row['E_energy'])<5e-13)

for row in two:
    test(f'grid_rule_and_sampling_{row["case"]}_{row["n"]}_{row["degree"]}', row['n'] in [8,16,32,64] and row['points']=={2:3,3:4,4:6,5:7}[row['degree']] and row['sampling_m']==16)
    if row['degree']==5:
        data=np.load(ROOT/'results'/f'solution2d_{row["case"]}_{int(row["n"])}_d5.npz')
        xy=data['xy']; U=data['U']; Uref=data['Uref']; x,y=xy.T
        exact=(np.sin(np.pi*x) if row['case']=='sine' else x**(5/3)-x)*np.sin(np.pi*y)
        test(f'saved_solution_nodal_diagnostics_{row["case"]}_{row["n"]}',
             abs(np.max(np.abs(U-exact))-row['nodal_gap'])<1e-15 and
             abs(np.max(np.abs(U-Uref))-row['quadrature_solution_gap'])<1e-15)
        boundary=(x==0)|(x==1)|(y==0)|(y==1)
        test(f'zero_dirichlet_values_{row["case"]}_{row["n"]}', np.all(U[boundary]==0) and np.all(Uref[boundary]==0))

for row in sampling:
    test(f'nested_sampling_{row["case"]}_{row["n"]}', row['m16']<=row['m32']<=row['m64'] and abs((row['m64']-row['m16'])/row['m64']-row['relative_16_64'])<1e-15)

validation=json.loads((ROOT/'results/validation.json').read_text())
test('stored_checks_satisfy_their_recorded_limits',len(validation['checks']) in (150,151) and all(c['passed'] and math.isfinite(c['value']) and c['value']<=c['limit'] for c in validation['checks']))
before=json.loads((ROOT/'evidence/step02_before/results_manifest.json').read_text(encoding='utf-8-sig'))
test('saved_experiment_files_unchanged',all(hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256'].lower() for r in before))

summary={
    'one_dimensional_final_rows':[r for r in one if r['n']==128],
    'two_dimensional_final_rows':[r for r in two if r['n']==64 and r['degree']==5],
    'power_2d_sampled_orders':[r['p_E_sample'] for r in two if r['case']=='power' and r['degree']==5 and r['n']>8],
    'power_Q5_to_S16_percent':[{'n':r['n'],'percent':100*r['quadrature_solution_gap']/r['E_sample']} for r in two if r['case']=='power' and r['degree']==5],
    'max_reference12_20_gap':max(r['ref12_ref20_solution_gap'] for r in refs),
    'max_1d_nodal_gap':max(r['nodal_gap'] for r in one),
    'max_1d_residual':max(r['residual_inf'] for r in one),
    'max_2d_residual':max(r['residual_inf'] for r in two),
    'sampling_changes_percent':[{'case':r['case'],'n':r['n'],'percent':100*r['relative_16_64']} for r in sampling],
}
out=ROOT/'qa/results_audit';out.mkdir(exist_ok=True)
report={'passed':all(c['passed'] for c in checks),'checks':checks,'checks_count':len(checks),'summary':summary,'new_FEM_runs':False,'notes':'Closed-form identities and stored-vector diagnostics independently recomputed; previous FEM runs reused.'}
(out/'audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'passed':report['passed'],'checks_count':len(checks),'summary':summary},indent=2))
assert report['passed']
