"""Audit saved dissertation data, or recompute it in an isolated run directory."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,csv,hashlib,json,math,shutil,subprocess,sys

ROOT=Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recompute',action='store_true',help='assemble and solve all 1D/2D systems before auditing')
    parser.add_argument('--output',type=Path,help='new output directory; existing directories are refused')
    args=parser.parse_args()
    expected=json.loads((ROOT/'reference_manifest.json').read_text())
    assert len(expected)==14 and all(digest(ROOT/r['path'])==r['sha256'] for r in expected),'Reference files changed'
    run=(args.output or ROOT/'runs'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')).resolve()
    run.mkdir(parents=True,exist_ok=False)
    shutil.copytree(ROOT/'code',run/'code',ignore=shutil.ignore_patterns('__pycache__'))
    for name in ['results','qa','tables','figures','logs','evidence/step02_before']:(run/name).mkdir(parents=True,exist_ok=True)
    def execute(name):
        with (run/'logs'/f'{name}.log').open('w',encoding='utf-8') as out:
            result=subprocess.run([sys.executable,'-X','utf8',str(run/'code'/name)],cwd=run,stdout=out,stderr=subprocess.STDOUT)
        if result.returncode:
            print((run/'logs'/f'{name}.log').read_text(encoding='utf-8'))
            raise RuntimeError(f'{name} failed')
        print(f'PASS {name}',flush=True)
    if args.recompute:execute('fem_experiments.py')
    else:
        for p in (ROOT/'reference_results').iterdir():shutil.copy2(p,run/'results'/p.name)
    baseline=[{'path':'results/'+p.name,'sha256':digest(p)} for p in sorted((run/'results').iterdir())]
    (run/'evidence/step02_before/results_manifest.json').write_text(json.dumps(baseline),encoding='utf-8')
    for name in ['audit_results.py','audit_rate_bounds.py','check_integral_limits.py','make_template_tables.py','make_readability_figures.py','make_profile.py']:execute(name)
    checks=[]
    def check(name,ok):
        checks.append({'name':name,'passed':bool(ok)})
    # Compare fresh values with the release reference; CSV row order is preserved.
    # Absolute tolerances accommodate roundoff-sized nodal/reference differences.
    for name in ['one_dimensional.csv','two_dimensional.csv','sampling_sensitivity.csv','reference_convergence.csv']:
        with (ROOT/'reference_results'/name).open(encoding='utf-8') as f:ref=list(csv.DictReader(f))
        with (run/'results'/name).open(encoding='utf-8') as f:actual=list(csv.DictReader(f))
        check(name+'_shape',len(ref)==len(actual) and set(ref[0])==set(actual[0]))
        for i,(a,b) in enumerate(zip(ref,actual)):
            for key,value in a.items():
                if key=='case' or value=='':ok=b[key]==value
                else:ok=math.isclose(float(value),float(b[key]),rel_tol=2e-6,abs_tol=1e-12)
                check(f'{name}:{i}:{key}',ok)
    records={}
    for name in ['results_audit/audit.json','rate_bound_audit.json','integral_limits_audit.json']:
        data=json.loads((run/'qa'/name).read_text())
        check(name,data['passed']);records[name]=len(data['checks'])
    # Recompute the displayed summary orders from raw errors, not stored rate fields.
    with (run/'results/one_dimensional.csv').open() as f:one=list(csv.DictReader(f))
    text=(run/'tables/order_summary.tex').read_text()
    body=text.split(r'\midrule',1)[1].split(r'\bottomrule',1)[0]
    table=[[v.strip() for v in line.removesuffix(r'\\').split('&')] for line in body.splitlines() if '&' in line]
    for cells,key in zip(table,['E_inf','E_D_inf','E_L2','E_energy']):
        for case,column in [('sine',2),('power',4)]:
            a=next(r for r in one if r['case']==case and r['n']=='64');b=next(r for r in one if r['case']==case and r['n']=='128')
            check('summary_'+case+'_'+key,abs(math.log2(float(a[key])/float(b[key]))-float(cells[column]))<=5.1e-6)
    check('outputs_present',all((run/'figures'/f'{n}.{ext}').is_file() for n in ['one_convergence','two_convergence','power_profile'] for ext in ['png','pdf']) and len(list((run/'tables').glob('*.tex')))==8)
    check('results_unchanged_by_audits',all(digest(run/r['path'])==r['sha256'] for r in baseline))
    check('release_reference_unchanged',all(digest(ROOT/r['path'])==r['sha256'] for r in expected))
    validation=json.loads((run/'results/validation.json').read_text())
    summary={'passed':all(c['passed'] for c in checks),'recomputed_FEM':args.recompute,'recorded_experiment_checks':len(validation['checks']),'audit_counts':records,'checks':checks,'python':sys.version.split()[0],'comparison_tolerances':{'relative':2e-6,'absolute':1e-12}}
    (run/'reproduction_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    assert summary['passed'],[c for c in checks if not c['passed']]
    print('REPRODUCTION PASSED; mode='+('recompute' if args.recompute else 'saved-data audit'))
    print('Output:',run)
if __name__=='__main__':main()
