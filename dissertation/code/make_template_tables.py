"""Compact presentation of saved results; no new numerical solves."""
from pathlib import Path
import csv,math,json
ROOT=Path(__file__).resolve().parents[1]
def read(name):
    return list(csv.DictReader((ROOT/'results'/name).open(encoding='utf-8')))
def sci(v):
    a,b=f'{float(v):.3e}'.split('e')
    return f'${a}\\times10^{{{int(b)}}}$'
def table(name,short,caption,label,headers,rows):
    s=[r'\begin{table}[htbp]',r'\centering\small',f'\\caption[{short}]{{{caption}}}\\label{{{label}}}',
       r'\begin{tabular}{@{}'+'r'*len(headers)+r'@{}}',r'\toprule',' & '.join(headers)+r' \\',r'\midrule']
    s+=[' & '.join(row)+r' \\' for row in rows]
    s += [r'\bottomrule',r'\end{tabular}',r'\end{table}']
    (ROOT/'tables'/f'{name}.tex').write_text('\n'.join(s)+'\n',encoding='utf-8')
one=read('one_dimensional.csv'); two=read('two_dimensional.csv'); sens=read('sampling_sensitivity.csv')
two.sort(key=lambda r: (['sine','power'].index(r['case']), int(r['n']), int(r['degree'])))
table('one_compact','Interval errors for both solutions','Interval errors from assembled systems with analytic loads. All norms are defined in \\eqref{eq:norms}; $h=1/n$.','tab:one',
      ['Case','$n$','$E_\\infty$','$E_{1,\\infty}$','$E_2$','$E_E$'],
      [[r['case'] if r['n']=='8' else '',r['n']]+[sci(r[k]) for k in ['E_inf','E_D_inf','E_L2','E_energy']] for r in one])
table('power_normalized','Normalized integral errors and analytic limits','Power-case integral errors divided by their proved asymptotic powers. The final row gives analytic limits, not another computed mesh.','tab:power-normalized',
      ['$n$','$E_2/h^2$','$E_E/h$'],
      [[r['n'],f"{float(r['E_L2'])/float(r['h'])**2:.6f}",f"{float(r['E_energy'])/float(r['h']):.6f}"] for r in one if r['case']=='power']+
      [['$\\infty$ (limit)','$\\sqrt{10}/18\\approx0.175682$','$5/9\\approx0.555556$']])
table('two_compact','Square-domain errors with degree-five quadrature','Square-domain errors with the degree-five load rule. $S_{16}$ is sampled, and $p_S$ is its observed order before rounding.','tab:two',
      ['Case','$n$','$E_2$','$E_E$','$S_{16}$','$p_S$'],
      [[r['case'] if r['n']=='8' else '',r['n']]+[sci(r[k]) for k in ['E_L2','E_energy','E_sample']]+['--' if not r['p_E_sample'] else f"{float(r['p_E_sample']):.4f}"] for r in two if r['degree']=='5'])
table('quadrature_compact','Power-case quadrature differences','Power-case differences $Q_d$ from the reference solution on the same mesh. Degree denotes polynomial exactness.','tab:quad',
      ['$n$','$Q_2$','$Q_3$','$Q_4$','$Q_5$'],
      [[str(n)]+[sci(next(r['quadrature_solution_gap'] for r in two if r['case']=='power' and int(r['n'])==n and int(r['degree'])==d)) for d in [2,3,4,5]] for n in [8,16,32,64]])
table('nodal_compact','Nodal discrepancy versus quadrature change','Nodal errors $D_5,D_{\\rm ref}$ against the exact solution, compared with the change $Q_5$ between the two discrete solutions.','tab:nodal',
      ['Case','$n$','$D_5$','$D_{\\rm ref}$','$Q_5$'],
      [[r['case'],r['n']]+[sci(r[k]) for k in ['nodal_gap','reference_nodal_gap','quadrature_solution_gap']] for r in two if r['degree']=='5' and r['n'] in ['8','64']])
table('sampling_compact','Sensitivity to the barycentric sampling density','Nested sampling of degree-five solutions. Change is $100(S_{64}-S_{16})/S_{64}$.','tab:sampling',
      ['Case','$n$','$S_{16}$','$S_{32}$','$S_{64}$','Change (\\%)'],
      [[r['case'],r['n']]+[sci(r[k]) for k in ['m16','m32','m64']]+[f"{100*float(r['relative_16_64']):.4f}"] for r in sens])
p=sorted([r for r in two if r['case']=='power' and r['degree']=='5'],key=lambda r:int(r['n']))
bounds=[]
for a,b in zip(p,p[1:]):
    s0,q0=float(a['E_sample']),float(a['quadrature_solution_gap']);s1,q1=float(b['E_sample']),float(b['quadrature_solution_gap'])
    bounds.append([a['n'],b['n'],f"{float(b['p_E_sample']):.4f}",f'{math.log2((s0-q0)/(s1+q1)):.4f}',f'{math.log2((s0+q0)/(s1-q1)):.4f}'])
table('rate_bounds_compact','Bounds for the reference sampled orders','Power-case sampled order and bounds on the reference order on the same sampled sets, from \\eqref{eq:quadrature-rate-bounds}. These do not bound the continuous-error order.','tab:quad-rate',
      ['Coarse $n$','Fine $n$','$p_S$','Lower','Upper'],bounds)
finest={case:max((r for r in one if r['case']==case),key=lambda r:int(r['n'])) for case in ['sine','power']}
summary=[]
for norm,key,sine_theory,power_theory in [
    ('$E_\\infty$','p_E_inf','$O(h^2)$','$C_*h^{5/3}$'),
    ('$E_{1,\\infty}$','p_E_D_inf','$O(h)$','$h^{2/3}$'),
    ('$E_2$','p_E_L2','$\\sim c_s h^2$','$\\sim (\\sqrt{10}/18)h^2$'),
    ('$E_E$','p_E_energy','$\\sim d_s h$','$\\sim (5/9)h$')]:
    summary.append([norm,sine_theory,f"{float(finest['sine'][key]):.5f}",power_theory,f"{float(finest['power'][key]):.5f}"])
table('order_summary','Theory and measured interval orders',
      'Interval theory and observed orders $p$ from $n=64$ to $128$. The smooth maximum bounds are \\eqref{eq:sharpbounds}; exact power maxima are \\eqref{eq:exactpower}; integral equivalents follow from Proposition~\\ref{prp:integral-limits}, with $c_s=\\pi^2/\\sqrt{240}$ and $d_s=\\pi^2/\\sqrt{24}$. All theoretical statements refer to exact Galerkin integration and uniform refinement.',
      'tab:order-summary',['Norm','Sine: theory','$p$','Power: theory','$p$'],summary)
print('Eight tables regenerated from saved CSV data; no FEM solves.')
