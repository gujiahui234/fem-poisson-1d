"""Independent dissertation experiments. Writes only to this new project.

Adapted from the project's earlier 1D/2D implementations; singular-source
reference quadrature, degree-labelled triangle rules and checks are new here.
"""
from pathlib import Path
import csv,json,math,platform,time
from itertools import permutations
import numpy as np
import scipy
from scipy.linalg import solve_banded
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.optimize import brentq
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
OUT.mkdir(exist_ok=True)
CHECKS=[]

def check(name,value,limit):
    passed=bool(np.isfinite(value) and value<=limit)
    CHECKS.append(dict(name=name,value=float(value),limit=float(limit),passed=passed))
    if not passed:raise AssertionError(CHECKS[-1])

def write_csv(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def add_rates(rows,group,keys):
    previous={}
    for r in rows:
        tag=tuple(r[k] for k in group)
        prev=previous.get(tag)
        for k in keys:r['p_'+k]=math.log(prev[k]/r[k])/math.log(prev['h']/r['h']) if prev else ''
        previous[tag]=r

def rule(degree):
    """Symmetric barycentric rules; weights sum to one (multiply by area)."""
    lam=[];weights=[]
    def orbit(a,b,w):
        p=sorted(set(permutations((a,b,b))))
        lam.extend(p);weights.extend([w]*len(p))
    if degree==2:orbit(2/3,1/6,1/3)
    elif degree==3:
        lam.append((1/3,)*3);weights.append(-27/48)
        orbit(0.6,0.2,25/48)
    elif degree==4:
        orbit(0.108103018168070,0.445948490915965,0.223381589678011)
        orbit(0.816847572980459,0.091576213509771,0.109951743655322)
    elif degree==5:
        lam.append((1/3,)*3);weights.append(0.225)
        orbit(0.059715871789770,0.470142064105115,0.132394152788506)
        orbit(0.797426985353087,0.101286507323456,0.125939180544827)
    else:raise ValueError(degree)
    return np.asarray(lam),np.asarray(weights)

def check_rules():
    specs=[]
    for d in range(2,6):
        lam,w=rule(d);err=0.
        for i in range(d+1):
            for j in range(d+1-i):
                exact=2*math.factorial(i)*math.factorial(j)/math.factorial(i+j+2)
                err=max(err,abs(w@(lam[:,1]**i*lam[:,2]**j)-exact))
        check(f'quadrature_degree_{d}_all_moments',err,3e-14)
        specs.append(dict(degree=d,points=len(w),barycentric=lam.tolist(),weights_normalized_to_one=w.tolist(),max_moment_error=err))
    (OUT/'quadrature_rules.json').write_text(json.dumps(specs,indent=2),encoding='utf-8')

def functions1(case):
    if case=='sine':return (lambda x:np.sin(np.pi*x),lambda x:np.pi*np.cos(np.pi*x),lambda x:np.pi**2*np.sin(np.pi*x))
    if case=='power':return (lambda x:x**(5/3)-x,lambda x:(5/3)*x**(2/3)-1,lambda x:-(10/9)*x**(-1/3))
    return (lambda x:x*(1-x)/2,lambda x:0.5-x,lambda x:np.ones_like(x))

def load1(case,n):
    """Analytic element loads, evaluated in 60-digit arithmetic to avoid cancellation."""
    import mpmath as mp
    mp.mp.dps=60
    F0={'power':lambda x:-mp.mpf(5)/3*x**(mp.mpf(2)/3),'sine':lambda x:-mp.pi*mp.cos(mp.pi*x),'quadratic':lambda x:x}[case]
    F1={'power':lambda x:-mp.mpf(2)/3*x**(mp.mpf(5)/3),'sine':lambda x:-mp.pi*x*mp.cos(mp.pi*x)+mp.sin(mp.pi*x),'quadratic':lambda x:x*x/2}[case]
    b=np.zeros(n+1);local=[]
    for i in range(n):
        a=mp.mpf(i)/n;c=mp.mpf(i+1)/n;h=c-a
        j0=F0(c)-F0(a);j1=F1(c)-F1(a)
        pair=np.array([float((c*j0-j1)/h),float((j1-a*j0)/h)])
        b[i:i+2]+=pair;local.append(pair)
    return b,np.asarray(local)

def measure1(case,n,U,q=32):
    u,du,_=functions1(case);x=np.linspace(0,1,n+1)
    z,w=leggauss(q);z=(z+1)/2;w=w/2
    e2=ee=em=dm=0.
    for i in range(n):
        a,b=x[i:i+2];h=b-a;slope=(U[i+1]-U[i])/h
        candidates=[a,b]
        if du(a)-slope>0 and du(b)-slope<0 or du(a)-slope<0 and du(b)-slope>0:
            candidates.append(brentq(lambda t:du(t)-slope,a,b,xtol=5e-15))
        em=max(em,max(abs(u(t)-U[i]-slope*(t-a)) for t in candidates))
        dm=max(dm,abs(du(a)-slope),abs(du(b)-slope))
        if case=='power' and i==0:xx=h*z**3;ww=3*h*w*z**2
        else:xx=a+h*z;ww=h*w
        e2+=ww@((u(xx)-U[i]-slope*(xx-a))**2)
        ee+=ww@((du(xx)-slope)**2)
    return dict(E_inf=em,E_D_inf=dm,E_L2=math.sqrt(e2),E_energy=math.sqrt(ee))

def run1():
    rows=[]
    for case in ['sine','power']:
        u,du,f=functions1(case)
        for n in [8,16,32,64,128]:
            b,local=load1(case,n);ab=np.zeros((3,n-1));ab[0,1:]=-n;ab[1]=2*n;ab[2,:-1]=-n
            U=np.zeros(n+1);U[1:-1]=solve_banded((1,1),ab,b[1:-1]);x=np.linspace(0,1,n+1)
            vals=measure1(case,n,U);v64=measure1(case,n,U,64)
            gap=max(abs(U-u(x)));res=max(abs(2*n*U[1:-1]-n*U[:-2]-n*U[2:]-b[1:-1]))
            check(f'1d_{case}_{n}_nodal',gap,3e-12)
            check(f'1d_{case}_{n}_error_integrals',max(abs(vals[k]-v64[k]) for k in vals),3e-12)
            # Independent scipy adaptive check of the first and last element loads.
            for i in [0,n-1]:
                for side in [0,1]:
                    def integrand(t):
                        if t==0:return 0.
                        r=t**3;xx=(i+r)/n
                        return f(xx)*(r if side else 1-r)*3*t*t/n
                    ref=quad(integrand,0,1,epsabs=2e-13,epsrel=2e-13)[0]
                    check(f'1d_{case}_{n}_load_{i}_{side}',abs(ref-local[i,side]),1e-12)
            row=dict(case=case,n=n,h=1/n,**vals,nodal_gap=gap,residual_inf=res,first_E_inf='',first_E_D_inf='')
            if case=='power':
                row['first_E_inf']=(2/5)*(3/5)**1.5*n**(-5/3);row['first_E_D_inf']=n**(-2/3)
                check(f'power_{n}_analytic_max',abs(vals['E_inf']/row['first_E_inf']-1),2e-8)
                check(f'power_{n}_analytic_derivative',abs(vals['E_D_inf']/row['first_E_D_inf']-1),2e-8)
                if n in [8,128]:
                    import mpmath as mp
                    mp.mp.dps=60;alpha=mp.mpf(5)/3;l2=mp.mpf(0);en=mp.mpf(0)
                    for i in range(n):
                        a=mp.mpf(i)/n;c=mp.mpf(i+1)/n;hh=c-a
                        slope=(c**alpha-a**alpha)/hh;intercept=a**alpha-slope*a
                        D=lambda k:c**k-a**k
                        l2+=D(2*alpha+1)/(2*alpha+1)-2*slope*D(alpha+2)/(alpha+2)-2*intercept*D(alpha+1)/(alpha+1)+slope**2*D(3)/3+slope*intercept*D(2)+intercept**2*hh
                        en+=alpha**2*D(2*alpha-1)/(2*alpha-1)-slope**2*hh
                    check(f'power_{n}_independent_L2',abs(vals['E_L2']-float(mp.sqrt(l2))),5e-13)
                    check(f'power_{n}_independent_energy',abs(vals['E_energy']-float(mp.sqrt(en))),5e-13)
            else:
                import mpmath as mp
                mp.mp.dps=60;t=mp.pi/n;S=2*(1-mp.cos(t))/t**2
                exactL=float(mp.sqrt((5+mp.cos(t))/6-S));exactE=float(mp.sqrt(mp.pi**2*(1-S)/2))
                check(f'sine_{n}_independent_L2',abs(vals['E_L2']-exactL),5e-13)
                check(f'sine_{n}_independent_energy',abs(vals['E_energy']-exactE),5e-13)
            rows.append(row)
        print('1D complete',case,flush=True)
    add_rates(rows,['case'],['E_inf','E_D_inf','E_L2','E_energy']);write_csv('one_dimensional.csv',rows)
    b,_=load1('quadratic',2);U=np.array([0,b[1]/4,0]);val=measure1('quadratic',2,U)
    check('two_element_U',abs(U[1]-1/8),1e-15);check('two_element_max',abs(val['E_inf']-1/32),1e-15)
    return rows

def setup2(n):
    xy=np.array([(i/n,j/n) for j in range(n+1) for i in range(n+1)])
    cells=np.array([j*(n+1)+i for j in range(n) for i in range(n)])
    tri=np.empty((2*n*n,3),int);tri[0::2]=np.stack([cells,cells+1,cells+n+2],axis=-1);tri[1::2]=np.stack([cells,cells+n+2,cells+n+1],axis=-1)
    pts=xy[tri];M=np.concatenate([np.ones((len(tri),3,1)),pts],axis=2);gr=np.linalg.inv(M)[:,1:,:].transpose(0,2,1);area=np.linalg.det(M)/2
    ke=area[:,None,None]*np.einsum('tix,tjx->tij',gr,gr)
    A=coo_matrix((ke.ravel(),(np.repeat(tri,3,axis=1).ravel(),np.tile(tri,(1,3)).ravel())),shape=(len(xy),len(xy))).tocsr()
    interior=np.flatnonzero(np.all((xy>0)&(xy<1),axis=1));Ai=A[interior][:,interior].tocsc()
    return dict(n=n,xy=xy,tri=tri,pts=pts,gr=gr,area=area,A=A,Ai=Ai,interior=interior,lu=splu(Ai))

def exact2(p,case):
    x,y=p[...,0],p[...,1]
    return (np.sin(np.pi*x) if case=='sine' else x**(5/3)-x)*np.sin(np.pi*y)

def gradient2(p,case):
    x,y=p[...,0],p[...,1]
    v=np.sin(np.pi*x) if case=='sine' else x**(5/3)-x
    dv=np.pi*np.cos(np.pi*x) if case=='sine' else (5/3)*x**(2/3)-1
    return np.stack([dv*np.sin(np.pi*y),np.pi*v*np.cos(np.pi*y)],axis=-1)

def force2(p,case):
    x,y=p[...,0],p[...,1]
    if case=='sine':return 2*np.pi**2*exact2(p,case)
    return (-(10/9)*x**(-1/3)+np.pi**2*(x**(5/3)-x))*np.sin(np.pi*y)

def ref_rule2(indices,n,q):
    """Cartesian vertical-section integration; x=h*t^3 in the first strip."""
    z,w=leggauss(q);z=(z+1)/2;w=w/2
    r,s=np.meshgrid(z,z,indexing='ij');ww=np.outer(w,w).ravel();r=r.ravel();s=s.ravel()
    cell=indices//2;i=cell%n;lower=(indices%2==0)
    rr=np.broadcast_to(r,(len(indices),len(r))).copy();dr=np.ones_like(rr)
    singular=(i==0);rr[singular]=r**3;dr[singular]=3*r*r
    lam=np.empty((len(indices),len(r),3))
    lam[lower]=np.stack([1-rr[lower],rr[lower]*(1-s),rr[lower]*s],axis=-1)
    lam[~lower]=np.stack([(1-rr[~lower])*(1-s),rr[~lower],(1-rr[~lower])*s],axis=-1)
    jac=np.where(lower[:,None],rr,1-rr)*dr/n**2
    return lam,jac*ww

def load2(m,case,degree=None,qref=16):
    b=np.zeros(len(m['xy']));tri=m['tri']
    if degree is not None:
        lam,w=rule(degree);X=np.einsum('qi,tij->tqj',lam,m['pts']);F=force2(X,case)
        loc=m['area'][:,None]*np.einsum('qi,tq,q->ti',lam,F,w)
        np.add.at(b,tri.ravel(),loc.ravel())
    else:
        for begin in range(0,len(tri),256):
            ids=np.arange(begin,min(begin+256,len(tri)));lam,w=ref_rule2(ids,m['n'],qref)
            X=np.einsum('tqi,tij->tqj',lam,m['pts'][ids]);loc=np.einsum('tqi,tq,tq->ti',lam,force2(X,case),w)
            np.add.at(b,tri[ids].ravel(),loc.ravel())
    return b

def solve2(m,b):
    U=np.zeros(len(m['xy']));U[m['interior']]=m['lu'].solve(b[m['interior']]);return U

def sample2(m,U,case,subdivisions=16):
    sm=subdivisions;lam=np.array([(1-(i+j)/sm,i/sm,j/sm) for i in range(sm+1) for j in range(sm+1-i)])
    value=0.
    for begin in range(0,len(m['tri']),128):
        tr=m['tri'][begin:begin+128];pts=m['pts'][begin:begin+128]
        X=np.einsum('qi,tij->tqj',lam,pts)
        value=max(value,float(np.max(abs(exact2(X,case)-np.einsum('qi,ti->tq',lam,U[tr])))))
    return value

def integral_errors2(m,U,case,q=16):
    e2=ee=0.;tri=m['tri'];grad=np.einsum('ti,tij->tj',U[tri],m['gr'])
    for begin in range(0,len(tri),256):
        ids=np.arange(begin,min(begin+256,len(tri)));lam,w=ref_rule2(ids,m['n'],q)
        X=np.einsum('tqi,tij->tqj',lam,m['pts'][ids]);approx=np.einsum('tqi,ti->tq',lam,U[tri[ids]])
        e=exact2(X,case)-approx;de=gradient2(X,case)-grad[ids,None,:]
        e2+=np.sum(w*e*e);ee+=np.sum(w*np.sum(de*de,axis=-1))
    return math.sqrt(e2),math.sqrt(ee)

def validate2(m):
    # A separate affine patch test uses nonzero boundary data.
    xy=m['xy'];u=1+2*xy[:,0]-3*xy[:,1];r=m['A']@u
    check(f'2d_{m["n"]}_affine_patch',max(abs(r[m['interior']])),2e-13)
    check(f'2d_{m["n"]}_total_area',abs(sum(m['area'])-1),3e-12)

def run2():
    rows=[];sens=[];refs=[]
    for n in [8,16,32,64]:
        m=setup2(n);validate2(m)
        for case in ['sine','power']:
            bref=load2(m,case,qref=20);b12=load2(m,case,qref=12)
            Uref=solve2(m,bref);U12=solve2(m,b12);refgap=max(abs(Uref-U12))
            check(f'2d_{case}_{n}_ref12_ref20',refgap,2e-11)
            nodalref=max(abs(Uref-exact2(m['xy'],case)))
            refs.append(dict(case=case,n=n,ref12_ref20_solution_gap=refgap,reference_nodal_gap=nodalref))
            for degree in range(2,6):
                b=load2(m,case,degree=degree);U=solve2(m,b);E2,EE=integral_errors2(m,U,case,q=20)
                maximum=sample2(m,U,case,16);delta=max(abs(U-Uref));res=max(abs(m['Ai']@U[m['interior']]-b[m['interior']]))
                check(f'2d_{case}_{n}_d{degree}_residual',res,5e-13)
                rows.append(dict(case=case,n=n,h=1/n,degree=degree,points=len(rule(degree)[1]),E_L2=E2,E_energy=EE,E_sample=maximum,nodal_gap=float(max(abs(U-exact2(m['xy'],case)))),quadrature_solution_gap=float(delta),reference_nodal_gap=float(nodalref),reference_12_20_gap=float(refgap),residual_inf=float(res),sampling_m=16))
                if degree==5:
                    l12,e12=integral_errors2(m,U,case,q=12)
                    check(f'2d_{case}_{n}_norm_integral12_20',max(abs(l12-E2),abs(e12-EE)),3e-11)
                    if n in [8,64]:
                        s32=sample2(m,U,case,32);s64=sample2(m,U,case,64)
                        sens.append(dict(case=case,n=n,m16=maximum,m32=s32,m64=s64,relative_16_64=(s64-maximum)/s64))
                        check(f'2d_{case}_{n}_nested_samples',max(0,maximum-s32,s32-s64),2e-14)
                    np.savez_compressed(OUT/f'solution2d_{case}_{n}_d5.npz',xy=m['xy'],tri=m['tri'],U=U,Uref=Uref)
            print('2D complete',case,n,flush=True)
    add_rates(rows,['case','degree'],['E_L2','E_energy','E_sample']);write_csv('two_dimensional.csv',rows);write_csv('sampling_sensitivity.csv',sens);write_csv('reference_convergence.csv',refs)
    return rows

def main():
    started=time.time();check_rules();one=run1();two=run2()
    result=dict(passed=all(r['passed'] for r in CHECKS),checks=CHECKS,elapsed_seconds=time.time()-started,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,rows_1d=len(one),rows_2d=len(two),quadrature_interpretation='2 through 5 means polynomial exactness degree, following the supplied problem sheet; point counts are 3,4,6,7.',scope='1D primary; 2D extension. No new general 2D theorem proved.')
    (OUT/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('PASS',len(CHECKS),'checks; seconds',round(result['elapsed_seconds'],2),flush=True)

if __name__=='__main__':main()
