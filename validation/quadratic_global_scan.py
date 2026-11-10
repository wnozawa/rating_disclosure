import numpy as np,csv
from pathlib import Path
from scipy.interpolate import CubicSpline
from numpy.polynomial.hermite import hermgauss
N=90; nodes,weights=hermgauss(N);z=nodes*np.sqrt(2);w=weights/np.sqrt(np.pi)
R=np.linspace(-9.5,9.5,601)

def one(x,eta,mult):
    lam=1.;V=1/x;g0=1/np.sqrt(2*np.pi*V);L=np.exp(-.5)/(np.sqrt(2*np.pi)*V)
    k=mult*L;t=1-eta;v=np.sqrt(1-eta*eta)
    R2=eta*R[:,None]+v*z[None,:];D=R[:,None]-R2
    def g(d):return g0*np.exp(-.5*d*d/V)
    f=np.sum(w*g(D),axis=1)/k
    for n in range(2000):
        spl=CubicSpline(R,f,bc_type='natural')
        new=np.sum(w*g(D+f[:,None]-spl(np.clip(R2,R[0],R[-1]))),axis=1)/k
        err=np.max(abs(new-f));f=.6*new+.4*f
        if err<1e-13:break
    spl=CubicSpline(R,f,bc_type='natural')
    WN=2*np.sum(w*spl(z)); WF=2*g0/(k*np.sqrt(1+2*x*t))
    return dict(x=x,eta=eta,k_over_bound=mult,W_F=WF,W_N=WN,F_minus_N=WF-WN,relative_gap=(WF-WN)/WF,iterations=n+1,residual=err)
rows=[]
for x in [.5,2,10]:
 for eta in [0,.5,.9]:
  for mult in [1.02,1.5]:
   a=one(x,eta,mult);rows.append(a);print(a,flush=True)
with open(Path(__file__).resolve().parent / 'quadratic_scan_regenerated.csv','w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)
