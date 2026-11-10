#!/usr/bin/env python3
"""Recompute every numerical display used in the revised Section 4/S7.

Run from any directory:
  OPENBLAS_NUM_THREADS=1 python3 numerics/reproduce.py
Use --render-only to rebuild figures/tables from the shipped CSVs.
"""
import argparse
import csv
import hashlib
import json
import math
import platform
from pathlib import Path
import sys
import numpy as np
import scipy
import matplotlib
from scipy.optimize import brentq
from scipy.special import ndtr
from numpy.polynomial.legendre import leggauss
from solver import Cost, SpectralSolver, independent_effort, kernel, c3

ROOT=Path(__file__).resolve().parent
DATA=ROOT/"data"
ETAS=np.r_[np.linspace(0,.90,47),np.linspace(.91,.995,19)]
CASES={f"E{E}":(Cost(E=float(E)),1.) for E in [3,4,5,6,10,20,100]}
CASES["quadratic"]=(Cost("quadratic"),1.)
CASES.update({f"c{c:g}":(Cost("cubic",c=c),1.) for c in [.01,.03,.1,.3,1.]})
CASES["noisy"]=(Cost(E=6.),.8)


def write_csv(name,rows):
    path=DATA/name
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def read_csv(name):
    with (DATA/name).open() as f:
        return list(csv.DictReader(f))


def root(solver,cost,lam=1.):
    def gap(eta):return solver.solve(cost,float(eta),lam)[0]["gap"]
    return float(brentq(gap,0.,.995,xtol=3e-11,rtol=2e-14))


def outcomes(solver,cost,eta):
    row,e=solver.solve(cost,eta)
    z,w=solver.z,solver.w
    en=solver.avg@e
    ef=cost.inverse(kernel(math.sqrt(2*(1-eta))*z,1.))
    row["U_F"]=1-2*float(w@cost.value(ef))
    row["U_N"]=1-2*float(w@cost.value(en))
    x,wx=leggauss(180)
    sd=math.sqrt(2*(1-eta))
    d=(x+1)*4.5*sd
    wd=wx*4.5*sd*2*kernel(d,sd**2)
    m=math.sqrt((1+eta)/2)*z
    hi=solver.interpolation(m[:,None]+d[None,:]/2)@e
    lo=solver.interpolation(m[:,None]-d[None,:]/2)@e
    row["p_H_N"]=float(np.sum(w[:,None]*wd[None,:]*ndtr(d[None,:]+hi-lo)))
    row["p_H_F"]=.5+math.atan(sd)/math.pi
    return row


def compute(resume=False):
    DATA.mkdir(exist_ok=True)
    s=SpectralSolver(n=161,nq=160,L=10.)
    for key,(cost,lam) in CASES.items():
        filename=f"policy_{key}.csv"
        if resume and (DATA/filename).exists():
            print("Reusing",filename,flush=True)
            continue
        rows=[]
        for eta in ETAS:
            r,_=s.solve(cost,float(eta),lam)
            rows.append(dict(case=key,**r))
        write_csv(filename,rows)
        print("Policy",key,"eta=0 effect",rows[0]["effect_percent"],"max off-grid",
              max(x["offgrid_residual"] for x in rows),flush=True)
    if not resume or not (DATA/"roots.csv").exists():
        rows=[]
        for key,(cost,lam) in CASES.items():
            data=read_csv(f"policy_{key}.csv")
            if float(data[0]["gap"])*float(data[-1]["gap"]) < 0:
                eta=root(s,cost,lam)
                rows.append(dict(case=key,kind=cost.kind,E=cost.E,c=cost.c,k=cost.k,lam=lam,
                                 V=1+2*(1-lam**2),eta_star=eta))
                print("Root",key,eta,flush=True)
        write_csv("roots.csv",rows)
    if not resume or not (DATA/"switching.csv").exists():
        rows=[]
        for E in [5,6,7,10,20,40,60,80,100]:
            eta=root(s,Cost(E=float(E)))
            first=1-math.sqrt(3)/E
            second=first-1.375/E**2
            third=second-c3()/E**3
            rows.append(dict(E=E,eta_star=eta,z_star=E*(1-eta),first=first,second=second,
                             third=third,err1=abs(eta-first),err2=abs(eta-second),err3=abs(eta-third)))
            print("Asymptotic root",E,eta,flush=True)
        write_csv("switching.csv",rows)
    if not resume or not (DATA/"switching_dense.csv").exists():
        # Finer grid for the third-order error display. Its signed error
        # changes sign near E=7, which the coarse grid alone cannot show.
        rows=[]
        for E in np.arange(5.,15.001,.25):
            eta=root(s,Cost(E=float(E)))
            third=1-math.sqrt(3)/E-1.375/E**2-c3()/E**3
            rows.append(dict(E=float(E),eta_star=eta,third=third,diff3=third-eta))
        write_csv("switching_dense.csv",rows)
        print("Dense asymptotic roots",len(rows),flush=True)
    if not resume or not (DATA/"channels.csv").exists():
        rows=[]
        etas=np.unique(np.r_[ETAS[ETAS<=.99],.5,.85,.9,.95,.98,.99])
        for eta in etas:
            r,_=s.solve(Cost(),float(eta))
            r.update(s.channels(Cost(),float(eta)))
            r["remainder"]=r["S"]-r["S2"]
            rows.append(r)
        write_csv("channels.csv",rows)
        print("Channels",len(rows),flush=True)
    if not resume or not (DATA/"matched_gap.csv").exists():
        rows=[]
        for E in [20,40,60,100]:
            for z in np.linspace(.25,3.,45):
                r,_=s.solve(Cost(E=float(E)),float(1-z/E))
                B=8/(2*math.pi)**1.5/.3**3
                rows.append(dict(E=E,z=z,eta=r["eta"],gap=r["gap"],
                                 scaled_gap=E**4*r["gap"]/B,limit=z*z*(3-z*z)))
            print("Matched gap",E,flush=True)
        write_csv("matched_gap.csv",rows)
    if not resume or not (DATA/"cost_levels.csv").exists():
        rows=[]
        for E in [40,80,120]:
            for k in [.27,.3,.4,.6]:
                eta=root(s,Cost(k=k,E=float(E)))
                predicted=1-math.sqrt(3)/E-1.375/E**2-c3(k=k)/E**3
                rows.append(dict(E=E,k=k,eta_star=eta,third=predicted))
            print("Cost-level roots",E,flush=True)
        write_csv("cost_levels.csv",rows)
    if not resume or not (DATA/"outcomes.csv").exists():
        rows=[outcomes(s,cost,eta) for cost,etas in [
            (Cost("quadratic"),[0.,.5]),(Cost(),[0.,.5,.9])] for eta in etas]
        write_csv("outcomes.csv",rows)
        print("Outcomes",len(rows),flush=True)
    validate(s)


def validate(s):
    fine=SpectralSolver(n=193,nq=220,L=12.)
    cases=[("quadratic",0.),("E6",0.),("E6",.5),("E6",.95),("E6",.99),
           ("E10",.785),("c0.1",0.),("c0.1",.3),("noisy",0.),("noisy",.9)]
    rows=[]
    for key,eta in cases:
        cost,lam=CASES[key]
        a,_=s.solve(cost,eta,lam)
        b,_=fine.solve(cost,eta,lam)
        wn,res=independent_effort(cost,eta,lam,ngrid=801,nq=220)
        rows.append(dict(case=key,eta=eta,gap=a["gap"],
                         offgrid_residual=a["offgrid_residual"],
                         resolution_gap_difference=abs(a["gap"]-b["gap"]),
                         independent_effort_difference=abs(a["W_N"]-wn),
                         independent_iteration_residual=res))
    write_csv("validation.csv",rows)
    roots=[]
    for E in [6,20,100]:
        a=root(s,Cost(E=float(E)))
        b=root(fine,Cost(E=float(E)))
        roots.append(dict(E=E,production=a,refined=b,difference=abs(a-b)))
    write_csv("root_validation.csv",roots)
    policy=[r for key in CASES for r in read_csv(f"policy_{key}.csv")]
    channel=read_csv("channels.csv")
    summary=dict(
        production=dict(collocation_nodes=161,hermite_nodes=160,own_rating_limit=10.),
        refined=dict(collocation_nodes=193,hermite_nodes=220,own_rating_limit=12.),
        independent=dict(method="uniform cubic-spline fixed point",rating_nodes=801,hermite_nodes=220),
        policy_evaluations=len(policy),channel_evaluations=len(channel),
        max_collocation_residual=max(float(r["collocation_residual"]) for r in policy),
        max_offgrid_residual=max(float(r["offgrid_residual"]) for r in policy),
        max_resolution_gap_difference=max(r["resolution_gap_difference"] for r in rows),
        max_independent_effort_difference=max(r["independent_effort_difference"] for r in rows),
        max_root_resolution_difference=max(r["difference"] for r in roots),
        max_accounting_error=max(abs(float(r["gap"])-(float(r["J"])-float(r["S"]))) for r in policy),
        checks_passed=True,
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                      matplotlib=matplotlib.__version__),
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in [ROOT/"solver.py",ROOT/"reproduce.py",ROOT/"render.py"]})
    assert summary["max_collocation_residual"]<1e-11,summary
    assert summary["max_offgrid_residual"]<1e-6,summary
    assert summary["max_resolution_gap_difference"]<2e-7,summary
    assert summary["max_independent_effort_difference"]<2e-7,summary
    assert summary["max_root_resolution_difference"]<2e-6,summary
    assert summary["max_accounting_error"]<2e-14,summary
    (DATA/"validation_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--render-only",action="store_true")
    p.add_argument("--resume",action="store_true",help="Reuse complete CSV stages from this unchanged run.")
    args=p.parse_args()
    if not args.render_only:
        compute(args.resume)
    from render import render_all
    render_all()


if __name__=="__main__":
    main()
