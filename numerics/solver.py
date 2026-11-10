"""Gaussian rating-contest equilibrium and local feedback channels.

The maintained primitives imply V = performance_variance + 2*(1-lam**2).
All integration weights below integrate N(0,1), not exp(-x*x).
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.special import roots_hermitenorm
from scipy.linalg import solve
from scipy.interpolate import CubicSpline


def normal_rule(n):
    x, w = roots_hermitenorm(n)
    return x, w / math.sqrt(2 * math.pi)


def kernel(x, V):
    return np.exp(-np.asarray(x)**2/(2*V))/math.sqrt(2*math.pi*V)


@dataclass(frozen=True)
class Cost:
    kind: str = "saturating"
    k: float = .3
    E: float = 6.
    c: float = 0.

    def marginal(self, e):
        if self.kind == "saturating":
            return self.k*e/(1-(e/self.E)**2)**2
        return self.k*e + self.c*e**2/2

    def second(self, e):
        if self.kind == "saturating":
            x=(e/self.E)**2
            return self.k*(1+3*x)/(1-x)**3
        return self.k+self.c*e

    def third(self, e):
        if self.kind == "saturating":
            x=(e/self.E)**2
            return 12*self.k*e/self.E**2*(1+x)/(1-x)**4
        return np.zeros_like(e)+self.c

    def value(self, e):
        if self.kind == "saturating":
            return self.k*e**2/(2*(1-(e/self.E)**2))
        return self.k*e**2/2+self.c*e**3/6

    def inverse(self, y):
        y=np.asarray(y)
        if self.kind != "saturating":
            return y/self.k if self.c == 0 else 2*y/(self.k+np.sqrt(self.k**2+2*self.c*y))
        e=np.minimum(y/self.k, .85*self.E)
        for _ in range(30):
            step=(self.marginal(e)-y)/self.second(e)
            e=np.clip(e-step,0,self.E*(1-1e-13))
            if np.max(np.abs(step)) < 2e-15:
                break
        return e


class SpectralSolver:
    """Chebyshev collocation, conditional Hermite quadrature, damped Newton.

    Opponent efforts use barycentric polynomial interpolation. Conditional
    rating arguments outside [-L,L] are clipped. Increasing L is included
    among the published resolution checks.
    """
    def __init__(self, n=97, nq=120, L=10.):
        self.n,self.nq,self.L=n,nq,L
        j=np.arange(n)
        self.r=L*np.cos(np.pi*j/(n-1))
        self.bw=(-1.)**j
        self.bw[[0,-1]]*=.5
        self.z,self.w=normal_rule(nq)
        self.avg=self.interpolation(self.z)

    def interpolation(self, x):
        shape=np.shape(x)
        x=np.clip(np.asarray(x).ravel(),-self.L,self.L)
        d=x[:,None]-self.r[None,:]
        hit=np.abs(d)<2e-13
        d[hit]=1.
        a=self.bw[None,:]/d
        a/=a.sum(axis=1)[:,None]
        rows=np.where(hit.any(axis=1))[0]
        a[rows]=hit[rows].astype(float)
        return a.reshape(shape+(self.n,))

    def solve(self,cost,eta,lam=1.,performance_variance=1.):
        V=performance_variance+2*(1-lam**2)
        assert 0 <= eta < 1 and 0 < lam <= 1 and V > 0
        assert cost.k > math.exp(-.5)/math.sqrt(2*math.pi)/V
        oth=eta*self.r[:,None]+math.sqrt(1-eta**2)*self.z[None,:]
        P=self.interpolation(oth)
        Z=lam*(self.r[:,None]-oth)
        b=kernel(lam*(1-eta)*self.r,V+lam**2*(1-eta**2))
        e=cost.inverse(b)
        for it in range(35):
            arg=Z+e[:,None]-np.einsum("iqj,j->iq",P,e)
            gg=kernel(arg,V)
            gp=-arg/V*gg
            F=cost.marginal(e)-gg@self.w
            err=float(np.max(np.abs(F)))
            if err < 3e-14:
                break
            A=np.diag(cost.second(e)-gp@self.w)+np.einsum("iq,q,iqj->ij",gp,self.w,P)
            step=solve(A,F,check_finite=False)
            alpha=1.
            for _ in range(22):
                trial=e-alpha*step
                if np.min(trial)>=0 and (cost.kind!="saturating" or np.max(trial)<cost.E):
                    argt=Z+trial[:,None]-np.einsum("iqj,j->iq",P,trial)
                    Ft=cost.marginal(trial)-kernel(argt,V)@self.w
                    if np.max(np.abs(Ft))<err:
                        break
                alpha*=.5
            e=trial
        else:
            raise RuntimeError(f"Newton failed: eta={eta}, residual={err}")
        en=self.avg@e
        WN=2*float(self.w@en)
        WF=2*float(self.w@cost.inverse(kernel(lam*math.sqrt(2*(1-eta))*self.z,V)))
        Q=2*float(self.w@cost.inverse(kernel(lam*(1-eta)*self.z,V+lam**2*(1-eta**2))))
        # FOC residual at sampled, non-collocation own ratings.
        rr=np.linspace(-5.5,5.5,143)
        ee=self.interpolation(rr)@e
        other=eta*rr[:,None]+math.sqrt(1-eta**2)*self.z[None,:]
        pe=np.einsum("iqj,j->iq",self.interpolation(other),e)
        off=cost.marginal(ee)-kernel(lam*(rr[:,None]-other)+ee[:,None]-pe,V)@self.w
        result=dict(eta=float(eta),lam=lam,performance_variance=performance_variance,V=V,
                    k=cost.k,E=cost.E if cost.kind=="saturating" else 0.,c=cost.c,
                    kind=cost.kind,W_F=WF,W_N=WN,J=Q-WF,S=Q-WN,
                    gap=WN-WF,effect_percent=100*(WN-WF)/WF,
                    collocation_residual=err,offgrid_residual=float(np.max(np.abs(off))),
                    iterations=it+1)
        return result,e

    def channels(self,cost,eta,lam=1.,performance_variance=1.):
        """Nested conditional integrals from the paper's three-channel formulas.

        q and u are evaluated analytically/by quadrature at every argument;
        the local channels do not use the solved equilibrium schedule.
        """
        V=performance_variance+2*(1-lam**2)
        r,w=self.z,self.w
        sig=math.sqrt(1-eta**2)
        def baseline(x):
            return cost.inverse(kernel(lam*(1-eta)*x,V+lam**2*(1-eta**2)))
        def response(x):
            x=np.asarray(x)
            other=eta*x[...,None]+sig*r
            qx=baseline(x)
            z=lam*(x[...,None]-other)
            gp=-z/V*kernel(z,V)
            return np.sum(gp*(qx[...,None]-baseline(other))*w,axis=-1)/cost.second(qx)
        other=eta*r[:,None]+sig*r[None,:]
        qr=baseline(r)
        qo=baseline(other)
        ur=response(r)
        uo=response(other)
        z=lam*(r[:,None]-other)
        gg=kernel(z,V)
        gp=-z/V*gg
        gpp=(z*z/V**2-1/V)*gg
        weight=w[:,None]*w[None,:]/cost.second(qr)[:,None]
        R=-2*float(np.sum(weight*gp*(ur[:,None]-uo)))
        D=-float(np.sum(weight*gpp*(qr[:,None]-qo)**2))
        C=float(np.sum(w*cost.third(qr)*ur**2/cost.second(qr)))
        return dict(response=R,density=D,cost_channel=C,S2=R+D+C)


def independent_effort(cost,eta,lam=1.,performance_variance=1.,ngrid=501,nq=160):
    """Independent uniform-grid cubic-spline fixed-point implementation."""
    V=performance_variance+2*(1-lam**2)
    z,w=normal_rule(nq)
    r=np.linspace(-10.,10.,ngrid)
    oth=eta*r[:,None]+math.sqrt(1-eta**2)*z[None,:]
    Z=lam*(r[:,None]-oth)
    e=cost.inverse(kernel(lam*(1-eta)*r,V+lam**2*(1-eta**2)))
    for it in range(900):
        spl=CubicSpline(r,e)
        er=cost.inverse(kernel(Z+e[:,None]-spl(np.clip(oth,-10,10)),V)@w)
        err=np.max(np.abs(er-e))
        if err<4e-14:
            e=er
            break
        e=.35*e+.65*er
    else:
        raise RuntimeError("independent iteration failed")
    return 2*float(w@CubicSpline(r,e)(np.clip(z,-10,10))),float(err)


def c3(k=.3,lam=1.,V=1.):
    g0=1/math.sqrt(2*math.pi*V)
    return math.sqrt(3)/(384*math.sqrt(V))*(5248*lam**4+1008*lam**2*V-27*V**2+576*V*(g0/k)**2)
