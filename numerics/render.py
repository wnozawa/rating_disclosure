"""Render publication figures and LaTeX tables exclusively from shipped CSVs."""
import csv
import json
import math
from io import BytesIO
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
DATA=ROOT/"data"
OUT=ROOT/"generated"
OUT.mkdir(exist_ok=True)
FIGURES=ROOT.parent/"figures"
FIGURES.mkdir(exist_ok=True)
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                     "pdf.fonttype":42,"ps.fonttype":42,"figure.dpi":140})
COLORS=["#17212b","#526375","#89939d","#101010","#6c7883","#a0a6ac"]
STYLES=["-","--","-.",":",(0,(5,2,1,2)),(0,(2,2))]


def read(name):
    with (DATA/name).open() as f:
        rows=list(csv.DictReader(f))
    for row in rows:
        for key,value in row.items():
            try:row[key]=float(value)
            except ValueError:pass
    return rows


def save(fig,name):
    # Finish the PDF in memory before replacing the public file. Readers
    # therefore never observe an unfinished stream during a rebuild.
    buf=BytesIO()
    fig.savefig(buf,format="pdf",bbox_inches="tight",metadata={"CreationDate":None,"ModDate":None})
    data=buf.getvalue()
    if not data.rstrip().endswith(b"%%EOF"):
        raise RuntimeError("Incomplete PDF output: "+name)
    temporary=FIGURES/(name+".tmp")
    temporary.write_bytes(data)
    temporary.replace(FIGURES/name)
    plt.close(fig)


def curve(ax,rows,y,**kw):
    ax.plot([r["eta"] for r in rows],[r[y] for r in rows],**kw)


def decorate(ax,ylabel):
    ax.axhline(0,color="#555555",lw=.7)
    ax.grid(alpha=.15)
    ax.set(xlim=(0,1),xlabel=r"Matched-rating correlation $\eta$",ylabel=ylabel)


def sci(x,digits=2):
    if x==0:return "0"
    exp=int(math.floor(math.log10(abs(x))))
    return rf"{x/10**exp:.{digits}f}\times10^{{{exp}}}"


def table(filename,caption,label,columns,header,rows,aliases=()):
    tex=[r"\begin{table}[!htbp]",r"\centering\small",
         r"\caption{"+caption+"}",r"\label{"+label+"}"]
    tex += [r"\label{"+a+"}" for a in aliases]
    tex += [r"\setlength{\tabcolsep}{5pt}",r"\begin{tabular}{@{}"+columns+r"@{}}",
            r"\toprule"," & ".join(header)+r" \\ \midrule"]
    # Typeset leading minus signs as math minus, not hyphens.
    rows=[[r"$-$"+c[1:] if (len(c)>1 and c[0]=="-" and c[1].isdigit()) else c for c in row] for row in rows]
    tex += [" & ".join(row)+r" \\" for row in rows]
    tex += [r"\bottomrule",r"\end{tabular}",r"\end{table}",""]
    (OUT/filename).write_text("\n".join(tex))


def render_all():
    policies={p.stem.removeprefix("policy_"):read(p.name) for p in DATA.glob("policy_*.csv")}
    roots={r["case"]:r["eta_star"] for r in read("roots.csv")}
    switches=read("switching.csv")
    channels=read("channels.csv")
    summary=json.loads((DATA/"validation_summary.json").read_text())
    # Main figure 1.
    fig,ax=plt.subplots(figsize=(6.7,3.6),layout="constrained")
    for i,(key,label) in enumerate([("E3",r"Saturating, $\bar{e}=3$"),("E6",r"Saturating, $\bar{e}=6$"),
                                   ("E10",r"Saturating, $\bar{e}=10$"),("quadratic","Quadratic")]):
        curve(ax,policies[key],"effect_percent",color=COLORS[i],ls=STYLES[i],lw=1.7,label=label)
    decorate(ax,r"$100(W_{\mathrm{N}}-W_{\mathrm{F}})/W_{\mathrm{F}}$ (%)")
    ax.legend(frameon=False,ncol=2,loc="upper center",bbox_to_anchor=(.5,1.22))
    save(fig,"S1_FIG1_relative_policy_gap.pdf")
    # Main figure 2: the crossing itself is included.
    fig,ax=plt.subplots(figsize=(6.7,3.65),layout="constrained")
    data=policies["E6"]
    for key,label,style,color in [("gap",r"$100(J-S)/W_{\mathrm{F}}$","-",COLORS[0]),
                                  ("J",r"$100J/W_{\mathrm{F}}$","--",COLORS[1]),
                                  ("S",r"$100S/W_{\mathrm{F}}$",":",COLORS[2])]:
        ax.plot([r["eta"] for r in data],[100*r[key]/r["W_F"] for r in data],
                label=label,ls=style,color=color,lw=1.8)
    ax.axvline(roots["E6"],ls=":",color="#999999",lw=.8)
    ax.scatter([roots["E6"]],[0],color=COLORS[0],s=20,zorder=5)
    ax.annotate(r"$\widehat{\eta}^{\dagger}\simeq"+f'{roots["E6"]:.3f}'+r"$",
                (roots["E6"],0),xytext=(12,-27),textcoords="offset points",fontsize=9)
    decorate(ax,"Percentage of full-disclosure effort")
    ax.legend(frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(.5,1.17))
    save(fig,"S2_FIG2_jensen_vs_feedback.pdf")
    # Supplemental zoom of the policy comparison.
    fig,ax=plt.subplots(figsize=(6.7,3.4),layout="constrained")
    for i,(key,label) in enumerate([("E6",r"Saturating, $\bar{e}=6$"),("E10",r"Saturating, $\bar{e}=10$"),
                                   ("quadratic","Quadratic")]):
        curve(ax,policies[key],"effect_percent",ls=STYLES[i],color=COLORS[i],label=label)
    decorate(ax,r"$100(W_{\mathrm{N}}-W_{\mathrm{F}})/W_{\mathrm{F}}$ (%)")
    ax.set(xlim=(.72,1),ylim=(-.7,.55))
    ax.legend(frameon=False)
    save(fig,"N4_SUPP_policy_zoom.pdf")
    # Supplemental endpoint comparison formerly in main figure 2.
    fig,axs=plt.subplots(1,2,figsize=(9,3.1),layout="constrained")
    for ax,key,title in zip(axs,["E6","quadratic"],[r"Saturating, $\bar{e}=6$","Quadratic"]):
        data=[r for r in policies[key] if .8<=r["eta"]<=.99]
        for q,label,style,color in [("gap",r"$100(J-S)/W_{\mathrm{F}}$","-",COLORS[0]),
                                     ("J",r"$100J/W_{\mathrm{F}}$","--",COLORS[1]),
                                     ("S",r"$100S/W_{\mathrm{F}}$",":",COLORS[2])]:
            ax.plot([r["eta"] for r in data],[100*r[q]/r["W_F"] for r in data],label=label,ls=style,color=color)
        decorate(ax,"Percentage of full-disclosure effort")
        ax.set(xlim=(.8,.99),title=title)
    axs[0].legend(frameon=False,fontsize=8)
    save(fig,"N4_SUPP_endpoint_components.pdf")
    # Supplemental signed local channels.
    fig,axs=plt.subplots(1,2,figsize=(9,3.2),layout="constrained")
    terms=[("response",r"Response $S^{\mathrm{R}}/S$"),("density",r"Density $S^{\mathrm{D}}/S$"),
           ("cost_channel",r"Cost $S^{\mathrm{C}}/S$"),("remainder",r"Remainder $(S-S_{\mathrm{loc}})/S$")]
    for i,(term,label) in enumerate(terms):
        for ax in axs:
            ax.plot([r["eta"] for r in channels],[100*r[term]/r["S"] for r in channels],
                    color=COLORS[i],ls=STYLES[i],label=label,lw=1.6)
    for ax in axs:decorate(ax,"Share of strategic loss $S$ (%)")
    axs[0].set(title="Signed contributions",xlim=(0,.99))
    axs[1].set(title="Smaller contributions near perfect matching",xlim=(.8,.99),ylim=(-.1,3.1))
    axs[0].legend(frameon=False,fontsize=8,loc="center left")
    save(fig,"S2_FIG3_signed_channels.pdf")
    # Supplemental cubic family.
    fig,ax=plt.subplots(figsize=(7,3.8),layout="constrained")
    for i,(key,label) in enumerate([("quadratic","0"),("c0.01","0.01"),("c0.03","0.03"),
                                   ("c0.1","0.1"),("c0.3","0.3"),("c1","1")]):
        curve(ax,policies[key],"effect_percent",ls=STYLES[i],color=COLORS[i],label=r"$\kappa="+label+"$",lw=1.6)
    decorate(ax,r"$100(W_{\mathrm{N}}-W_{\mathrm{F}})/W_{\mathrm{F}}$ (%)")
    ax.legend(frameon=False,ncol=3)
    save(fig,"S4_APP_FIG_cubic_policy_effects.pdf")
    # Additional noisy-rating example.
    fig,ax=plt.subplots(figsize=(6.7,3.4),layout="constrained")
    for i,(key,label) in enumerate([("E6",r"$\lambda=1$, $\sigma^2=1$"),("noisy",r"$\lambda=0.8$, $\sigma^2=1.72$")]):
        curve(ax,policies[key],"effect_percent",ls=STYLES[i],color=COLORS[i],label=label,lw=1.7)
    decorate(ax,r"$100(W_{\mathrm{N}}-W_{\mathrm{F}})/W_{\mathrm{F}}$ (%)")
    ax.legend(frameon=False)
    save(fig,"N4_SUPP_noisy_ratings.pdf")
    # Normalized root refinement.
    fig,ax=plt.subplots(figsize=(6.7,3.6),layout="constrained")
    selected=[r for r in switches if r["E"]>=20]
    dense=np.linspace(20,100,300)
    from solver import c3  # c3() is the coefficient a_2(k,lambda,sigma) in the paper
    ax.plot(dense,np.zeros_like(dense)+math.sqrt(3),":",color=COLORS[2],label="Leading limit")
    ax.plot(dense,math.sqrt(3)+1.375/dense,"--",color=COLORS[1],label="Second order")
    ax.plot(dense,math.sqrt(3)+1.375/dense+c3()/dense**2,"-",color=COLORS[0],label="Third order")
    ax.scatter([r["E"] for r in selected],[r["z_star"] for r in selected],facecolors="white",
               edgecolors=COLORS[0],label="Numerical equilibrium",zorder=5)
    ax.set(xlabel=r"Effort bound $\bar{e}$",ylabel=r"$z^{\dagger}(\bar{e})=\bar{e}(1-\widehat{\eta}^{\dagger})$")
    ax.grid(alpha=.15);ax.legend(frameon=False,ncol=2,fontsize=8)
    save(fig,"S4_FIG4_normalized_switch.pdf")
    # Absolute errors on log axes. Marker fill shows the sign: filled if the
    # approximation exceeds the numerical switch, open if it falls short.
    # The third-order curve uses the finer grid up to E=15, so its dip near
    # E=7 is drawn as the sign change it is rather than as a straight kink.
    dense=read("switching_dense.csv")
    from matplotlib.lines import Line2D
    fig,ax=plt.subplots(figsize=(6.7,3.5),layout="constrained")
    for i,(order,label) in enumerate([("first","First order"),("second","Second order"),
                                      ("third","Third order")]):
        EE=[r["E"] for r in switches]
        diff=[r[order]-r["eta_star"] for r in switches]
        if order=="third":
            line=sorted([(r["E"],abs(r["diff3"])) for r in dense]+
                        [(E,abs(d)) for E,d in zip(EE,diff) if E>dense[-1]["E"]])
        else:
            line=[(E,abs(d)) for E,d in zip(EE,diff)]
        ax.plot([a for a,_ in line],[b for _,b in line],ls=STYLES[i],color=COLORS[i],
                lw=1.4,label=label,zorder=2)
        for E,d in zip(EE,diff):
            ax.plot([E],[abs(d)],marker="o",ms=4.5,color=COLORS[i],
                    markerfacecolor=COLORS[i] if d>0 else "white",zorder=3)
    for a,b in zip(dense,dense[1:]):
        if a["diff3"]*b["diff3"]<0:
            Ez=a["E"]-a["diff3"]*(b["E"]-a["E"])/(b["diff3"]-a["diff3"])
            ax.axvline(Ez,ls=":",color="#999999",lw=.8,zorder=0)
            ax.annotate("third-order\nsign change",(Ez,2e-6),xytext=(6,0),
                        textcoords="offset points",fontsize=8,va="bottom",color="#555555")
    ax.set_xscale("log");ax.set_yscale("log")
    ax.set(xlabel=r"Effort bound $\bar{e}$",ylabel=r"Absolute error in $\widehat{\eta}^{\dagger}$")
    handles,labels=ax.get_legend_handles_labels()
    handles+=[Line2D([],[],ls="none",marker="o",ms=4.5,color="#444444"),
              Line2D([],[],ls="none",marker="o",ms=4.5,color="#444444",markerfacecolor="white")]
    labels+=[r"Approximation above $\widehat{\eta}^{\dagger}$",
             r"Approximation below $\widehat{\eta}^{\dagger}$"]
    ax.grid(alpha=.15);ax.legend(handles,labels,frameon=False,fontsize=8,ncol=2)
    save(fig,"S3_APP_FIG_absolute_errors.pdf")
    matched=read("matched_gap.csv")
    fig,ax=plt.subplots(figsize=(6.7,3.6),layout="constrained")
    for i,E in enumerate([20,40,60,100]):
        rr=[r for r in matched if r["E"]==E]
        ax.plot([r["z"] for r in rr],[r["scaled_gap"] for r in rr],
                ls=STYLES[i],color=COLORS[i],label=rf"$\bar{{e}}={E}$")
    zz=np.linspace(.25,3,300)
    ax.plot(zz,zz*zz*(3-zz*zz),"--",color="black",lw=1,label=r"Limit $z^2(3\sigma^2-z^2)$")
    ax.axhline(0,color="#888",lw=.6);ax.grid(alpha=.15)
    ax.set(xlabel=r"$z=\bar{e}(1-\eta)$",ylabel=r"$\bar{e}^{4}(W_{\mathrm{N}}-W_{\mathrm{F}})/B_0$")
    ax.legend(frameon=False,ncol=2,fontsize=8)
    save(fig,"S5_APP_FIG_scaled_policy_gap_recomputed.pdf")
    # Main summary table: representative saturation/cubic costs plus noisy ratings.
    labels={"E3":r"Saturating, $\bar e=3$","E6":r"Saturating, $\bar e=6$","E10":r"Saturating, $\bar e=10$",
            "quadratic":"Quadratic","c0.01":r"Cubic, $\kappa=0.01$","c0.1":r"Cubic, $\kappa=0.1$",
            "c0.3":r"Cubic, $\kappa=0.3$","noisy":r"Saturating, $\bar e=6$"}
    def row(key):
        a=policies[key][0]
        pattern=r"$\mathrm F\to\mathrm N$" if key in roots else (r"$\mathrm N$ on grid" if a["gap"]>0 else r"$\mathrm F$ on grid")
        return [labels.get(key,key),f'{a["lam"]:g}',f'{a["effect_percent"]:+.2f}'+r"\%",
                f'{roots[key]:.3f}' if key in roots else "--",pattern]
    table("N4_TABLE_policy_summary.tex",
          r"Representative numerical policy comparisons ($k=0.3$, $\tau^2=1$). "
          r"The effect is $100(W_{\mathrm N}-W_{\mathrm F})/W_{\mathrm F}$ at $\eta=0$; positive values favor concealment. "
          r"$\widehat\eta^\dagger$ denotes a numerically refined sign change. Grid patterns refer to 66 evaluated correlations through $\eta=0.995$.",
          "tab:n6-roots","lcrrl",["Cost",r"$\lambda$",r"Effect at $\eta=0$",r"$\widehat\eta^\dagger$","Grid pattern"],
          [row(key) for key in labels],aliases=["tab:cubic-benchmarks"])
    # Complete supplementary policy tables.
    for filename,label,keys,title in [
        ("S1_TABLE1_generated.tex","tab:numerics-saturating",
         ["E3","E4","E5","E6","E10","E20","E100","quadratic"],"Saturating and quadratic costs"),
        ("S4_TABLE_cubic_benchmarks.tex","tab:numerics-cubic",
         ["quadratic","c0.01","c0.03","c0.1","c0.3","c1"],"Cubic costs")]:
        rr=[]
        for key in keys:
            a=policies[key][0]
            name=("Quadratic" if key=="quadratic" else
                  (rf"$\bar e={a['E']:g}$" if key.startswith("E") else rf"$\kappa={a['c']:g}$"))
            rr.append([name,f'{a["effect_percent"]:+.2f}'+r"\%",
                       f'{roots[key]:.6f}' if key in roots else "--",
                       r"$\mathrm F\to\mathrm N$" if key in roots else (r"$\mathrm N$ on grid" if a["gap"]>0 else r"$\mathrm F$ on grid")])
        table(filename,title+r" ($k=0.3$, $\lambda=\tau^2=1$). The effect is "
              r"$100(W_{\mathrm N}-W_{\mathrm F})/W_{\mathrm F}$ at $\eta=0$. Patterns describe the evaluated grid.",
              label,"lrrl",["Cost","Effect",r"$\widehat\eta^\dagger$","Grid pattern"],rr)
    selected=[min(channels,key=lambda r:abs(r["eta"]-v)) for v in [0,.5,.85,.9,.95,.98,.99]]
    rr=[[f'{r["eta"]:.3f}',"$"+sci(r["S"],4)+"$"]+
        [f'{100*r[q]/r["S"]:.4f}' for q,_ in terms] for r in selected]
    table("theorem1_contribution_table.tex",
          r"Signed percentages of the strategic loss $S$ for $\bar e=6$. The first three contributions sum to $S_{\mathrm{loc}}$; "
          r"the last is its equilibrium approximation remainder. Baseline $k=0.3$, $\lambda=\tau^2=1$.",
          "tab:theorem1contributions","crrrrr",[r"$\eta$","$S$",r"$100S^{\mathrm R}/S$",
          r"$100S^{\mathrm D}/S$",r"$100S^{\mathrm C}/S$",r"$100(S-S_{\mathrm{loc}})/S$"],rr)
    table("S3_APP_TABLE_switching_generated.tex",
          r"Numerical switching correlations and asymptotic approximations ($k=0.3$, "
          r"$\lambda=\tau^2=1$). Coefficients are analytic, with no fitting to the roots.",
          "tab:thirdswitch","rrrrrr",[r"$\bar e$",r"Numerical $\widehat\eta^\dagger$","First order","Second order","Third order",r"$|\mathrm{err}_3|$"],
          [[f'{r["E"]:g}']+[f'{r[q]:.6f}' for q in ["eta_star","first","second","third"]]+
           ["$"+sci(r["err3"])+"$"] for r in switches])
    outcomes=read("outcomes.csv")
    table("N4_TABLE_outcomes.tex",
          r"Alternative evaluation criteria at selected equilibria ($k=0.3$, $\lambda=\tau^2=1$).",
          "tab:r4-welfare","lrrrrrr",["Cost",r"$\eta$",r"$W_{\mathrm N}-W_{\mathrm F}$",r"$p_{H,\mathrm F}$",r"$p_{H,\mathrm N}$",r"$U_{\mathrm F}$",r"$U_{\mathrm N}$"],
          [[("Quadratic" if r["kind"]=="quadratic" else r"$\bar e=6$"),f'{r["eta"]:g}',f'{r["gap"]:+.4f}']+
           [f'{r[q]:.4f}' for q in ["p_H_F","p_H_N","U_F","U_N"]] for r in outcomes])
    diag=[["Policy grid evaluations",f'{summary["policy_evaluations"]:,}'],
          ["Channel evaluations",str(summary["channel_evaluations"])],
          ["Maximum collocation FOC residual","$"+sci(summary["max_collocation_residual"])+"$"],
          ["Maximum sampled off-grid FOC residual","$"+sci(summary["max_offgrid_residual"])+"$"],
          ["Maximum policy-gap change under refinement","$"+sci(summary["max_resolution_gap_difference"])+"$"],
          ["Maximum independent-solver effort difference","$"+sci(summary["max_independent_effort_difference"])+"$"],
          ["Maximum switching-root change under refinement","$"+sci(summary["max_root_resolution_difference"])+"$"]]
    table("N4_TABLE_validation.tex",r"Diagnostics for the regenerated numerical displays. "
          r"Resolution and independent-solver comparisons use the stated validation cases; "
          r"off-grid residuals are sampled diagnostics.", "tab:s5-numerical-audit","lr",["Check","Result"],diag)
    # Macros ensure the narrative and its underlying data remain synchronized.
    sw={r["E"]:r for r in switches}
    values={
        "NfourEffectThree":f'{policies["E3"][0]["effect_percent"]:.2f}',
        "NfourEffectSix":f'{abs(policies["E6"][0]["effect_percent"]):.2f}',
        "NfourEffectQuad":f'{abs(policies["quadratic"][0]["effect_percent"]):.2f}',
        "NfourRootSix":f'{roots["E6"]:.3f}',"NfourRootTen":f'{roots["E10"]:.3f}',
        "NfourRootSixFull":f'{roots["E6"]:.6f}',
        "NfourCubicSmallRoot":f'{roots["c0.01"]:.3f}',
        "NfourCubicMiddleRoot":f'{roots["c0.1"]:.3f}',
        "NfourNoisyRoot":f'{roots["noisy"]:.3f}' if "noisy" in roots else "--",
        "NfourNoisyEffect":f'{policies["noisy"][0]["effect_percent"]:+.2f}',
        "NfourLeadingSix":f'{sw[6]["first"]:.3f}',
        "NfourThirdSix":f'{sw[6]["third"]:.3f}',
        "NfourThirdErrorTwenty":sci(sw[20]["err3"]),
        "NfourRefinementGap":sci(summary["max_resolution_gap_difference"]),
        "NfourIndependentGap":sci(summary["max_independent_effort_difference"]),
        "NfourAtwo":f'{c3():.4f}'}
    (OUT/"N4_NUMERIC_VALUES.tex").write_text("% Generated by numerics/render.py.\n"+
        "".join("\\newcommand{\\"+k+"}{"+v+"}\n" for k,v in values.items()))
    print("Rendered 10 figures, 7 tables, and synchronized numeric macros.",flush=True)


if __name__=="__main__":
    render_all()
