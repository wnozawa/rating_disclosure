# Numerical replication

This directory reproduces the numerical figures, policy and channel tables,
root comparisons, and outcome calculations in the paper and supplement.
The analytical proofs are in `rating_disclosure.tex` at the package root.

Run from the package root:

~~~sh
python3 -m pip install -r numerics/requirements.txt
OPENBLAS_NUM_THREADS=1 python3 numerics/reproduce.py
sh build.sh
~~~

To regenerate displays from the included CSVs without recomputing equilibria:

~~~sh
OPENBLAS_NUM_THREADS=1 python3 numerics/reproduce.py --render-only
~~~

The full run is deterministic and normally takes several minutes. It overwrites
its CSVs and diagnostics. The optional `--resume` flag reuses completed stages
and is intended only for an interrupted run with unchanged code and settings.
The separate 18-point quadratic scan has its own workflow in `validation/`.

## Outputs and the integrated manuscript

Figure PDFs are written directly to `figures/`. Generated LaTeX table and
numeric-macro audit fragments are written to `numerics/generated/`, created
on demand. These fragments are not compilation inputs: all manuscript content,
including edited captions and tables, is already embedded in
`rating_disclosure.tex`. After changing numerical parameters, compare or merge
the generated values with the embedded tables and macros before compiling.
The distributed PNGs are 300-dpi previews of the included PDF figures; the
numerical renderer regenerates the PDF figures.

## Code and settings

- `solver.py`: effort costs, Chebyshev collocation, conditional Hermite
  quadrature, analytic Newton Jacobian, independent spline fixed-point solver,
  local channels, and the analytical third-order switching coefficient.
- `reproduce.py`: 924 policy evaluations in 14 configurations, 69 channel
  evaluations, numerical switching roots (including a finer grid
  $\bar e=5,5.25,\dots,15$ for the third-order error display), matched-scale gaps, cost-level
  comparisons, five outcome comparisons, resolution checks, and
  independent-solver comparisons.
- `render.py`: displays generated from the CSVs; only its output paths were
  adapted for this integrated-source package.

Production settings are 161 collocation nodes on [-10, 10] and 160 normal
quadrature nodes. Resolution checks use 193 nodes, 220 quadrature nodes, and
[-12, 12]. The independent comparison uses 801 spline nodes and 220 quadrature
nodes. Supplement S5 gives the method and residual definitions.

Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, and Matplotlib 3.10.8 were used
for the November 6 replication. The November 10 revision changed the display
labels in `render.py` to the paper's notation and added the finer-grid
switching roots (`data/switching_dense.csv`) used in the absolute-error
figure. A full rerun with Python 3.12.3 and the same pinned libraries
reproduced every previously shipped CSV in `data/` byte for byte.
`requirements.txt` pins library versions.
`data/validation_summary.json` records the checks and current source hashes.
`data/display_manifest.json` records matching display hashes and the rendering
environment. The figure PDFs embed font subsets, so their bytes can also depend
on the installed fontTools version; the rendered images do not. These
floating-point checks measure numerical sensitivity rather than certify
continuum equilibrium errors or a global policy ranking.

## Code names and paper notation

File names, CSV columns, and function names keep the identifiers used when the
data were generated. Displays use the paper's notation.

| Code / CSV | Paper |
| --- | --- |
| `E` (saturating-cost effort bound; files `policy_E*.csv`) | $\bar e$ |
| `c` (cubic coefficient; files `policy_c*.csv`) | $\kappa$ in $c(e)=ke^2/2+\kappa e^3/6$ |
| `V` | $\sigma^2$ |
| `performance_variance` | $\tau^2$ |
| `lam` | $\lambda$ |
| `W_F`, `W_N` | $W_{\mathrm F}$, $W_{\mathrm N}$ |
| `gap` | $W_{\mathrm N}-W_{\mathrm F}=J-S$ |
| `eta_star`, `z_star` | $\widehat\eta^\dagger$, $z^\dagger$ |
| `diff3` in `switching_dense.csv` | third-order approximation minus $\widehat\eta^\dagger$ |
| `response`, `density`, `cost_channel` | $S^{\mathrm R}$, $S^{\mathrm D}$, $S^{\mathrm C}$ |
| `S2`, `remainder` | $S_{\mathrm{loc}}$, $S-S_{\mathrm{loc}}$ |
| `B` in `reproduce.py` (scales `scaled_gap` in `matched_gap.csv`) | $B_0$ |
| `solver.c3()` | $a_2(k,\lambda,\sigma)$ |
| `p_H_F`, `p_H_N`, `U_F`, `U_N` | $p_{H,\mathrm F}$, $p_{H,\mathrm N}$, $U_{\mathrm F}$, $U_{\mathrm N}$ |

## Figure-to-data mapping

| PDF figure in figures/ | CSVs in numerics/data/ |
| --- | --- |
| S1_FIG1_relative_policy_gap.pdf | policy_E3.csv, policy_E6.csv, policy_E10.csv, policy_quadratic.csv |
| S2_FIG2_jensen_vs_feedback.pdf | policy_E6.csv, roots.csv |
| S4_FIG4_normalized_switch.pdf | switching.csv |
| S3_APP_FIG_absolute_errors.pdf | switching.csv, switching_dense.csv |
| S5_APP_FIG_scaled_policy_gap_recomputed.pdf | matched_gap.csv |
| N4_SUPP_policy_zoom.pdf | policy_E6.csv, policy_E10.csv, policy_quadratic.csv |
| S4_APP_FIG_cubic_policy_effects.pdf | cubic and quadratic policy CSVs, roots.csv |
| N4_SUPP_noisy_ratings.pdf | policy_E6.csv, policy_noisy.csv, roots.csv |
| S2_FIG3_signed_channels.pdf | channels.csv |
| N4_SUPP_endpoint_components.pdf | policy_E6.csv, policy_quadratic.csv |

The supplementary cost-level discussion uses `cost_levels.csv`. Outcome
comparisons use `outcomes.csv`. Validation tables use `validation.csv`,
`root_validation.csv`, and `validation_summary.json`. LaTeX assigns numbering
independently of the stable figure filenames.
