# Quadratic-cost scan

This directory reproduces the 18 configurations reported in Supplement S4.
Run from the package root:

~~~sh
OPENBLAS_NUM_THREADS=1 python3 validation/quadratic_global_scan.py
~~~

`quadratic_scan.csv` contains the reported data. The script writes
`quadratic_scan_regenerated.csv`. The grid is x in {0.5, 2, 10}, eta in
{0, 0.5, 0.9}, and k divided by the curvature bound in {1.02, 1.5},
with lambda = 1. Full-disclosure effort is evaluated analytically.
The equilibrium under concealment uses the spline fixed-point method
described in Supplement S4.

The November 6 replication reproduced all 18 reported gaps exactly in the
recorded environment. The smallest positive gap was approximately
0.0002217654. The iteration residual is a discrete diagnostic rather than a
certified continuum error bound.
