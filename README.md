# A Centered Volterra--Fredholm Realization of the Riemann Xi-Function

This repository contains the manuscript and reproducibility materials for
*A Centered Volterra--Fredholm Realization of the Riemann Xi-Function* by
Jingwei Liu.

Repository: <https://github.com/AlexxxxxLiu/btheta-volterra-fredholm>

The paper proves a regularized Fredholm-determinant realization of the Riemann
Xi-function and develops its operator-theoretic closure. It explicitly does
not claim a proof of the Riemann hypothesis; spectral reality remains the
identified open problem.

## Files

- `main.tex`: complete manuscript.
- `references.bib`: bibliography.
- `verify_btheta.py`: independent numerical regression checks.
- `requirements.txt`: pinned Python dependencies for the checks.
- `btheta-volterra-fredholm-reproducibility.zip`: submission-ready ancillary
  archive containing the source and verification files.
- `main.pdf`: compiled manuscript.

## Build

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

## Reproducibility check

Create an isolated environment in this directory and run the checks:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python verify_btheta.py
```

The program requires no external data. It reports high-precision Fourier
normalization errors and a convergence table for midpoint Nystrom
discretizations. The checks are diagnostic only and are not used as evidence
for any theorem.

The archive `btheta-volterra-fredholm-reproducibility.zip` can be uploaded as
supplementary material without rebuilding it from the working tree.

## Scope

All numerical checks are diagnostic. The mathematical claims in the manuscript
are supported by the analytic arguments in the paper, not by computation.
