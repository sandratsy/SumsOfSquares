# SumsOfSquares
MATLAB and Mathematica code for reproducing the results found in the following paper:
```
Analysis of Optimization Algorithms via Sums-of-Squares
```
  
## This repository contains code to:
- Numerically verify the optimal contraction factors for all algorithms using YALMIP (*_YALMIP.m)
- Numerically verify the optimal contraction factors for one example algorithm using CVX
  - Gradient Descent with Exact Line Search (GDELS_CVX.m)
- Symbolically verify that the SOS certificates hold for all algorithms (*_verify.m)
- Symbolically verify that the matrix Q is positive semidefinite for two algorithms
  - Gradient Descent with Exact Line Search (GDELS_verify_Q_PSD.nb)
  - Proximal Gradient Method with Exact Line Search (PGMELS_verify_Q_PSD.nb)

## Requirements:
1. Scripts labelled 'CVX' require CVX and an SDP solver
2. Scripts labelled 'YALMIP' require YALMIP and an SDP solver
3. Scripts labelled 'verify' require MATLAB's symbolic computation toolbox

Each script can be executed independently.
