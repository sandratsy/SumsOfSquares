# SumsOfSquares
MATLAB and Mathematica code for reproducing the results found in the following paper, as well as the supplementary material the paper references to
```
A Unified Framework for the Convergence Analysis of Optimization Algorithms via Sum-of-Squares
```
  
## This repository contains code to:
- Numerically verify the optimal contraction factors for three algorithms using CVX
  - Gradient Descent with Exact Line Search (GDELS_CVX.m)
  - Proximal Gradient Method with constant step size (PGM_const_CVX.m)
  - Proximal Gradient Method with Exact Line Search (PGMELS_CVX.m)
- Numerically verify the optimal contraction factors for the same three algorithms using YALMIP (*_YALMIP.m)
- Symbolically verify that the set of linear equalities hold for all three algorithms (*_verify_lin_eq.m)
- Symbolically verify that the matrix Q is positive semidefinite for all three algorithms (*_verify_Q_PSD.nb)

## Requirements:
1. Scripts labelled 'CVX' require CVX and an SDP solver
2. Scripts labelled 'YALMIP' require YALMIP and an SDP solver
3. Scripts labelled 'verify_lin_eq' require MATLAB's symbolic computation toolbox

Each script can be executed independently.

## Additional comments:
- While the optimal contraction factors given by CVX and YALMIP are the same, the SOS multipliers need not be, as the SOS decomposition may not be unique. The results given in the paper correspond to the SOS multipliers given by CVX, not YALMIP.
- The optimal Q's for all three algorithms are sparse, with only a submatrix within each Q being non-empty. The Mathematica code directly checks for the PSD-ness of the submatrix, instead of the full Q.
