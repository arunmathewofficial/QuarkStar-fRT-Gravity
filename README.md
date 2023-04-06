# Quark Star in f(RT) Gravity

This is the source code for simulating a Quark star in f(R,T) = R + alpha R^2 + omega RT gravity, with an equation of state given by the bag model: p = k(epsilon - 4B).

The source code utilizes numerical methods to solve the Starobinsky model and f(R,T) = R + alpha R^2 + omega RT gravity. The code mainly performs two tasks:

1. Determines the accurate central value of the Ricci scalar for a given central density by solving the corresponding stellar model.

2. Uses this value as the initial guess for the central value of the Ricci scalar for the f(R,T) = R + alpha R^2 + omega RT model, the corresponding stelar model is solved for the same central density. The perturbation around the Starobinsky model is employed to solve the field equations of the f(R,T) = R + alpha R^2 + omega RT model, resulting in a new R_c for the f(R,T) model. However, this new value may not be consistent with the boundary conditions. To ensure convergence, the new value is taken as the initial value for the next simulation, and this process is repeated iteratively until a consistent value of R_c for the f(R,T) model is obtained for the given central density.

The results of the source code are published in: https://link.springer.com/article/10.1140%2Fepjc%2Fs10052-020-8130-4

Author: Arun Mathew
Institute: New Numerical Lab, Dept. of Physics, IIT Guwahati.
Date: 2019
