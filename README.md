# finite-volume-Euler

Implemented Godunov's method for 1D Euler equations.

Exact Riemann solvers.


Inputs :

Main :
1. Domain Length
2. Cell Count
3. Boundary Conditions (0 - Transparent)
4. Initial Conditions (Starting cell, ending cell, values)
5. Evaluation Time

Test :
1. Domain Length
2. Cell Count
3. Test Case (SOD / 123 / BLAST_LEFT / BLAST_RIGHT / BLAST_2)

Constants:
1. Gamma
2. CFL Coefficient
3. Riemann Solver (0 - Exact)
4. Pressure Approximation Method (TR / PV)
5. Pressure Approximation Tolerance
6. Plotting factor for cells
7. Parallel Processing (Off / On)