This repository contains a solver for the steady-state two-dimensional conduction problem on a unit square domain, discretized using the finite volume method. The mesh is uniform with Ni=Nj=128 volumes, and lexicographic orderingby lines of constant i is used. Dirichlet boundary conditions are applied: T=0at the left, right, and bottom walls, and T=1 at the top wall.

Key features:

   Discretization leads to a system of equations in the form Ax=b.

   Linear system solved with conjugate gradient (CG), preconditioned CG with Jacobi, ILU, and SIP preconditioners.

   Residual 2-norm stopping criterion: 10−6.
