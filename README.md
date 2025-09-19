This repository contains a solver for the steady-state two-dimensional conduction problem (∂2T∂x2+∂2T∂y2=0∂x2∂2T+∂y2∂2T=0) on a unit square domain, discretized using the finite volume method. The mesh is uniform with Ni=Nj=128Ni=Nj=128 volumes, and lexicographic ordering by lines of constant ii is used. Dirichlet boundary conditions are applied: T=0T=0 at the left, right, and bottom walls, and T=1T=1 at the top wall.

Key features:

    Discretization leads to a system of equations in the form Ax=bAx=b.

    Linear system solved with conjugate gradient (CG), preconditioned CG with Jacobi, ILU, and SIP preconditioners.

    Residual 22-norm stopping criterion: 10−610−6.
