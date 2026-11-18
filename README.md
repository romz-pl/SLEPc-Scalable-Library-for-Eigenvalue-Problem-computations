# SLEPc: Scalable Library for Eigenvalue Problem Computations

## 1. What SLEPc is

[**SLEPc**](https://slepc.upv.es) is an open-source library for the solution of large-scale eigenvalue problems and related numerical linear algebra tasks on parallel computers. It is developed at the Universitat Politècnica de València (UPV), Spain, by a group led by José E. Román, Carmen Campos, Eloy Romero, Andrés Tomás, and Jose Roman's collaborators, among others. It is released under a permissive BSD-style license.

SLEPc is built on top of PETSc (the Portable, Extensible Toolkit for Scientific Computation). This relationship is central to its design. PETSc supplies the distributed vector and matrix data structures, the linear solvers, the preconditioners, the parallel communication layer based on MPI, and the options database. SLEPc adds the eigensolvers, the spectral transformation machinery, and the other solver families described below. A SLEPc program therefore looks much like a PETSc program: objects are created, configured from the command line or from code, solved, queried, and destroyed.

The name signals its target: the eigenvalue problems that arise in science and engineering are often too large to treat with dense, direct methods, so SLEPc focuses on iterative, projection-based methods that compute a few eigenpairs of very large sparse (or matrix-free) operators, using many processes or accelerators.

## 2. Problem classes covered

SLEPc is not limited to the standard eigenvalue problem. Its scope includes several families of problems, each with its own top-level solver object.

**Standard eigenvalue problem (EPS).** Find scalars $\lambda$ and nonzero vectors $x$ such that

$$
A x = \lambda x,
$$

where $A$ is a large, typically sparse, matrix that may be real or complex, symmetric (Hermitian) or non-symmetric, and may be given explicitly or only through its action on vectors (matrix-free).

**Generalized eigenvalue problem (EPS).** Find $\lambda$ and $x$ such that

$$
A x = \lambda B x,
$$

where $B$ is often symmetric positive definite, as in finite element mass matrices. Structural vibration and many other physical models lead to this form.

**Singular value decomposition (SVD).** Compute a few largest or smallest singular triplets $(\sigma, u, v)$ of a rectangular matrix $A$, satisfying $A v = \sigma u$ and $A^{H} u = \sigma v$. This is useful for low-rank approximation, principal component analysis on large data, and model reduction. SLEPc's SVD solver is a separate object (SVD) with its own set of algorithms.

**Polynomial eigenvalue problem (PEP).** Find $\lambda$ and $x$ such that

$$
P(\lambda) x = \left( \sum_{k=0}^{d} \lambda^{k} A_k \right) x = 0,
$$

for matrix coefficients $A_k$. The most important special case is the quadratic eigenvalue problem,

$$
\left( \lambda^{2} M + \lambda C + K \right) x = 0,
$$

which describes damped mechanical systems, with $M$, $C$, and $K$ the mass, damping, and stiffness matrices. Such problems are typically solved either by linearization into a larger generalized problem or by dedicated methods that work with the polynomial structure directly.

**Nonlinear eigenvalue problem (NEP).** Find $\lambda$ and $x$ such that

$$
T(\lambda) x = 0,
$$

where $T(\lambda)$ is a matrix-valued function of $\lambda$ that is not necessarily polynomial. Such problems arise in waveguide modelling, boundary element formulations, and delay differential equations. The NEP solver offers Newton-type and residual-inverse-iteration methods as well as nonlinear Arnoldi variants.

**Matrix functions (MFN).** Compute the action $f(A) v$ of a matrix function, such as the exponential, the square root, or the inverse square root, on a vector, without forming $f(A)$ explicitly. Time integration of linear systems $\dot{y} = A y$ is a classical application, since the solution is $y(t) = \exp(tA) y_0$.

## 3. Architecture and main objects

SLEPc is organized around a small set of user-facing solver objects and several internal components that they share.

The user-facing objects are EPS for standard and generalized eigenproblems, SVD for singular values, PEP for polynomial problems, NEP for nonlinear problems, and MFN for matrix functions. Each of these has a create, set-operators, configure, solve, query, and destroy cycle.

Internally, several components are reused across solvers. The spectral transformation object ST encapsulates how the problem is transformed before the eigensolver sees it, most notably shift-and-invert, Cayley, and polynomial transformations, along with the associated linear solves. The direct solver object DS manages small dense problems that arise inside the iteration, such as the projected matrix in a Krylov method. The basis vector object BV stores and orthogonalizes the Krylov or subspace basis, handling reorthogonalization strategies and distributed storage. The region object RG defines a region of the complex plane or an interval of interest, allowing filtering of eigenvalues. The function object FN represents the functions used by the NEP and MFN solvers.

This modular structure means that choosing a different eigensolver, a different spectral transformation, or a different orthogonalization strategy is typically a matter of changing one option rather than rewriting the program.

## 4. Eigensolver algorithms

SLEPc implements a broad catalogue of methods. The most important for the standard and generalized eigenproblems are the following.

**Krylov-Schur** is the default EPS method. It is a restarted Arnoldi method that retains a Schur-decomposed projection at restart, which gives robust and efficient restarts, particularly for selecting specific eigenvalues. It works for both symmetric and non-symmetric problems and is usually the best first choice.

**Arnoldi and Lanczos** are the classical Krylov projection methods. Lanczos exploits symmetry to use short recurrences, and both can be combined with restarts and with reorthogonalization to control loss of orthogonality.

**Subspace iteration** and the **power method** are simple, robust methods that are useful for teaching, for very well separated eigenvalues, or as building blocks. The Rayleigh-Ritz procedure is used to extract approximations.

**Davidson and Jacobi-Davidson** methods build subspaces using approximate inverses of the operator, which makes them effective with good preconditioners, especially for the interior eigenvalues of large problems.

**Generalized Davidson and LOBPCG** (locally optimal block preconditioned conjugate gradient) are block methods suited to symmetric problems where a good preconditioner is available. LOBPCG is particularly attractive for computing many eigenpairs simultaneously.

**Contour integral methods (CISS)** compute all eigenvalues inside a specified region of the complex plane by evaluating contour integrals of the resolvent. They are naturally parallel across quadrature points and are useful when the number of eigenvalues in a region is unknown but moderate.

**LAPACK-based dense solver** is available for small problems or as a reference, and is typically used to validate results on test cases.

For the SVD, SLEPc offers cross-product (applying the eigensolver to $A^{H}A$ or $A A^{H}$), cyclic matrix formulations, and Lanczos bidiagonalization variants including thick-restart Lanczos. For the polynomial problem, methods include linearization combined with Krylov-Schur, the quadratic Arnoldi method, the second-order Arnoldi and its variants, and Jacobi-Davidson style approaches. For the nonlinear problem, methods include successive linear problems, residual inverse iteration, nonlinear Arnoldi, and Newton-type approaches with interpolation.

The choice of algorithm depends on symmetry, the portion of the spectrum wanted, the cost of applying the operator, the availability of preconditioners, and the memory budget. SLEPc deliberately exposes these alternatives so that users can experiment.

## 5. Spectral transformations and interior eigenvalues

Krylov methods naturally converge first to eigenvalues at the periphery of the spectrum, such as those of largest magnitude. To compute eigenvalues near a target $\sigma$, one uses a spectral transformation. The most important is shift-and-invert, which applies the eigensolver to the operator

$$
(A - \sigma B)^{-1} B,
$$

whose largest eigenvalues correspond to the eigenvalues of the original problem closest to $\sigma$. Each application of this operator requires solving a linear system with $A - \sigma B$. In SLEPc this is handled through the PETSc KSP and PC infrastructure, so one can use direct factorizations (MUMPS, SuperLU_DIST, PARDISO, or MKL PARDISO, depending on the PETSc configuration) or preconditioned iterative solvers.

Shift-and-invert is extremely powerful but has a cost: the factorization of $A - \sigma B$ must be computed and held in memory, and its quality strongly affects performance. Other transformations include the Cayley transform, which is useful when eigenvalues are clustered in a particular region, and polynomial transformations, which can accelerate convergence without factorization. Users often need to tune the shift, choose an appropriate linear solver, and balance factorization cost against the number of iterations.

Selecting the wanted eigenvalues is controlled by the "which" setting. Options include largest and smallest magnitude, largest and smallest real part, largest and smallest imaginary part, eigenvalues of target magnitude, target real part, or target imaginary part, and a user-defined selection. A region can also be specified to restrict eigenvalues to an interval or a set in the complex plane.

## 6. Convergence control and diagnostics

Eigensolvers in SLEPc are controlled by a tolerance and a maximum number of iterations, and the user specifies how many eigenpairs are wanted and how large the working subspace may grow. Convergence is typically assessed through the residual norm

$$
\| A x - \lambda B x \| / \| x \|,
$$

or through a related estimate that can be computed cheaply during the iteration. After solving, the program queries the number of converged eigenpairs, retrieves eigenvalues and eigenvectors, and can inspect the error estimates. The tool also provides monitoring options that print the progress of convergence at each iteration, which is invaluable when diagnosing a solver that fails to converge.

Problems can be ill-conditioned, eigenvalues can be clustered, and non-normal operators can have sensitive spectra, so SLEPc provides options for balancing, for the orthogonalization procedure, for restart parameters, and for the size of the projected problem. Understanding these is a substantial part of using the library effectively.

## 7. Parallel computing and hardware support

SLEPc uses MPI for distributed-memory parallelism, inheriting the distribution of vectors and matrices from PETSc. A matrix is partitioned by rows across processes, and operations such as matrix-vector products, inner products, and norms are performed with local computation followed by collective communication. Basis vectors in the BV object are distributed in the same way, so the memory and computational work scale with the number of processes.

Several design choices support scalability. Orthogonalization of the basis can be performed with classical Gram-Schmidt with reorthogonalization, iterated classical Gram-Schmidt, modified Gram-Schmidt, or a tall-skinny QR approach, which reduces the number of global synchronization points. The CISS method is naturally parallel over quadrature points, which can be exploited through subcommunicators. Because the small projected problems are dense and replicated, their cost is kept modest by limiting the subspace size.

SLEPc also inherits GPU support through PETSc. Depending on the build, vectors and matrices can reside on NVIDIA or AMD GPUs, and PETSc's Kokkos backend provides a further route to heterogeneous architectures. The effectiveness of GPU offloading depends on the operator, the spectral transformation, and whether the linear solves can also run on the device.

Scalar types matter. SLEPc can be built for real or complex arithmetic. Complex builds are necessary for many physical problems, such as quantum mechanics and damped systems, and they are often used when eigenvalues are not real. The choice of scalar type is made at PETSc configuration time and affects the available algorithms and their behavior.

## 8. Typical application areas

SLEPc is used wherever large eigenproblems arise. In structural and mechanical engineering, modal analysis computes natural frequencies and mode shapes of finite element models, a generalized symmetric problem with a shift-and-invert approach being the standard workflow. In fluid dynamics, linear stability analysis of flows requires the eigenvalues of the linearized Navier-Stokes operator, typically non-symmetric and with interior eigenvalues of interest. In electronic structure and quantum chemistry, the Kohn-Sham and related Hamiltonians have large sparse or matrix-free representations, and many of their lowest eigenstates are needed. In acoustics and electromagnetics, the problems are often nonlinear or polynomial in frequency, making the NEP and PEP solvers relevant. In network and data science, the leading eigenvectors of graph Laplacians and adjacency matrices support spectral clustering, community detection, and ranking algorithms. In control theory and model reduction, singular value computations underpin balanced truncation and related methods.

## 9. Interfaces and ecosystem

SLEPc is primarily a C library, with Fortran interfaces for legacy and scientific codes. Python users typically employ slepc4py, which wraps the C API and integrates with petsc4py, allowing rapid prototyping while retaining the performance of the underlying compiled code. SLEPc is also used through the command line, since its options database is accessible through standard arguments such as the eigensolver type, the spectral transformation type, the number of eigenvalues wanted, the target, and the tolerance. This makes it easy to try different configurations without recompiling.

Because it is built on PETSc, SLEPc is widely embedded in higher-level finite element and simulation frameworks. Many codes built on PETSc can call SLEPc with little additional effort, and it is a common choice in scientific packages that need eigenvalues from a discretized partial differential equation.

## 10. A minimal C example

The following sketch illustrates the typical structure of a SLEPc program that computes the eigenvalues of a generalized problem $A x = \lambda B x$ closest to a target. Error checking is omitted for brevity; in real code every PETSc and SLEPc call should be wrapped in the standard error-handling macros.

```c
#include <slepceps.h>

int main(int argc, char **argv)
{
  Mat         A, B;
  EPS         eps;
  PetscInt    nconv;
  PetscScalar kr, ki;
  Vec         xr, xi;

  SlepcInitialize(&argc, &argv, NULL, NULL);

  /* Assemble or load the matrices A and B (distributed across MPI ranks) */
  /* ... MatCreate, MatSetValues, MatAssemblyBegin/End ... */

  EPSCreate(PETSC_COMM_WORLD, &eps);
  EPSSetOperators(eps, A, B);
  EPSSetProblemType(eps, EPS_GHEP);          /* generalized Hermitian */
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  EPSSetTarget(eps, 1.0);                    /* eigenvalues near 1.0 */
  EPSSetDimensions(eps, 6, PETSC_DEFAULT, PETSC_DEFAULT);
  EPSSetFromOptions(eps);                    /* allow command-line tuning */

  EPSSolve(eps);

  EPSGetConverged(eps, &nconv);
  for (PetscInt i = 0; i < nconv; i++) {
    EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
    /* use the eigenvalue kr and eigenvector xr */
  }

  EPSDestroy(&eps);
  SlepcFinalize();
  return 0;
}
```

Run with options such as `-eps_type krylovschur -st_type sinvert -st_ksp_type preonly -st_pc_type lu -eps_tol 1e-8` to select the algorithm, the spectral transformation, and the direct solver for the shift-and-invert step.

## 11. Installation and build

SLEPc is most conveniently built together with PETSc. The PETSc configure step can download and build SLEPc automatically, with an option such as `--download-slepc`, or SLEPc can be configured separately against an existing PETSc installation by setting the PETSC_DIR and PETSC_ARCH environment variables, then running the SLEPc configure script, followed by make and a verification target such as `make check`. Users should match the SLEPc version to the PETSc version, since the two are developed in tandem and releases are coordinated. Optional external packages, including MUMPS, SuperLU_DIST, Hypre, and ARPACK-based alternatives, are enabled through PETSc's configuration. The current release and documentation should be checked on the official SLEPc website, slepc.upv.es, since the version numbers and recommended build procedures change over time.

## 12. Practical guidance

Several points consistently determine success with SLEPc. The first is to identify the problem class correctly: a symmetric generalized problem, a quadratic polynomial, and a general nonlinear problem call for very different solvers. The second is to decide which eigenvalues are actually needed, since computing the largest eigenvalues is much cheaper than computing an interior cluster, which generally requires shift-and-invert and a factorization. The third is to choose the linear solver for the spectral transformation carefully, because for shift-and-invert it often dominates total cost; direct factorizations are robust but memory-hungry, while iterative solvers require good preconditioners. The fourth is to set the subspace dimensions and tolerances sensibly, as too small a subspace causes slow convergence and too large a subspace consumes memory and orthogonalization time. The fifth is to monitor convergence, using the built-in monitors, before trusting results, and to verify eigenpairs by computing residuals independently.

## 13. Relation to other eigensolver libraries

SLEPc occupies a particular niche among eigensolver libraries. ARPACK, in its ARPACK-NG form, is a widely used reverse-communication Arnoldi and Lanczos implementation that is serial or shared-memory and often serves as a back end for other tools; SLEPc offers a much broader and more parallel framework. Trilinos Anasazi provides a comparable block-eigensolver library within the Trilinos ecosystem. PRIMME specializes in symmetric and Hermitian problems with strong preconditioning and block methods. FEAST computes all eigenvalues within an interval using contour integration and is well suited to dense or sparse symmetric problems. For dense problems at scale, ScaLAPACK and ELPA are standard, and they do not target the large sparse regime that SLEPc addresses. The choice among these depends on the problem class, the existing software stack, and whether the application is already built on PETSc.

## 14. Summary

SLEPc is a comprehensive, parallel, and extensible library for eigenvalue and singular value computations on large sparse or matrix-free operators. Its strengths are its breadth of problem classes, from standard and generalized eigenproblems through SVD, polynomial and nonlinear problems, and matrix functions; its wide catalogue of algorithms and spectral transformations; its tight integration with PETSc, which brings scalable linear algebra, preconditioning, and hardware support; and its modular architecture, which lets users switch methods through options rather than code changes. Its main demands are that the user understand the numerical structure of the problem and the cost of the spectral transformation, and that the linear solver and subspace parameters be chosen with care. For researchers and engineers who need a few eigenpairs of very large operators on modern parallel machines, SLEPc is among the most complete and mature tools available.
