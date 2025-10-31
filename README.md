# Lightweight Numerical Optimization Tools (LNOT)
A list basic tools to solve non-constrainted optimization programs.

## Solving a Quadratic Program

We want to solve the following problem:
\f$\min_{x\in\mathbb{R}^n} \frac{1}{2}(x, Ax) + (b, x).\f$

LNOT's linear solver uses matrix-free Hessian operators \f$d\in\mathbb{R}^n \to Ad\in\mathbb{R}^n\f$. 
An hessian operator is a *functor* that define a `void operaor() (const Scalar* d, Scalar* Ad)`.
Writting such operator might be cumbersome. We thus provide the following class:
```cpp
LNOT::SymmetricDenseMatrixOp<double> Aop(A, N);
```
We build a *solver*, here the Conjugate Gradient:
```cpp
LNOT::ConjugateGradient<double> cg;
cg.setOutput(stdout); // the convergence history will be written in stdout
```
And solve the problem:
```cpp
cg.solve(Aop, b, N, x);
```
Each linear solver implements both a `solve` and a `solveWithGuess` methods that with the follwoing inputs:

- A, a Hessian operator
- b, a gradient vector
- x0 (for `solveWithGuess`) a first guess
- N the number of dimensions of the problem
- x the vector in which the solution will be written

Available solvers : 

- `ConjugateGradient`
- `LanczosSolver` *c.f.* Saad [1]
- `CoupledLanczosSolver` derived from Saad [1]

## Solving a Trust Region Subproblem (TRS)


We want to solve the following problem:
\f$\min_{x\in\mathbb{R}^n} \frac{1}{2}(x, Ax) + (b, x) \text{s.t.} \|x\|_2\leq\Delta.\f$

We build another type of solver based on Lanczos's iteration Gould *et. al.* [2]:
```cpp
LNOT::LanczosTRSSolver<double> lanczosTrs;
```
And solve the problem: 
```cpp
lanczosTrs.solve(Aop, b, N, delta, x);
```
Note that TRS solvers do not have a `solveWithGuess` method. The inputs are essentially the same than linear solvers but with an additional parameter delta, the size of the TR.

Available solvers : 

- `TruncatedConjugateGradient`[3]
- `LanczosTRSSolver` [2]
- `CoupledLanczosSolver` addapted from [2] 

## Solving a general unconstrained problem

We now want to solve the following problem \f$\min_{x\in\mathbb{R}^n} J(x)\f$ where \f$J\f$ is our cost function.
LNOT's general solvers relly on an *Oracle*, an object that for a given \f$x\f$ provides \f$J(x)\f$, \f$\nabla J(x)\f$ and \f$\nabla^2 J(x)\f$. Note that some oracle may not be able to provide \f$\nabla^2 J(x)\f$.

### What is an Oracle ? 

An oracle provides the following functions:

- `Size getNDims() const` the number of dimensions of the problem
- `void setCurrentPoint(const Scalar* x)` do whatever is necessary to compute \f$J(x)\f$, \f$\nabla J(x)\f$ and \f$\nabla^2 J(x)\f$. *e.g.* for an inverse problem the oracle might discretize the forward and adjoint equation with the model parameter \f$x\f$.
- `Scalar getValue() const` returns \f$J(x)\f$
- `void getGradient(Scalar* g) const` [optional] writes \f$\nabla J(x)\f$ in `g`
- `void getHessianProd(const Scalar* d, Scalar* Hd) const` [optional] writes \f$\nabla^2 J(x)d\f$ in `Hd`
- `void applyPrecond(const Scalar* d, Scalar* invBd) const` [optional] writes \f$B^{-1}d\f$ in `invBd` with \f$B\approx\nabla^2 J(x)\f$ in some norm.

If a large amout of pre-computation is required to evaluate \f$J\f$ and its derivatives, we strongly recommand to write a custom Oracle. 
On the other hand, for a functions such as Rosenbrok's function [4], we can simply define three functors:

```cpp
constexpr Size N = 10;

const auto func = [](const double* x)
{ 
	double y(0);
	for (Size i=0; i!=Size(N-1); ++i)
	{
		y += 100*(x[i+1] - x[i]*x[i])*(x[i+1] - x[i]*x[i]) + (1 - x[i])*(1 - x[i]);
	}
	return y;
};
const auto grad = [](const double* x, double* g)
{
	g[0] = -400*x[0]*(x[1] - x[0]*x[0]) + 2*(x[0] - 1);
	for (Size i=1; i!=Size(N-1); ++i)
	{
		g[i] = 200*(x[i] - x[i-1]*x[i-1]) - 400*x[i]*(x[i+1] - x[i]*x[i]) + 2*(x[i] - 1);
	}
	g[N-1] = 200*(x[N-1] - x[N-2]*x[N-2]);
};
const auto hessOp = [](const double* x, const double* d, double* Hd) 
{
	std::fill(Hd, Hd + N, 0);

	Hd[0] += (2 - 400*(x[1] - 3*x[0]*x[0]))*d[0];
	Hd[1] += -400*x[0]*d[0];
	
	for (Size i=1; i!=Size(N-1); ++i)
	{
		Hd[i-1] += -400*x[i-1]*d[i]; 
		Hd[i]   += (202 - 400*(x[i+1] - 3*x[i]*x[i]))*d[i];
		Hd[i+1] += -400*x[i]*d[i];
	}
	
	Hd[N-2] += -400*x[N-2]*d[N-1];
	Hd[N-1] += 200*d[N-1];
};
```

The solver will wrap theses three functors in an *OracleWrapper*.

### Solving the problem with a Newton/quasi-Newton method

A Newton/quasi-Newton method first computes a descent direction as \f$d_k = -B_k \nabla J(x_k)\f$, where \f$B_k\f$ is an approximation of the Hessian of \f$J\f$, and then computes the step lenght \f$\alpha_k\in\mathbb{R}\f$ with a *LineSearch* method.

Lets solve our problem with Newton's method:
```cpp
using GC = LNOT::ConjugateGradient<double>; // we will compute \f$d_k\f$ with the Conjgate Gradient 
using LS  = LNOT::BacktrackingLineSearch<double>; // we will compute \f$\alpha_k\f$ with a Backtracking Linesearch
LNOT::NewtonSolver<CG,LS> newtonSolver;
newtonSolver1.solve(func, grad, hessOp, N, x);
```
We can also solve the problem witn the L-BFGS method [5]:
```cpp
LNOT::LBFGS<LS> lbfgs1(5); // 5 is the number of steps used to estimate \f$B_k^{-1}\f$
lbfgs1.solve(func, grad, N, x); // here we do not need hessOp
```

Available Solvers:

- `LBFGS` [5]
- `NonLinearConjugateGradient` [8, 9]

Available LineSearches:

- `BacktrackingLineSearch` \f$\alpha\f$ is multiplied repeatedly by \f$\tau<1\f$ ultil \f$\alpha\f$ satisfies the Armijo-Goldstein condition
- `BisectionLineSearch`Use a bisection method to compute an \f$\alpha\f$ that satifies Wolfe's conditions

### Solving the problem with a Trust Region Method (TRM)

A TRM first fixes the step length, and searches a step \f$s_k\f$ within an *trust region* : \f$s_k = \arg\min_{s\in\mathbb{R}^n} (s, Hs) + (g,s) \text{s.t.} \|s\|\leq\Delta_k\f$. 
```cpp
using LanczosTrs = LNOT::LanczosTRSSolver<double>;
LNOT::NewtonTrustRegionSolver<LanczosTrs> trNewtonSolver;
trNewtonSolver.solve(func, grad, hessOp, N, x);
```
We can also solve the problem with both Symmetric Rank 1 (SR1) update of the Hessian matrix or Limited-memory SR1 (L-SR1) methods [6, 7]:
```
LNOT::SR1TrustRegionSolver<LanczosTrs> sr1TR2;
sr1TR2.solve(func, grad, N, x);
```

# Reference
[1] 4. Basic Iterative Methods. (2003). In Other Titles in Applied Mathematics. Iterative Methods for Sparse Linear Systems (pp. 103–128). doi:10.1137/1.9780898718003.ch4

[2] Gould, N. I. M., Lucidi, S., Roma, M., & Toint, P. L. (1999). Solving the Trust-Region Subproblem using the Lanczos Method. SIAM Journal on Optimization, 9(2), 504–525. doi:10.1137/S1052623497322735

[3] Steihaug, T. (1983). The Conjugate Gradient Method and Trust Regions in Large Scale Optimization. SIAM Journal on Numerical Analysis, 20(3), 626–637. doi:10.1137/0720042

[4] H. H. Rosenbrock, An Automatic Method for Finding the Greatest or Least Value of a Function, The Computer Journal, Volume 3, Issue 3, 1960, Pages 175–184, https://doi.org/10.1093/comjnl/3.3.175

[5] Liu, D.C., Nocedal, J. On the limited memory BFGS method for large scale optimization. Mathematical Programming 45, 503–528 (1989). https://doi.org/10.1007/BF01589116

[6] (2006). Quasi-Newton Methods. In: Numerical Optimization. Springer Series in Operations Research and Financial Engineering. Springer, New York, NY. https://doi.org/10.1007/978-0-387-40065-5_6

[7] Lu, X. (1996). A study of the limited memory SR1 method in practice. University of Colorado at Boulder.

[8] R. Fletcher, C. M. Reeves, Function minimization by conjugate gradients, The Computer Journal, Volume 7, Issue 2, 1964, Pages 149–154, https://doi.org/10.1093/comjnl/7.2.149

[9] Predrag S. Stanimirović, Branislav Ivanov, Haifeng Ma, Dijana Mosić. A survey of gradient methods for solving nonlinear optimization[J]. Electronic Research Archive, 2020, 28(4): 1573-1624. doi: 10.3934/era.2020115
