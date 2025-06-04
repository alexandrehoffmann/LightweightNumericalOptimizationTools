# LightweightNumericalOptimizationTools
A list basic tools to solve non-constrainted optimization programs

## Solving a Quadratic Program

We want to solve the following problem:
$\min_{x\in\mathbb{R}^n} \frac{1}{2}(x, Ax) + (b, x).$

LNOT's linear solver uses matrix-free Hessian operators $d\in\mathbb{R}^n \to Ad\in\mathbb{R}^n$. 
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
$\min_{x\in\mathbb{R}^n} \frac{1}{2}(x, Ax) + (b, x) \text{s.t.} \|x\|_2\leq\Delta.$

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

We now want to solve the following problem $\min_{x\in\mathbb{R}^n} J(x)$ where $J$ is our cost function.
LNOT's general solvers relly on an *Oracle*, an object that for a given $x$ provides $J(x)$, $\nabla J(x)$ and $\nabla^2 J(x)$. Note that some oracle may not be able to provide $\nabla^2 J(x)$.

### What is an Oracle

An oracle provides the following functions:

- `Size getNDims() const` the number of dimensions of the problem
- `void setCurrentPoint(const Scalar* x)` do whatever is necessary to compute $J(x)$, $\nabla J(x)$ and $\nabla^2 J(x)$. *e.g.* for an inverse problem the oracle might discretize the forward and adjoint equation with the model parameter $x$.
- `Scalar getValue() const` returns $J(x)$
- `void getGradient(Scalar* g) const` [optional] writes $\nabla J(x)$ in `g`
- `void getHessianProd(const Scalar* d, Scalar* Hd) const ` [optional] writes $\nabla J(x)d$ in `Hd`

If a large amout of pre-computation is required to evaluate $J$ and its derivatives, we strongly recommand to write a custom Oracle. 
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

A Newton/quasi-Newton method first computes a descent direction as $d_k = -B_k \nabla J(x_k)$, where $B_k$ is an approximation of the Hessian of $J$, and then computes the step lenght $\alpha_k\in\mathbb{R}$ with a *LineSearch* method.

Lets solve our problem with Newton's method:
```cpp
LNOT::ConjugateGradient<double> cg; // we will compute $d_k$ with the Conjgate Gradient 
LNOT::BacktrackingLineSearch<double> ls; // we will compute $\alpha_k$ with a Backtracking Linesearch
auto newtonSolver = LNOT::makeNewtonSolver(cg, ls);
newtonSolver1.solve(func, grad, hessOp, N, x);
```
We can also solve the problem witn the L-BFGS method [5]:
```cpp
auto lbfgs1 = LNOT::makeLBFGS(ls, 5); // 5 is the number of steps used to estimate $B_k^{-1}$
lbfgs1.solve(func, grad, N, x); // here we do not need hessOp
```

Available LineSearches:

- `BacktrackingLineSearch` $\alpha$ is multiplied repeatedly $\tau<1$ ultil $\alpha$ satisfies the Armijo-Goldstein condition
- `BisectionLineSearch`Use a bisection method to compute an $\alpha$ that satifies Wolfe's conditions

### Solving the problem with a Trust Region Method (TRM)

A TRM first fixes the step length, and searches a step $s_k$ within an *trust region* : $s_k = \arg\min_{s\in\mathbb{R}^n} (s, Hs) + (g,s) \text{s.t.} \|s\|\leq\Delta_k$. 

# Reference
[1] 4. Basic Iterative Methods. (2003). In Other Titles in Applied Mathematics. Iterative Methods for Sparse Linear Systems (pp. 103–128). doi:10.1137/1.9780898718003.ch4
[2] Gould, N. I. M., Lucidi, S., Roma, M., & Toint, P. L. (1999). Solving the Trust-Region Subproblem using the Lanczos Method. SIAM Journal on Optimization, 9(2), 504–525. doi:10.1137/S1052623497322735
[3] Steihaug, T. (1983). The Conjugate Gradient Method and Trust Regions in Large Scale Optimization. SIAM Journal on Numerical Analysis, 20(3), 626–637. doi:10.1137/0720042
[4] H. H. Rosenbrock, An Automatic Method for Finding the Greatest or Least Value of a Function, The Computer Journal, Volume 3, Issue 3, 1960, Pages 175–184, https://doi.org/10.1093/comjnl/3.3.175
[5] Liu, D.C., Nocedal, J. On the limited memory BFGS method for large scale optimization. Mathematical Programming 45, 503–528 (1989). https://doi.org/10.1007/BF01589116
