#ifndef LNOT_NEWTON_SOLVER_IMPL_HPP
#define LNOT_NEWTON_SOLVER_IMPL_HPP

#include <LNOT/GeneralSolvers/NewtonSolver.hpp>
#include <LNOT/LineSearches/BisectionLineSearch.hpp>
#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>
#include <LNOT/LinearSolvers/ConjugateGradient.hpp>
#include <LNOT/LinearSolvers/LanczosSolver.hpp>
#include <LNOT/LinearSolvers/CoupledLanczosSolver.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

//// explicit template instanciations ////

// explicit instanciation for float
extern template class NewtonSolver<ConjugateGradient<float>,    BisectionLineSearch<float> >;
extern template class NewtonSolver<LanczosSolver<float>,        BisectionLineSearch<float> >;
extern template class NewtonSolver<CoupledLanczosSolver<float>, BisectionLineSearch<float> >;

extern template class NewtonSolver<ConjugateGradient<float>,    BacktrackingLineSearch<float> >;
extern template class NewtonSolver<LanczosSolver<float>,        BacktrackingLineSearch<float> >;
extern template class NewtonSolver<CoupledLanczosSolver<float>, BacktrackingLineSearch<float> >;

// explicit instanciation for double
extern template class NewtonSolver<ConjugateGradient<double>,    BisectionLineSearch<double> >;
extern template class NewtonSolver<LanczosSolver<double>,        BisectionLineSearch<double> >;
extern template class NewtonSolver<CoupledLanczosSolver<double>, BisectionLineSearch<double> >;

extern template class NewtonSolver<ConjugateGradient<double>,    BacktrackingLineSearch<double> >;
extern template class NewtonSolver<LanczosSolver<double>,        BacktrackingLineSearch<double> >;
extern template class NewtonSolver<CoupledLanczosSolver<double>, BacktrackingLineSearch<double> >;

// explicit instanciation for long double
extern template class NewtonSolver<ConjugateGradient<long double>,    BisectionLineSearch<long double> >;
extern template class NewtonSolver<LanczosSolver<long double>,        BisectionLineSearch<long double> >;
extern template class NewtonSolver<CoupledLanczosSolver<long double>, BisectionLineSearch<long double> >;

extern template class NewtonSolver<ConjugateGradient<long double>,    BacktrackingLineSearch<long double> >;
extern template class NewtonSolver<LanczosSolver<long double>,        BacktrackingLineSearch<long double> >;
extern template class NewtonSolver<CoupledLanczosSolver<long double>, BacktrackingLineSearch<long double> >;

//// method implementations ////


template<typename LinSolver, typename LineSearch> 
void NewtonSolver<LinSolver,LineSearch>::clearWorkSpaceImpl()
{
	if (m_gk != nullptr) { delete[] m_gk; m_gk = nullptr; }
	if (m_sk != nullptr) { delete[] m_sk; m_sk = nullptr; }
	m_workCapacity = 0;
}

template<typename LinSolver, typename LineSearch> template<SecondOrderOracle_concept Oracle, typename ABool> 
void NewtonSolver<LinSolver,LineSearch>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
{
	using Oracle_Size = typename Oracle::Size;
	
	const Oracle_Size size = oracle.getNDims();
	
	if (m_workCapacity < size)
	{
		clearWorkSpaceImpl();
		m_workCapacity = size;
		m_gk = new Scalar[m_workCapacity];
		m_sk = new Scalar[m_workCapacity];
	}
	if (not solveInPlace) { std::fill(x, x + size, 0); }
	
	const auto Hk    = [&oracle] (const Scalar* d, Scalar* Hd)    -> void { oracle.getHessianProd(d, Hd);  };
	const auto invBk = [&oracle] (const Scalar* d, Scalar* invBd) -> void { oracle.applyPrecond(d, invBd); };
	
	m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk);
	
	m_fx = oracle.getValue();
	m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	
	const Scalar relTol2 = m_relTol*m_relTol*m_squaredNormGrad;
	const Scalar absTol2 = m_absTol*m_absTol;
	
	const FPComparator<Scalar> cmp;
	
	if (m_out != nullptr) { fmt::print(m_out, "#Newton method\n#Iteration f(x) residual relative_tol absolute_tol\n"); }
	
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{				
		if (m_out) { fmt::print(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}\n", m_nIt, m_fx, m_squaredNormGrad, relTol2, absTol2); }
		if (m_squaredNormGrad < relTol2 or m_squaredNormGrad < absTol2) { m_info = Info::SUCCESS; break; }
		m_linsSolver.solve(Hk, invBk, m_gk, size, m_sk); 
		if (m_linsSolver.getInfo() == LinSolver::Info::NEGATIVE_CURVATURE) 
		{ 
			#pragma omp simd
			for (Size i=0; i!=size; ++i) { m_sk[i] = -m_gk[i]; }
		}
		m_innerIts.push_back(m_linsSolver.getIterations());
		
		const Scalar alpha = m_lineSearch.solve(x, m_fx, m_gk, m_sk, oracle);
		
		if (not cmp.isDefPositive(alpha)) { m_info = Info::BREAKDOWN; break; }
		
		BasicLinalg::axpy(alpha, m_sk, size, x);
		
		oracle.setCurrentPoint(x);
		oracle.getGradient(m_gk);
		
		m_fx = oracle.getValue();
		m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	}
}

} // namespace LNOT

#endif // LNOT_NEWTON_SOLVER_IMPL_HPP
