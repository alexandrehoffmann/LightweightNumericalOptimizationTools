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

//// method implementations ////


template<typename LinSolver, typename LineSearch> 
void NewtonSolver<LinSolver,LineSearch>::clearWorkSpace()
{
	if (m_gk != nullptr) { delete[] m_gk; m_gk = nullptr; }
	if (m_sk != nullptr) { delete[] m_sk; m_sk = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename LinSolver, typename LineSearch> template<SecondOrderOracle_concept Oracle, bool solveInPlace> 
void NewtonSolver<LinSolver,LineSearch>::solve_impl(Oracle& oracle, std::bool_constant<solveInPlace>, Scalar* x) requires (not Oracle::hasApplyPrecond)
{
	const Size size = oracle.getNDims();
	
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_gk = new Scalar[Base::m_workCapacity];
		m_sk = new Scalar[Base::m_workCapacity];
	}
	if constexpr (not solveInPlace) { std::fill(x, x + size, 0); }
	
	auto Hk = [&oracle] (const Scalar* d, Scalar* Hd) -> void {	oracle.getHessianProd(d, Hd); };
	
	Base::m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk);
	
	Base::m_fx = oracle.getValue();
	Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	
	const Scalar tol2 = Base::m_tol*Base::m_tol*std::max(Scalar(1), Base::m_squaredNormGrad);
	
	if (Base::m_out != nullptr) { fmt::print(Base::m_out, "#Newton method\n#Iteration f(x) residual tol\n"); }
	
	Base::m_info = Info::FAILURE;
	for (Base::m_nIt=0;Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{				
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e} {:10.2e}\n", Base::m_nIt, Base::m_fx, Base::m_squaredNormGrad, tol2); }
		if (Base::m_squaredNormGrad < tol2) { Base::m_info = Info::SUCCESS; break; }
		m_linsSolver.solve(Hk, m_gk, size, m_sk); 
		if (m_linsSolver.getInfo() == LinSolver::Info::NEGATIVE_CURVATURE) 
		{ 
			#pragma omp simd
			for (Size i=0; i!=size; ++i) { m_sk[i] = -m_gk[i]; }
		}
		Base::m_innerIts.push_back(m_linsSolver.getIterations());
		
		const Scalar alpha = m_lineSearch.solve(x, Base::m_fx, m_gk, m_sk, oracle);
		
		if (alpha < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::BREAKDOWN; break; }
		
		BasicLinalg::axpy(alpha, m_sk, size, x);
		
		oracle.setCurrentPoint(x);
		oracle.getGradient(m_gk);
		
		Base::m_fx = oracle.getValue();
		Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	}
}

template<typename LinSolver, typename LineSearch> template<SecondOrderOracle_concept Oracle, bool solveInPlace> 
void NewtonSolver<LinSolver,LineSearch>::solve_impl(Oracle& oracle, std::bool_constant<solveInPlace>, Scalar* 	x) requires (Oracle::hasApplyPrecond)
{
	const Size size = oracle.getNDims();
	
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_gk = new Scalar[Base::m_workCapacity];
		m_sk = new Scalar[Base::m_workCapacity];
	}
	if constexpr (not solveInPlace) { std::fill(x, x + size, 0); }
	
	auto Hk    = [&oracle] (const Scalar* d, Scalar* Hd)    -> void {	oracle.getHessianProd(d, Hd);  };
	auto invBk = [&oracle] (const Scalar* d, Scalar* invBd) -> void {	oracle.applyPrecond(d, invBd); };
	
	Base::m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk);
	
	Base::m_fx = oracle.getValue();
	Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	
	const Scalar tol2 = Base::m_tol*Base::m_tol*std::max(Scalar(1), Base::m_squaredNormGrad);
	
	if (Base::m_out != nullptr) { fmt::print(Base::m_out, "#Newton method\n#Iteration f(x) residual tol\n"); }
	
	Base::m_info = Info::FAILURE;
	for (Base::m_nIt=0;Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{				
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e} {:10.2e}\n", Base::m_nIt, Base::m_fx, Base::m_squaredNormGrad, tol2); }
		if (Base::m_squaredNormGrad < tol2) { Base::m_info = Info::SUCCESS; break; }
		m_linsSolver.solve(Hk, invBk, m_gk, size, m_sk); 
		if (m_linsSolver.getInfo() == LinSolver::Info::NEGATIVE_CURVATURE) 
		{ 
			#pragma omp simd
			for (Size i=0; i!=size; ++i) { m_sk[i] = -m_gk[i]; }
		}
		Base::m_innerIts.push_back(m_linsSolver.getIterations());
		
		const Scalar alpha = m_lineSearch.solve(x, Base::m_fx, m_gk, m_sk, oracle);
		
		if (alpha < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::BREAKDOWN; break; }
		
		BasicLinalg::axpy(alpha, m_sk, size, x);
		
		oracle.setCurrentPoint(x);
		oracle.getGradient(m_gk);
		
		Base::m_fx = oracle.getValue();
		Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	}
}

} // namespace LNOT

#endif // LNOT_NEWTON_SOLVER_IMPL_HPP
