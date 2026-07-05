#ifndef LNOT_BFGS_IMPL_HPP
#define LNOT_BFGS_IMPL_HPP

#include <LNOT/GeneralSolvers/BFGS.hpp>
#include <LNOT/Containers/CircularBuffer.hpp>
#include <LNOT/LineSearches/BisectionLineSearch.hpp>
#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>
#include <LNOT/LineSearches/NoLineSearch.hpp>

namespace LNOT
{

//// explicit template instanciations ////

extern template class BFGS< BisectionLineSearch<float> >;
extern template class BFGS< BacktrackingLineSearch<float> >;
extern template class BFGS< NoLineSearch<float> >;

extern template class BFGS< BisectionLineSearch<double> >;
extern template class BFGS< BacktrackingLineSearch<double> >;
extern template class BFGS< NoLineSearch<double> >;

extern template class BFGS< BisectionLineSearch<long double> >;
extern template class BFGS< BacktrackingLineSearch<long double> >;
extern template class BFGS< NoLineSearch<long double> >;

//// method implementations ////

template<typename LineSearch, typename ConvergenceCriterion> template<CFirstOrderOracle Oracle, typename ABool> 
void BFGS<LineSearch, ConvergenceCriterion>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(isBool<ABool>)
{
	using Oracle_Size = typename Oracle::Size;
	
	const Oracle_Size size = oracle.getNDims();
	
	if (not m_gk or m_workCapacity < size)
	{
		m_workCapacity = size;

		m_gk    = std::make_unique<Scalar[]>(m_workCapacity);
		m_gkp1  = std::make_unique<Scalar[]>(m_workCapacity);
		m_sk    = std::make_unique<Scalar[]>(m_workCapacity);
		m_yk    = std::make_unique<Scalar[]>(m_workCapacity);
		m_uk    = std::make_unique<Scalar[]>(m_workCapacity);
		m_invBk = std::make_unique<Scalar[]>(m_workCapacity*m_workCapacity);
	}
	if (not solveInPlace) { std::fill(x, x + size, 0); }
	
	std::fill(m_invBk.get(), m_invBk.get() + size*size, 0);
	for (Size i=0; i!=size; ++i) { m_invBk[i + i*size] = 1; }
	
	m_innerIts.clear();
	
	constexpr FPComparator<Scalar> cmp;
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk.get());
	
	m_criterion.init(m_gk.get(), size);
	
	m_fx       = oracle.getValue();
	m_residual = m_criterion.getInitialResidual(); 
	
	const Scalar relTol = m_criterion.getRelTol(m_relTol);
	const Scalar absTol = m_criterion.getAbsTol(m_absTol);
	
	if (m_out != nullptr) { fmt::println(m_out, "#BFGS method\n#Iteration f(x) residual relative_tol absolute_tol"); }
	
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_fx, m_residual, relTol, absTol); std::fflush(m_out); }
		if (m_residual < relTol or m_residual < absTol) { m_info = Info::SUCCESS; break; }
		
		BasicLinalg::symMatrixVectorProd(StorageOrder::ROW_MAJOR, UpLo::LOWER, Scalar(-1), m_invBk.get(), m_gk.get(), size, BIC::fixed<bool,false>, m_sk.get());
		m_innerIts.push_back(1);
		
		////
		const Scalar step_len = m_lineSearch.solve(x, m_fx, m_gk.get(), m_sk.get(), oracle);
		
		if (not cmp.isDefPositive(step_len*BasicLinalg::norm(m_sk.get(), size))) { m_info = Info::BREAKDOWN; break; }
		
		BasicLinalg::scal(step_len, size, m_sk.get());
		BasicLinalg::axpy(Scalar(1), m_sk.get(), size, x);
		
		oracle.setCurrentPoint(x);
		oracle.getGradient(m_gkp1.get());
		
		#pragma omp simd
		for (Size j=0; j!=size; ++j) 
		{ 
			m_yk[j] = m_gkp1[j] - m_gk[j];
		}
		
		const Scalar sk_dot_yk = BasicLinalg::inner(m_sk.get(), m_yk.get(), size);
		
		if (sk_dot_yk > 0)
		{			
			BasicLinalg::symMatrixVectorProd(StorageOrder::ROW_MAJOR, UpLo::LOWER, Scalar(1), m_invBk.get(), m_yk.get(), size, BIC::fixed<bool,false>, m_uk.get());
			
			const Scalar yk_dot_uk = BasicLinalg::inner(m_yk.get(), m_uk.get(), size);
			
			BasicLinalg::symRk1Update(StorageOrder::ROW_MAJOR, UpLo::LOWER, (sk_dot_yk + yk_dot_uk) / (sk_dot_yk*sk_dot_yk), m_sk.get(), size, m_invBk.get());
			BasicLinalg::symRk2Update(StorageOrder::ROW_MAJOR, UpLo::LOWER, Scalar(-1) / sk_dot_yk, m_uk.get(), m_sk.get(), size, m_invBk.get());
		}
		
		////
		
		std::copy(m_gkp1.get(), m_gkp1.get() + size, m_gk.get());
		m_fx = oracle.getValue();
		m_residual = m_criterion.getResidual(m_gk.get(), size);
	}
}

} // namespace LNOT

#endif // LNOT_BFGS_IMPL_HPP
