#ifndef LNOT_LBFGS_IMPL_HPP
#define LNOT_LBFGS_IMPL_HPP

#include <LNOT/GeneralSolvers/LBFGS.hpp>
#include <LNOT/Containers/CircularBuffer.hpp>
#include <LNOT/LineSearches/BisectionLineSearch.hpp>
#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>
#include <LNOT/LineSearches/NoLineSearch.hpp>

#include <span>

namespace LNOT
{

//// explicit template instanciations ////

extern template class LBFGS< BisectionLineSearch<float> >;
extern template class LBFGS< BacktrackingLineSearch<float> >;
extern template class LBFGS< NoLineSearch<float> >;

extern template class LBFGS< BisectionLineSearch<double> >;
extern template class LBFGS< BacktrackingLineSearch<double> >;
extern template class LBFGS< NoLineSearch<double> >;

extern template class LBFGS< BisectionLineSearch<long double> >;
extern template class LBFGS< BacktrackingLineSearch<long double> >;
extern template class LBFGS< NoLineSearch<long double> >;

//// method implementations ////

template<typename LineSearch> template<CFirstOrderOracle Oracle, typename ABool> 
void LBFGS<LineSearch>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
{
	using Oracle_Size         = typename Oracle::Size;
	using CircularBuffer_Size = typename CircularBuffer<Scalar>::size_type;
	
	const Oracle_Size size = oracle.getNDims();
	
	if (not m_gk or m_workCapacity < size)
	{
		m_workCapacity = size;

		m_gk   = std::make_unique<Scalar[]>(m_workCapacity);
		m_gkp1 = std::make_unique<Scalar[]>(m_workCapacity);
		m_dk   = std::make_unique<Scalar[]>(m_workCapacity);
		m_S    = std::make_unique<Scalar[]>(m_memory*m_workCapacity);
		m_Y    = std::make_unique<Scalar[]>(m_memory*m_workCapacity);
	}
	if (not solveInPlace) { std::fill(x, x + size, 0); }
	
	CircularBuffer<Scalar> rho(m_memory);
	
	std::vector<Scalar> alpha(m_memory);
	
	m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk.get());
	
	m_fx = oracle.getValue();
	m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk.get(), size);
	
	const Scalar relTol2 = m_relTol*m_relTol*m_squaredNormGrad;
	const Scalar absTol2 = m_absTol*m_absTol;
	
	const FPComparator<Scalar> cmp;
	
	if (m_out != nullptr) { fmt::println(m_out, "#L-BFGS method\n#Iteration f(x) residual relative_tol absolute_tol"); }
	
	Size curr_i{};
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_fx, m_squaredNormGrad, relTol2, absTol2); std::fflush(m_out); }
		if (m_squaredNormGrad < relTol2 or m_squaredNormGrad < absTol2) { m_info = Info::SUCCESS; break; }
		
		//// two loop recursion
		#pragma omp simd
		for (Size j=0; j!=size; ++j) { m_dk[j] = -m_gk[j]; }
		rho.reverseForeach([&alpha, this, size](const CircularBuffer_Size i, const Scalar& rho_i)
		{
			alpha[i] = rho_i*BasicLinalg::inner(m_S.get() + i*size, m_dk.get(), size);
			BasicLinalg::axpy(-alpha[i], m_Y.get() + i*size, size, m_dk.get());
		});
		const Size offset = Size(size*rho.firstInsertedIndex());
		const Scalar gamma = rho.empty() ? Scalar(1) : BasicLinalg::inner(m_S.get() + offset, m_Y.get() + offset, size) / BasicLinalg::squaredNorm(m_Y.get() + offset, size);
		BasicLinalg::scal(gamma, size, m_dk.get());
		rho.foreach([this, size, &alpha](const CircularBuffer_Size i, const Scalar& rho_i)
		{
			const Scalar beta = rho_i*BasicLinalg::inner(m_Y.get() + i*size, m_dk.get(), size);
			BasicLinalg::axpy((alpha[i] - beta), m_S.get() + i*size, size, m_dk.get());
		});
		m_innerIts.push_back(1);
		////
		const Scalar step_len = m_lineSearch.solve(x, m_fx, m_gk.get(), m_dk.get(), oracle);
		
		if (not cmp.isDefPositive(step_len*BasicLinalg::norm(m_dk.get(), size))) { m_info = Info::BREAKDOWN; break; }
		
		BasicLinalg::axpy(step_len, m_dk.get(), size, x);
		oracle.setCurrentPoint(x);
		oracle.getGradient(m_gkp1.get());
		
		#pragma omp simd
		for (Size j=0; j!=size; ++j) 
		{ 
			m_S[j + curr_i*size] = step_len*m_dk[j];
			m_Y[j + curr_i*size] = m_gkp1[j] - m_gk[j];
		}
		
		const Scalar yk_dot_sk = BasicLinalg::inner(m_Y.get() + curr_i*size, m_S.get() + curr_i*size, size);
		
		if (yk_dot_sk > 0)
		{
			++curr_i; if (curr_i == m_memory) { curr_i = 0; }
			rho.push(1. / yk_dot_sk);	
		}
		
		////
		
		std::copy(m_gkp1.get(), m_gkp1.get() + size, m_gk.get());
		m_fx = oracle.getValue();
		m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk.get(), size);
	}
}
	
} // namespace LNOT

#endif // LNOT_LBFGS_IMPL_HPP
