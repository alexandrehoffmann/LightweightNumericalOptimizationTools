#ifndef LNOT_LBFGS_IMPL_HPP
#define LNOT_LBFGS_IMPL_HPP

#include <LNOT/GeneralSolvers/LBFGS.hpp>
#include <LNOT/Containers/CircularBuffer.hpp>
#include <LNOT/LineSearches/BisectionLineSearch.hpp>
#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>

#include <span>

namespace LNOT
{

//// explicit template instanciations ////

extern template class LBFGS< BisectionLineSearch<float> >;
extern template class LBFGS< BacktrackingLineSearch<float> >;

extern template class LBFGS< BisectionLineSearch<double> >;
extern template class LBFGS< BacktrackingLineSearch<double> >;

//// method implementations ////

template<typename LineSearch> 
void LBFGS<LineSearch>::clearWorkSpaceImpl()
{
	if (m_gk   != nullptr) { delete[] m_gk;   m_gk   = nullptr; }
	if (m_gkp1 != nullptr) { delete[] m_gkp1; m_gkp1 = nullptr; }
	if (m_dk   != nullptr) { delete[] m_dk;   m_dk   = nullptr; }
	if (m_S    != nullptr) { delete[] m_S;    m_S    = nullptr; }
	if (m_Y    != nullptr) { delete[] m_Y;    m_Y    = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename LineSearch> template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
void LBFGS<LineSearch>::solveImpl(Oracle& oracle, std::bool_constant<solveInPlace>, Scalar* x)
{
	using CircularBuffer_size = typename CircularBuffer<Scalar>::size_type;
	
	const Size size = oracle.getNDims();
	
	if (Base::m_workCapacity < size)
	{
		clearWorkSpaceImpl();
		Base::m_workCapacity = size;

		m_gk   = new Scalar[Base::m_workCapacity];
		m_gkp1 = new Scalar[Base::m_workCapacity];
		m_dk   = new Scalar[Base::m_workCapacity];
		m_S    = new Scalar[m_memory*Base::m_workCapacity];
		m_Y    = new Scalar[m_memory*Base::m_workCapacity];
	}
	if constexpr (not solveInPlace) { std::fill(x, x + size, 0); }
	
	CircularBuffer<Scalar> rho(m_memory);
	
	std::vector<Scalar> alpha(m_memory);
	
	Base::m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk);
	
	Base::m_fx = oracle.getValue();
	Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	
	const Scalar tol2 = Base::m_tol*Base::m_tol*std::max(Scalar(1), Base::m_squaredNormGrad);
	
	if (Base::m_out != nullptr) { fmt::print(Base::m_out, "#L-BFGS method\n#Iteration f(x) residual tol\n"); }
	
	Size curr_i = 0;
	Base::m_info = Info::FAILURE;
	for (Base::m_nIt=0;Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e} {:10.2e}\n", Base::m_nIt, Base::m_fx, Base::m_squaredNormGrad, tol2); }
		if (Base::m_squaredNormGrad < tol2) { Base::m_info = Info::SUCCESS; break; }
		
		//// two loop recursion
		#pragma omp simd
		for (Size j=0; j!=size; ++j) { m_dk[j] = -m_gk[j]; }
		rho.reverseForeach([&alpha, this, size](const CircularBuffer_size i, const Scalar& rho_i)
		{
			alpha[i] = rho_i*BasicLinalg::inner(m_S + i*size, m_dk, size);
			BasicLinalg::axpy(-alpha[i], m_Y + i*size, size, m_dk);
		});
		const Size offset = Size(size*rho.firstInsertedIndex());
		const Scalar gamma = rho.empty() ? Scalar(1) : BasicLinalg::inner(m_S + offset, m_Y + offset, size) / BasicLinalg::squaredNorm(m_Y + offset, size);
		BasicLinalg::scal(gamma, size, m_dk);
		rho.foreach([this, size, &alpha](const CircularBuffer_size i, const Scalar& rho_i)
		{
			const Scalar beta = rho_i*BasicLinalg::inner(m_Y + i*size, m_dk, size);
			BasicLinalg::axpy((alpha[i] - beta), m_S + i*size, size, m_dk);
		});
		Base::m_innerIts.push_back(1);
		////
		const Scalar step_len = m_lineSearch.solve(x, Base::m_fx, m_gk, m_dk, oracle);
		
		if (step_len < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::BREAKDOWN; break; }
		
		BasicLinalg::axpy(step_len, m_dk, size, x);
		oracle.setCurrentPoint(x);
		oracle.getGradient(m_gkp1);
		
		#pragma omp simd
		for (Size j=0; j!=size; ++j) 
		{ 
			m_S[j + curr_i*size] = step_len*m_dk[j];
			m_Y[j + curr_i*size] = m_gkp1[j] - m_gk[j];
		}
		
		const Scalar yk_dot_sk = BasicLinalg::inner(m_Y + curr_i*size, m_S + curr_i*size, size);
		
		if (yk_dot_sk > 0)
		{
			++curr_i; if (curr_i == m_memory) { curr_i = 0; }
			rho.push(1. / yk_dot_sk);	
		}
		
		////
		
		std::copy(m_gkp1, m_gkp1 + size, m_gk);
		Base::m_fx = oracle.getValue();
		Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	}
}
	
} // namespace LNOT

#endif // LNOT_LBFGS_IMPL_HPP
