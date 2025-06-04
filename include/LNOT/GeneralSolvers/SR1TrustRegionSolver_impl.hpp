#ifndef LNOT_SR1_TRUST_REGION_SOLVER_IMPL_HPP
#define LNOT_SR1_TRUST_REGION_SOLVER_IMPL_HPP

#include <LNOT/GeneralSolvers/SR1TrustRegionSolver.hpp>
#include <LNOT/SymmetricDenseMatrixOp.hpp>
#include <LNOT/TRSSolvers/TruncatedConjugateGradient.hpp>
#include <LNOT/TRSSolvers/LanczosTRSSolver.hpp>
#include <LNOT/TRSSolvers/CoupledLanczosTRSSolver.hpp>

namespace LNOT
{

//// explicit template instanciations ////

extern template class SR1TrustRegionSolver< TruncatedConjugateGradient<float> >;
extern template class SR1TrustRegionSolver< LanczosTRSSolver<float> >;
extern template class SR1TrustRegionSolver< CoupledLanczosTRSSolver<float> >;

extern template class SR1TrustRegionSolver< TruncatedConjugateGradient<double> >;
extern template class SR1TrustRegionSolver< LanczosTRSSolver<double> >;
extern template class SR1TrustRegionSolver< CoupledLanczosTRSSolver<double> >;

//// method implementations ////

template<typename TRSSolver>
void SR1TrustRegionSolver<TRSSolver>::clearWorkSpace()
{
	if (m_gk     != nullptr) { delete[] m_gk;     m_gk     = nullptr; }
	if (m_gkp1   != nullptr) { delete[] m_gkp1;   m_gkp1   = nullptr; }
	if (m_uk     != nullptr) { delete[] m_uk;     m_uk     = nullptr; }
	if (m_sk     != nullptr) { delete[] m_sk;     m_sk     = nullptr; }
	if (m_Bk     != nullptr) { delete[] m_Bk;     m_Bk     = nullptr; }
	if (m_xTrial != nullptr) { delete[] m_xTrial; m_xTrial = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename TRSSolver> template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
void SR1TrustRegionSolver<TRSSolver>::solve_impl(Oracle& oracle, std::bool_constant<solveInPlace>, Scalar* x)
{
	const Scalar sr1DropTol = std::sqrt( std::numeric_limits<Scalar>::epsilon() );
	const Size size = oracle.getNDims();
	
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_gk     = new Scalar[Base::m_workCapacity];
		m_gkp1   = new Scalar[Base::m_workCapacity];
		m_uk     = new Scalar[Base::m_workCapacity];
		m_sk     = new Scalar[Base::m_workCapacity];
		m_Bk     = new Scalar[Base::m_workCapacity*Base::m_workCapacity];
		m_xTrial = new Scalar[Base::m_workCapacity];
	}
	if constexpr (not solveInPlace) { std::fill(x, x + size, 0); }	
	
	std::fill(m_Bk, m_Bk + size*size, 0);
	for (Size i=0; i!=size; ++i) { m_Bk[i + i*size] = 1; }
	
	SymmetricDenseMatrixOp<Scalar> BkOp(m_Bk, size);
	
	Base::m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk);
	
	Base::m_fx = oracle.getValue();
	Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	
	Scalar delta = std::pow(10.0, std::floor(std::log10(std::sqrt(Scalar(size)))));
	
	const Scalar tol2 = Base::m_tol*Base::m_tol*std::max(Scalar(1), Base::m_squaredNormGrad);
	
	if (Base::m_out != nullptr) { fmt::print(Base::m_out, "#SR1 Trust region method\n#Iteration f(x) delta residual tol\n"); }
	
	Base::m_info = Info::FAILURE;
	for (Base::m_nIt=0;Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{				
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}\n", Base::m_nIt, Base::m_fx, delta, Base::m_squaredNormGrad, tol2); }
		if (Base::m_squaredNormGrad < tol2) { Base::m_info = Info::SUCCESS; break; }
		
		m_trsSolver.solve(BkOp, m_gk, size, delta, m_sk); 
		
		Base::m_innerIts.push_back(m_trsSolver.getIterations());
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_xTrial[i] = x[i] + m_sk[i]; }
		
		oracle.setCurrentPoint(m_xTrial);
		oracle.getGradient(m_gkp1);
		BkOp(m_sk, m_uk);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_uk[i] = m_gkp1[i] - m_gk[i] - m_uk[i]; } // u_k = y_k - Bks_k
		
		const Scalar invRho = BasicLinalg::inner(m_sk, m_uk, size);
		
		if (std::abs(invRho) > sr1DropTol*BasicLinalg::squaredNorm(m_sk, size)*BasicLinalg::squaredNorm(m_uk, size))
		{
			BasicLinalg::symRk1Update(BkOp.getStorageOrder(), BkOp.getUplo(), 1. / invRho, m_uk, size, m_Bk);
		}
		
		const Scalar fxTrial   = oracle.getValue();
		const Scalar pred      = -m_trsSolver.getModelReduction();
		const Scalar ared      = Base::m_fx - fxTrial;
		const Scalar normS     = BasicLinalg::norm(m_sk, size);
		const Scalar tol_delta = std::max(delta*m_trsSolver.getTolTR(), std::numeric_limits<Scalar>::epsilon());
		
		const bool isStepFeasible       = std::isfinite(fxTrial);
		const bool isStepSuccessful     = ared > 0 and ared > pred*TRSBase::m_etaSuccessful;
		const bool isStepVerySuccessful = ared > 0 and ared > pred*TRSBase::m_etaVerySuccessful;
		const bool isStepAccepted       = isStepFeasible and ared > 0 and ared > pred*TRSBase::m_etaAccept;
		const bool normS_eq_delta       = std::abs(normS - delta) < tol_delta;
		
		if      (not (isStepSuccessful and isStepFeasible)) { delta *= TRSBase::m_gammaDecrease; }
		else if (isStepVerySuccessful and normS_eq_delta)   { delta *= TRSBase::m_gammaIncrease; }
		
		if (isStepAccepted) 
		{ 
			std::copy(m_xTrial, m_xTrial + size, x); 
			std::copy(m_gkp1,   m_gkp1   + size, m_gk); 
			Base::m_fx = fxTrial; 
			Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
		} 
		else
		{
			oracle.setCurrentPoint(x);
		}
		
		if (delta < std::numeric_limits<Scalar>::epsilon()) { Base::m_info = Info::BREAKDOWN; break; }    
	}
}

} // namespace LNOT

#endif // LNOT_SR1_TRUST_REGION_SOLVER_IMPL_HPP
