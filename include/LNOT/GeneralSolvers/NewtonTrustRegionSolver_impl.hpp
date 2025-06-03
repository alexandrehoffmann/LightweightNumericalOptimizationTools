#ifndef LNOT_NEWTON_TRUST_REGION_SOLVER_IMPL_HPP
#define LNOT_NEWTON_TRUST_REGION_SOLVER_IMPL_HPP

#include <LNOT/GeneralSolvers/NewtonTrustRegionSolver.hpp>

namespace LNOT
{
	
template<typename TRSSolver>
void NewtonTrustRegionSolver<TRSSolver>::clearWorkSpace()
{
	if (m_gk     != nullptr) { delete[] m_gk;     m_gk     = nullptr; }
	if (m_sk     != nullptr) { delete[] m_sk;     m_sk     = nullptr; }
	if (m_xTrial != nullptr) { delete[] m_xTrial; m_xTrial = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename TRSSolver> template<SecondOrderOracle_concept Oracle, bool solveInPlace> 
void NewtonTrustRegionSolver<TRSSolver>::solve_impl(Oracle& oracle, std::bool_constant<solveInPlace>, Scalar* x)
{
	const Size size = oracle.getNDims();
	
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_gk     = new Scalar[Base::m_workCapacity];
		m_sk     = new Scalar[Base::m_workCapacity];
		m_xTrial = new Scalar[Base::m_workCapacity];
	}
	if constexpr (not solveInPlace) { std::fill(x, x + size, 0); }
	
	auto Hk = [&oracle] (const Scalar* __restrict__ d, Scalar* __restrict__ Hd) -> void {	oracle.getHessianProd(d, Hd); };
	
	Base::m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk);
	
	Base::m_fx = oracle.getValue();
	Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	
	Scalar delta = std::pow(10.0, std::floor(std::log10(std::sqrt(Scalar(size)))));
	
	const Scalar tol2 = Base::m_tol*Base::m_tol*std::max(Scalar(1), Base::m_squaredNormGrad);
	
	if (Base::m_out != nullptr) { fmt::print(Base::m_out, "#Newton Trust region method\n#Iteration f(x) delta residual tol\n"); }
	
	Base::m_info = Info::FAILURE;
	for (Base::m_nIt=0;Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{				
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}\n", Base::m_nIt, Base::m_fx, delta, Base::m_squaredNormGrad, tol2); }
		if (Base::m_squaredNormGrad < tol2) { Base::m_info = Info::SUCCESS; break; }
		
		m_trsSolver.solve(Hk, m_gk, size, delta, m_sk); 
		
		Base::m_innerIts.push_back(m_trsSolver.getIterations());
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_xTrial[i] = x[i] + m_sk[i]; }
		
		oracle.setCurrentPoint(m_xTrial);
		
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
			oracle.getGradient(m_gk);	
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

#endif // LNOT_NEWTON_TRUST_REGION_SOLVER_IMPL_HPP
