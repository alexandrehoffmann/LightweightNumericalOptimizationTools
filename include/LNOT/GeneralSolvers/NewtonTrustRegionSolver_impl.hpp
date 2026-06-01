#ifndef LNOT_NEWTON_TRUST_REGION_SOLVER_IMPL_HPP
#define LNOT_NEWTON_TRUST_REGION_SOLVER_IMPL_HPP

#include <LNOT/GeneralSolvers/NewtonTrustRegionSolver.hpp>
#include <LNOT/TRSSolvers/TruncatedConjugateGradient.hpp>
#include <LNOT/TRSSolvers/LanczosTRSSolver.hpp>
#include <LNOT/FloatingPoint/FPComparator.hpp>

namespace LNOT
{

//// explicit template instanciations ////

extern template class NewtonTrustRegionSolver< TruncatedConjugateGradient<float> >;
extern template class NewtonTrustRegionSolver< LanczosTRSSolver<float> >;

extern template class NewtonTrustRegionSolver< TruncatedConjugateGradient<double> >;
extern template class NewtonTrustRegionSolver< LanczosTRSSolver<double> >;

extern template class NewtonTrustRegionSolver< TruncatedConjugateGradient<long double> >;
extern template class NewtonTrustRegionSolver< LanczosTRSSolver<long double> >;

//// method implementations ////

template<typename TRSSolver> template<CSecondOrderOracle Oracle, typename ABool> 
void NewtonTrustRegionSolver<TRSSolver>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
{
	using AdlMath::sqrt;
	using AdlMath::floor;
	using AdlMath::log10;
	using AdlMath::abs;
	using AdlMath::pow;
	
	using Oracle_Size = typename Oracle::Size;
	
	const Oracle_Size size = oracle.getNDims();
	
	if (not m_gk or m_workCapacity < size)
	{
		m_workCapacity = size;
		m_gk     = std::make_unique<Scalar[]>(m_workCapacity);
		m_sk     = std::make_unique<Scalar[]>(m_workCapacity);
		m_xTrial = std::make_unique<Scalar[]>(m_workCapacity);
	}
	if (not solveInPlace) { std::fill(x, x + size, 0); }
	
	const auto Hk    = [&oracle] (const Scalar* d, Scalar* Hd)    -> void { oracle.getHessianProd(d, Hd);  };
	const auto invBk = [&oracle] (const Scalar* d, Scalar* invBd) -> void { oracle.applyPrecond(d, invBd); };
	
	m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk.get());
	
	m_fx = oracle.getValue();
	m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk.get(), size);
	
	Scalar delta = pow(Scalar(10.0), floor(log10(sqrt(Scalar(size)))));
	
	const Scalar relTol2 = m_relTol*m_relTol*m_squaredNormGrad;
	const Scalar absTol2 = m_absTol*m_absTol;
	
	const FPComparator<Scalar> cmp;
	const FPComparator<Scalar> cmpTr(m_trsSolver.getRelTolTR(), m_trsSolver.getAbsTolTR());
	
	if (m_out != nullptr) { fmt::println(m_out, "#Newton Trust region method\n#Iteration f(x) delta residual relative_tol absolute_tol"); }
	
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{				
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_fx, delta, m_squaredNormGrad, relTol2, absTol2); std::fflush(m_out); }
		if (m_squaredNormGrad < relTol2 or m_squaredNormGrad < absTol2) { m_info = Info::SUCCESS; break; }
		
		const Scalar normS = m_trsSolver.solve(Hk, invBk, m_gk.get(), size, delta, m_sk.get()); 
		
		m_innerIts.push_back(m_trsSolver.getIterations());
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_xTrial[i] = x[i] + m_sk[i]; }
		
		oracle.setCurrentPoint(m_xTrial.get());
		
		const Scalar fxTrial = oracle.getValue();
		const Scalar pred    = -m_trsSolver.getModelReduction();
		const Scalar ared    = m_fx - fxTrial;
		
		const bool isStepFeasible       = oracle.isFeasible();
		const bool isStepSuccessful     = TRMBase::isStepSuccessful(ared, pred, cmp);
		const bool isStepVerySuccessful = TRMBase::isStepVerySuccessful(ared, pred, cmp);
		const bool isStepAccepted       = TRMBase::isStepAccepted(ared, pred, cmp) and isStepFeasible;
		const bool normS_eq_delta       = cmpTr.isApproxEq(normS, delta);
		
		if      (not (isStepSuccessful and isStepFeasible)) { delta *= m_gammaDecrease; }
		else if (isStepVerySuccessful and normS_eq_delta)   { delta *= m_gammaIncrease; }
		
		delta = std::min(delta, m_deltaMax);
		
		if (isStepAccepted) 
		{ 
			std::copy(m_xTrial.get(), m_xTrial.get() + size, x); 
			oracle.getGradient(m_gk.get());	
			m_fx = fxTrial; 
			m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk.get(), size);
		} 
		else
		{
			oracle.setCurrentPoint(x);
		}
		
		if (not cmp.isDefPositive(delta)) { m_info = Info::BREAKDOWN; break; }    
	}
}
		
} // namespace LNOT

#endif // LNOT_NEWTON_TRUST_REGION_SOLVER_IMPL_HPP
