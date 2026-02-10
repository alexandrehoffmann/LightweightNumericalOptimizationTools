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

template<typename TRSSolver>
void NewtonTrustRegionSolver<TRSSolver>::clearWorkSpaceImpl()
{
	if (m_gk     != nullptr) { delete[] m_gk;     m_gk     = nullptr; }
	if (m_sk     != nullptr) { delete[] m_sk;     m_sk     = nullptr; }
	if (m_xTrial != nullptr) { delete[] m_xTrial; m_xTrial = nullptr; }
	m_workCapacity = 0;
}

template<typename TRSSolver> template<CSecondOrderOracle Oracle, typename ABool> 
void NewtonTrustRegionSolver<TRSSolver>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
{
	using std::sqrt;
	using std::floor;
	using std::log10;
	
	using Oracle_Size = typename Oracle::Size;
	
	const Oracle_Size size = oracle.getNDims();
	
	if (m_workCapacity < size)
	{
		clearWorkSpaceImpl();
		m_workCapacity = size;
		m_gk     = new Scalar[m_workCapacity];
		m_sk     = new Scalar[m_workCapacity];
		m_xTrial = new Scalar[m_workCapacity];
	}
	if (not solveInPlace) { std::fill(x, x + size, 0); }
	
	const auto Hk    = [&oracle] (const Scalar* d, Scalar* Hd)    -> void { oracle.getHessianProd(d, Hd);  };
	const auto invBk = [&oracle] (const Scalar* d, Scalar* invBd) -> void { oracle.applyPrecond(d, invBd); };
	
	m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk);
	
	m_fx = oracle.getValue();
	m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	
	Scalar delta = pow(10.0, floor(log10(sqrt(Scalar(size)))));
	
	const Scalar relTol2 = m_relTol*m_relTol*m_squaredNormGrad;
	const Scalar absTol2 = m_absTol*m_absTol;
	
	const FPComparator<Scalar> cmp;
	const FPComparator<Scalar> cmpTr(m_trsSolver.getRelTolTR(), m_trsSolver.getAbsTolTR());
	
	if (m_out != nullptr) { fmt::print(m_out, "#Newton Trust region method\n#Iteration f(x) delta residual relative_tol absolute_tol\n"); }
	
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{				
		if (m_out) { fmt::print(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e}\n", m_nIt, m_fx, delta, m_squaredNormGrad, relTol2, absTol2); }
		if (m_squaredNormGrad < relTol2 or m_squaredNormGrad < absTol2) { m_info = Info::SUCCESS; break; }
		
		m_trsSolver.solve(Hk, invBk, m_gk, size, delta, m_sk); 
		
		m_innerIts.push_back(m_trsSolver.getIterations());
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_xTrial[i] = x[i] + m_sk[i]; }
		
		oracle.setCurrentPoint(m_xTrial);
		
		const Scalar fxTrial   = oracle.getValue();
		const Scalar pred      = -m_trsSolver.getModelReduction();
		const Scalar ared      = m_fx - fxTrial;
		const Scalar normS     = BasicLinalg::norm(m_sk, size);
		
		const bool isStepFeasible       = oracle.isFeasible();
		const bool isStepSuccessful     = TRMBase::isStepSuccessful(ared, pred, cmp);
		const bool isStepVerySuccessful = TRMBase::isStepVerySuccessful(ared, pred, cmp);
		const bool isStepAccepted       = TRMBase::isStepAccepted(ared, pred, cmp) and isStepFeasible;
		const bool normS_eq_delta       = cmpTr.isApproxEq(normS, delta);
		
		if      (not (isStepSuccessful and isStepFeasible)) { delta *= m_gammaDecrease; }
		else if (isStepVerySuccessful and normS_eq_delta)   { delta *= m_gammaIncrease; }
		
		if (isStepAccepted) 
		{ 
			std::copy(m_xTrial, m_xTrial + size, x); 
			oracle.getGradient(m_gk);	
			m_fx = fxTrial; 
			m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
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
