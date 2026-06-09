#ifndef LNOT_SR1_TRUST_REGION_SOLVER_IMPL_HPP
#define LNOT_SR1_TRUST_REGION_SOLVER_IMPL_HPP

#include <LNOT/GeneralSolvers/SR1TrustRegionSolver.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/TRSSolvers/TruncatedConjugateGradient.hpp>
#include <LNOT/TRSSolvers/LanczosTRSSolver.hpp>

namespace LNOT
{

//// explicit template instanciations ////

extern template class SR1TrustRegionSolver< TruncatedConjugateGradient<float> >;
extern template class SR1TrustRegionSolver< LanczosTRSSolver<float> >;

extern template class SR1TrustRegionSolver< TruncatedConjugateGradient<double> >;
extern template class SR1TrustRegionSolver< LanczosTRSSolver<double> >;

extern template class SR1TrustRegionSolver< TruncatedConjugateGradient<long double> >;
extern template class SR1TrustRegionSolver< LanczosTRSSolver<long double> >;

//// method implementations ////

template<typename TRSSolver, typename ConvergenceCriterion> template<CFirstOrderOracle Oracle, typename ABool> 
void SR1TrustRegionSolver<TRSSolver, ConvergenceCriterion>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
{
	using AdlMath::sqrt;
	using AdlMath::floor;
	using AdlMath::log10;
	using AdlMath::abs;
	using AdlMath::pow;
	
	using Oracle_Size = typename Oracle::Size;
	
	const Scalar sr1DropTol = sqrt( NumTraits<Scalar>::epsilon );
		
	const Oracle_Size size = oracle.getNDims();
	
	if (not m_gk or m_workCapacity < size)
	{
		m_workCapacity = size;
		m_gk     = std::make_unique<Scalar[]>(m_workCapacity);
		m_gkp1   = std::make_unique<Scalar[]>(m_workCapacity);
		m_uk     = std::make_unique<Scalar[]>(m_workCapacity);
		m_sk     = std::make_unique<Scalar[]>(m_workCapacity);
		m_Bk     = std::make_unique<Scalar[]>(m_workCapacity*m_workCapacity);
		m_xTrial = std::make_unique<Scalar[]>(m_workCapacity);
	}
	if (not solveInPlace) { std::fill(x, x + size, 0); }	
	
	std::fill(m_Bk.get(), m_Bk.get() + size*size, 0);
	for (Size i=0; i!=size; ++i) { m_Bk[i + i*size] = 1; }
	
	SymmetricDenseMatrixOp<Scalar> BkOp(m_Bk.get(), size);
	
	m_innerIts.clear();
	
	Scalar delta = pow(Scalar(10.0), floor(log10(sqrt(Scalar(size)))));
	
	constexpr ConvergenceCriterion criterion;
	constexpr FPComparator<Scalar> cmp;
	const     FPComparator<Scalar> cmpTr(m_trsSolver.getRelTolTR(), m_trsSolver.getAbsTolTR());
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk.get());
	
	m_fx = oracle.getValue();
	m_residual = criterion.getResidual(m_gk.get(), size); 
	      
	const Scalar relTol   = criterion.getRelTol(m_relTol, m_residual);
	const Scalar absTol   = criterion.getAbsTol(m_absTol);
	
	if (m_out != nullptr) { fmt::println(m_out, "#SR1 Trust region method\n#Iteration f(x) delta residual relative_tol absolute_tol"); }
	
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{				
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_fx, delta, m_residual, relTol, absTol); std::fflush(m_out); }
		if (m_residual < relTol or m_residual < absTol) { m_info = Info::SUCCESS; break; }
		
		const Scalar normS = m_trsSolver.solve(BkOp, m_gk.get(), size, delta, m_sk.get()); 
		
		m_innerIts.push_back(m_trsSolver.getIterations());
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_xTrial[i] = x[i] + m_sk[i]; }
		
		oracle.setCurrentPoint(m_xTrial.get());
		oracle.getGradient(m_gkp1.get());
		BkOp(m_sk.get(), m_uk.get());
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_uk[i] = m_gkp1[i] - m_gk[i] - m_uk[i]; } // u_k = y_k - Bks_k
		
		const Scalar invRho = BasicLinalg::inner(m_sk.get(), m_uk.get(), size);
		
		if (abs(invRho) > sr1DropTol*BasicLinalg::squaredNorm(m_sk.get(), size)*BasicLinalg::squaredNorm(m_uk.get(), size))
		{
			BasicLinalg::symRk1Update(BkOp.getStorageOrder(), BkOp.getUplo(), Scalar(1) / invRho, m_uk.get(), size, m_Bk.get());
		}
		
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
			std::copy(m_gkp1.get(),   m_gkp1.get()   + size, m_gk.get()); 
			m_fx = fxTrial; 
			m_residual = criterion.getResidual(m_gk.get(), size);
		} 
		else
		{
			oracle.setCurrentPoint(x);
		}
		
		if (not cmp.isDefPositive(delta)) { m_info = Info::BREAKDOWN; break; }    
	}
}

} // namespace LNOT

#endif // LNOT_SR1_TRUST_REGION_SOLVER_IMPL_HPP
