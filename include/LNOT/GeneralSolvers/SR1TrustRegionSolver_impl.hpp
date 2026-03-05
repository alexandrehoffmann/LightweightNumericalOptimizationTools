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

template<typename TRSSolver>
void SR1TrustRegionSolver<TRSSolver>::clearWorkSpaceImpl()
{
	if (m_gk     != nullptr) { delete[] m_gk;     m_gk     = nullptr; }
	if (m_gkp1   != nullptr) { delete[] m_gkp1;   m_gkp1   = nullptr; }
	if (m_uk     != nullptr) { delete[] m_uk;     m_uk     = nullptr; }
	if (m_sk     != nullptr) { delete[] m_sk;     m_sk     = nullptr; }
	if (m_Bk     != nullptr) { delete[] m_Bk;     m_Bk     = nullptr; }
	if (m_xTrial != nullptr) { delete[] m_xTrial; m_xTrial = nullptr; }
	m_workCapacity = 0;
}

template<typename TRSSolver> template<CFirstOrderOracle Oracle, typename ABool> 
void SR1TrustRegionSolver<TRSSolver>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
{
	using AdlMath::sqrt;
	using AdlMath::floor;
	using AdlMath::log10;
	using AdlMath::abs;
	using AdlMath::pow;
	
	using Oracle_Size = typename Oracle::Size;
	
	const Scalar sr1DropTol = sqrt( NumTraits<Scalar>::epsilon );
		
	const Oracle_Size size = oracle.getNDims();
	
	if (m_workCapacity < size)
	{
		clearWorkSpaceImpl();
		m_workCapacity = size;
		m_gk     = new Scalar[m_workCapacity];
		m_gkp1   = new Scalar[m_workCapacity];
		m_uk     = new Scalar[m_workCapacity];
		m_sk     = new Scalar[m_workCapacity];
		m_Bk     = new Scalar[m_workCapacity*m_workCapacity];
		m_xTrial = new Scalar[m_workCapacity];
	}
	if (not solveInPlace) { std::fill(x, x + size, 0); }	
	
	std::fill(m_Bk, m_Bk + size*size, 0);
	for (Size i=0; i!=size; ++i) { m_Bk[i + i*size] = 1; }
	
	SymmetricDenseMatrixOp<Scalar> BkOp(m_Bk, size);
	
	m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk);
	
	m_fx = oracle.getValue();
	m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	
	Scalar delta = pow(Scalar(10.0), floor(log10(sqrt(Scalar(size)))));
	
	const Scalar relTol2 = m_relTol*m_relTol*m_squaredNormGrad;
	const Scalar absTol2 = m_absTol*m_absTol;
	
	const FPComparator<Scalar> cmp;
	const FPComparator<Scalar> cmpTr(m_trsSolver.getRelTolTR(), m_trsSolver.getAbsTolTR());
	
	if (m_out != nullptr) { fmt::println(m_out, "#SR1 Trust region method\n#Iteration f(x) delta residual relative_tol absolute_tol"); }
	
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{				
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_fx, delta, m_squaredNormGrad, relTol2, absTol2); std::fflush(m_out); }
		if (m_squaredNormGrad < relTol2 or m_squaredNormGrad < absTol2) { m_info = Info::SUCCESS; break; }
		
		const Scalar normS = m_trsSolver.solve(BkOp, m_gk, size, delta, m_sk); 
		
		m_innerIts.push_back(m_trsSolver.getIterations());
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_xTrial[i] = x[i] + m_sk[i]; }
		
		oracle.setCurrentPoint(m_xTrial);
		oracle.getGradient(m_gkp1);
		BkOp(m_sk, m_uk);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_uk[i] = m_gkp1[i] - m_gk[i] - m_uk[i]; } // u_k = y_k - Bks_k
		
		const Scalar invRho = BasicLinalg::inner(m_sk, m_uk, size);
		
		if (abs(invRho) > sr1DropTol*BasicLinalg::squaredNorm(m_sk, size)*BasicLinalg::squaredNorm(m_uk, size))
		{
			BasicLinalg::symRk1Update(BkOp.getStorageOrder(), BkOp.getUplo(), Scalar(1) / invRho, m_uk, size, m_Bk);
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
		
		if (isStepAccepted) 
		{ 
			std::copy(m_xTrial, m_xTrial + size, x); 
			std::copy(m_gkp1,   m_gkp1   + size, m_gk); 
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

#endif // LNOT_SR1_TRUST_REGION_SOLVER_IMPL_HPP
