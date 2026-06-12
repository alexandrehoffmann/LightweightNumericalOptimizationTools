#ifndef LNOT_LSR1_TRUST_REGION_SOLVER_IMPL_HPP
#define LNOT_LSR1_TRUST_REGION_SOLVER_IMPL_HPP

#include <LNOT/GeneralSolvers/LSR1TrustRegionSolver.hpp>
#include <LNOT/Containers/CircularBuffer.hpp>
#include <LNOT/TRSSolvers/TruncatedConjugateGradient.hpp>
#include <LNOT/TRSSolvers/LanczosTRSSolver.hpp>

namespace LNOT
{
	
//// explicit template instanciations ////

extern template class LSR1TrustRegionSolver< TruncatedConjugateGradient<float> >;
extern template class LSR1TrustRegionSolver< LanczosTRSSolver<float> >;

extern template class LSR1TrustRegionSolver< TruncatedConjugateGradient<double> >;
extern template class LSR1TrustRegionSolver< LanczosTRSSolver<double> >;

extern template class LSR1TrustRegionSolver< TruncatedConjugateGradient<long double> >;
extern template class LSR1TrustRegionSolver< LanczosTRSSolver<long double> >;

//// method implementations ////

template<typename TRSSolver, typename ConvergenceCriterion> template<CFirstOrderOracle Oracle, typename ABool> 
void LSR1TrustRegionSolver<TRSSolver, ConvergenceCriterion>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
{
	using AdlMath::sqrt;
	using AdlMath::floor;
	using AdlMath::log10;
	using AdlMath::pow;
	using AdlMath::abs;
	using AdlMath::isfinite;
	
	using Oracle_Size         = typename Oracle::Size;
	using CircularBuffer_Size = typename CircularBuffer<Scalar>::size_type;
	
	const Scalar sr1DropTol = sqrt( NumTraits<Scalar>::epsilon );
	const Oracle_Size size = oracle.getNDims();
	
	if (not m_gk or m_workCapacity < size)
	{
		m_workCapacity = size;
		m_gk     = std::make_unique<Scalar[]>(m_workCapacity);
		m_gkp1   = std::make_unique<Scalar[]>(m_workCapacity);
		m_xTrial = std::make_unique<Scalar[]>(m_workCapacity);
		m_P      = std::make_unique<Scalar[]>(m_workCapacity*m_memory);
		m_Y      = std::make_unique<Scalar[]>(m_workCapacity*m_memory);
		m_S      = std::make_unique<Scalar[]>(m_workCapacity*m_memory);
	}
	if constexpr (not solveInPlace) { std::fill(x, x + size, 0); }	
	
	CircularBuffer<Scalar> invRho(m_memory);
	std::vector<bool>      isVectorKept(m_memory);
	
	Size curr_idx{};
	auto BkOp = [this, &curr_idx, &invRho, &isVectorKept, size](const Scalar* d, Scalar* Bd) -> void
	{
		// from https://home.cs.colorado.edu/~richard/lu_dissertation.pdf
		const Size prev_idx = (curr_idx == 0) ? m_memory-1 : curr_idx-1;
		const Scalar* skm1 = m_S.get() + prev_idx*size;
		const Scalar* ykm1 = m_Y.get() + prev_idx*size;
		//~ const Scalar protoGamma0 = invRho.empty() ? 1 : BasicLinalg::inner(ykm1, skm1, size) / BasicLinalg::squaredNorm(skm1, size);
		//~ const Scalar protoGamma0 = invRho.empty() ? 1 : BasicLinalg::norm(ykm1, size) / BasicLinalg::norm(skm1, size);
		const Scalar protoGamma0 = invRho.empty() ? 1 : BasicLinalg::squaredNorm(ykm1, size) / BasicLinalg::inner(ykm1, skm1, size);
		const Scalar gamma0 = isfinite(protoGamma0) ? protoGamma0 : 1;
		
		#pragma omp simd
		for (Size k=0; k!=size; ++k) { Bd[k] = gamma0*d[k]; }
		invRho.foreach([this, &isVectorKept, size, d, Bd](const CircularBuffer_Size i, const Scalar& invRho_i) { if (isVectorKept[i])
		{
			const Scalar* pi   = m_P.get() + i*size;
			const Scalar alpha = BasicLinalg::inner(pi, d, size) / invRho_i;
			BasicLinalg::axpy(alpha, pi, size, Bd);
		}});
	};

	m_innerIts.clear();
		
	Scalar delta = pow(Scalar(10.0), floor(log10(sqrt(Scalar(size)))));
	
	constexpr FPComparator<Scalar> cmp;
	const     FPComparator<Scalar> cmpTr(m_trsSolver.getRelTolTR(), m_trsSolver.getAbsTolTR());
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk.get());
	
	m_criterion.init(m_gk.get(), size);
	
	m_fx = oracle.getValue();
	m_residual = m_criterion.getResidual(m_gk.get(), size); 
	      
	const Scalar relTol   = m_criterion.getRelTol(m_relTol);
	const Scalar absTol   = m_criterion.getAbsTol(m_absTol);
	
	if (m_out != nullptr) { fmt::println(m_out, "#L-SR1 Trust region method\n#Iteration f(x) delta usedVectors residual relative_tol absolute_tol"); }
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{				
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_fx, delta, std::ranges::count(isVectorKept, true), m_residual, relTol, absTol); std::fflush(m_out); }
		if (m_residual < relTol or m_residual < absTol) { m_info = Info::SUCCESS; break; }
		// building the Bk matrix from the last saved vectors
		// c.f. https://optimization-online.org/wp-content/uploads/2015/10/5167.pdf
		std::fill(std::ranges::begin(isVectorKept), std::ranges::end(isVectorKept), false);
		invRho.foreach([this, &BkOp, &isVectorKept, &sr1DropTol, size](CircularBuffer_Size i, Scalar& invRho_i)
		{
			const Scalar* si = m_S.get() + i*size;
			const Scalar* yi = m_Y.get() + i*size;
			Scalar* pi = m_P.get() + i*size;
			BkOp(si, pi);
			#pragma omp simd
			for (Size k=0; k!=size; ++k) { pi[k] = yi[k] - pi[k]; }
			
			invRho_i = BasicLinalg::inner(pi, si, size);
			isVectorKept[i] = abs(invRho_i) > sr1DropTol*BasicLinalg::norm(pi, size)*BasicLinalg::norm(si, size);
		});
		// now resume as usual TR method
		const Scalar normS = m_trsSolver.solve(BkOp, m_gk.get(), size, delta, m_S.get() + curr_idx*size); 
		
		m_innerIts.push_back(m_trsSolver.getIterations());
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_xTrial[i] = x[i] + m_S[i + curr_idx*size]; }
		
		oracle.setCurrentPoint(m_xTrial.get());
		oracle.getGradient(m_gkp1.get());
		
		const Scalar fxTrial   = oracle.getValue();
		const Scalar pred      = -m_trsSolver.getModelReduction();
		const Scalar ared      = m_fx - fxTrial;
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Y[i + curr_idx*size] = m_gkp1[i] - m_gk[i]; } // u_k = y_k - Bks_k
		invRho.push(Scalar{});
		
		++curr_idx; 
		if (curr_idx == m_memory) { curr_idx = 0; }
		
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
			m_residual = m_criterion.getResidual(m_gk.get(), size);
		} 
		else
		{
			oracle.setCurrentPoint(x);
		}
		
		if (not cmp.isDefPositive(delta)) { m_info = Info::BREAKDOWN; break; }    
	}

}
	
} // namespace LNOT

#endif // LNOT_LSR1_TRUST_REGION_SOLVER_IMPL_HPP
