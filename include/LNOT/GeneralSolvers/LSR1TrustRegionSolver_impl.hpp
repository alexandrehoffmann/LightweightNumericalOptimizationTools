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

template<typename TRSSolver>
void LSR1TrustRegionSolver<TRSSolver>::clearWorkSpaceImpl()
{
	if (m_gk     != nullptr) { delete[] m_gk;     m_gk     = nullptr; }
	if (m_gkp1   != nullptr) { delete[] m_gkp1;   m_gkp1   = nullptr; }
	if (m_xTrial != nullptr) { delete[] m_xTrial; m_xTrial = nullptr; }
	if (m_P      != nullptr) { delete[] m_P;      m_P      = nullptr; }
	if (m_Y      != nullptr) { delete[] m_Y;      m_Y      = nullptr; }
	if (m_S      != nullptr) { delete[] m_S;      m_S      = nullptr; }
	m_workCapacity = 0;
}

template<typename TRSSolver> template<CFirstOrderOracle Oracle, typename ABool> 
void LSR1TrustRegionSolver<TRSSolver>::solveImpl(Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
{
	using std::sqrt;
	using std::floor;
	using std::log10;
	using std::pow;
	using std::abs;
	using std::isfinite;
	
	using Oracle_Size         = typename Oracle::Size;
	using CircularBuffer_Size = typename CircularBuffer<Scalar>::size_type;
	
	const Scalar sr1DropTol = sqrt( std::numeric_limits<Scalar>::epsilon() );
	const Oracle_Size size = oracle.getNDims();
	
	if (m_workCapacity < size)
	{
		clearWorkSpaceImpl();
		m_workCapacity = size;
		m_gk     = new Scalar[m_workCapacity];
		m_gkp1   = new Scalar[m_workCapacity];
		m_xTrial = new Scalar[m_workCapacity];
		m_P      = new Scalar[m_workCapacity*m_memory];
		m_Y      = new Scalar[m_workCapacity*m_memory];
		m_S      = new Scalar[m_workCapacity*m_memory];
	}
	if constexpr (not solveInPlace) { std::fill(x, x + size, 0); }	
	
	CircularBuffer<Scalar> invRho(m_memory);
	std::vector<bool>      isVectorKept(m_memory);
	
	Size curr_idx{};
	auto BkOp = [this, &curr_idx, &invRho, &isVectorKept, size](const Scalar* d, Scalar* Bd) -> void
	{
		// from https://home.cs.colorado.edu/~richard/lu_dissertation.pdf
		const Size prev_idx = (curr_idx == 0) ? m_memory-1 : curr_idx-1;
		const Scalar* skm1 = m_S + prev_idx*size;
		const Scalar* ykm1 = m_Y + prev_idx*size;
		//~ const Scalar protoGamma0 = invRho.empty() ? 1 : BasicLinalg::inner(ykm1, skm1, size) / BasicLinalg::squaredNorm(skm1, size);
		//~ const Scalar protoGamma0 = invRho.empty() ? 1 : BasicLinalg::norm(ykm1, size) / BasicLinalg::norm(skm1, size);
		const Scalar protoGamma0 = invRho.empty() ? 1 : BasicLinalg::squaredNorm(ykm1, size) / BasicLinalg::inner(ykm1, skm1, size);
		const Scalar gamma0 = isfinite(protoGamma0) ? protoGamma0 : 1;
		
		#pragma omp simd
		for (Size k=0; k!=size; ++k) { Bd[k] = gamma0*d[k]; }
		invRho.foreach([this, &isVectorKept, size, d, Bd](const CircularBuffer_Size i, const Scalar& invRho_i) { if (isVectorKept[i])
		{
			const Scalar* pi   = m_P + i*size;
			const Scalar alpha = BasicLinalg::inner(pi, d, size) / invRho_i;
			BasicLinalg::axpy(alpha, pi, size, Bd);
		}});
	};

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
	
	if (m_out != nullptr) { fmt::print(m_out, "#L-SR1 Trust region method\n#Iteration f(x) delta usedVectors residual relative_tol absolute_tol\n"); }
	m_info = Info::FAILURE;
	for (m_nIt=0;m_nIt!=m_maxIt; ++m_nIt)
	{				
		if (m_out) { fmt::print(m_out, "{} {:10.2e} {:10.2e} {} {:10.2e} {:10.2e} {:10.2e}\n", m_nIt, m_fx, delta, std::ranges::count(isVectorKept, true), m_squaredNormGrad, relTol2, absTol2); }
		if (m_squaredNormGrad < relTol2 or m_squaredNormGrad < absTol2) { m_info = Info::SUCCESS; break; }
		// building the Bk matrix from the last saved vectors
		// c.f. https://optimization-online.org/wp-content/uploads/2015/10/5167.pdf
		std::fill(isVectorKept.begin(), isVectorKept.end(), false);
		invRho.foreach([this, &BkOp, &isVectorKept, &sr1DropTol, size](CircularBuffer_Size i, Scalar& invRho_i)
		{
			const Scalar* si = m_S + i*size;
			const Scalar* yi = m_Y + i*size;
			Scalar* pi = m_P + i*size;
			BkOp(si, pi);
			#pragma omp simd
			for (Size k=0; k!=size; ++k) { pi[k] = yi[k] - pi[k]; }
			
			invRho_i = BasicLinalg::inner(pi, si, size);
			isVectorKept[i] = abs(invRho_i) > sr1DropTol*BasicLinalg::norm(pi, size)*BasicLinalg::norm(si, size);
		});
		// now resume as usual TR method
		const Scalar normS = m_trsSolver.solve(BkOp, m_gk, size, delta, m_S + curr_idx*size); 
		
		m_innerIts.push_back(m_trsSolver.getIterations());
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_xTrial[i] = x[i] + m_S[i + curr_idx*size]; }
		
		oracle.setCurrentPoint(m_xTrial);
		oracle.getGradient(m_gkp1);
		
		const Scalar fxTrial   = oracle.getValue();
		const Scalar pred      = -m_trsSolver.getModelReduction();
		const Scalar ared      = m_fx - fxTrial;
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Y[i + curr_idx*size] = m_gkp1[i] - m_gk[i]; } // u_k = y_k - Bks_k
		invRho.push(0);
		
		++curr_idx; 
		if (curr_idx == m_memory) { curr_idx = 0; }
		
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

#endif // LNOT_LSR1_TRUST_REGION_SOLVER_IMPL_HPP
