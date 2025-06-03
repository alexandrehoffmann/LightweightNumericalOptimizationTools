#ifndef LNOT_LSR1_TRUST_REGION_SOLVER_IMPL_HPP
#define LNOT_LSR1_TRUST_REGION_SOLVER_IMPL_HPP

#include <LNOT/GeneralSolvers/LSR1TrustRegionSolver.hpp>
#include <LNOT/Containers/CircularBuffer.hpp>

namespace LNOT
{

template<typename TRSSolver>
void LSR1TrustRegionSolver<TRSSolver>::clearWorkSpace()
{
	if (m_gk     != nullptr) { delete[] m_gk;     m_gk     = nullptr; }
	if (m_gkp1   != nullptr) { delete[] m_gkp1;   m_gkp1   = nullptr; }
	if (m_Bksk   != nullptr) { delete[] m_Bksk;   m_Bksk   = nullptr; }
	if (m_xTrial != nullptr) { delete[] m_xTrial; m_xTrial = nullptr; }
	if (m_U      != nullptr) { delete[] m_U;      m_U      = nullptr; }
	if (m_Y      != nullptr) { delete[] m_Y;      m_Y      = nullptr; }
	if (m_S      != nullptr) { delete[] m_S;      m_S      = nullptr; }
	Base::m_workCapacity = 0;
}

template<typename TRSSolver> template<FirstOrderOracle_concept Oracle, bool solveInPlace> 
void LSR1TrustRegionSolver<TRSSolver>::solve_impl(Oracle& oracle, std::bool_constant<solveInPlace>, Scalar* x)
{
	using CircularBuffer_size = typename CircularBuffer<Scalar>::size_type;
	
	const Scalar sr1DropTol = std::sqrt( std::numeric_limits<Scalar>::epsilon() );
	const Size size = oracle.getNDims();
	
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_gk     = new Scalar[Base::m_workCapacity];
	  m_gkp1   = new Scalar[Base::m_workCapacity];
	  m_Bksk   = new Scalar[Base::m_workCapacity];
	  m_xTrial = new Scalar[Base::m_workCapacity];
	  m_U      = new Scalar[Base::m_workCapacity*m_memory];
	  m_Y      = new Scalar[Base::m_workCapacity*m_memory];
	  m_S      = new Scalar[Base::m_workCapacity*m_memory];
	}
	if constexpr (not solveInPlace) { std::fill(x, x + size, 0); }	
	
	CircularBuffer<Scalar> invRho(m_memory);
	std::vector<bool>      isVectorKept(m_memory);

	auto BkOp = [this, &invRho, &isVectorKept, size](const Scalar* d, Scalar* Bd) -> void
	{
		std::copy(d, d + size, Bd);
		invRho.foreach([this, &isVectorKept, size, d, Bd](const CircularBuffer_size i, const Scalar& invRho_i) { if (isVectorKept[i])
		{
			const Scalar alpha = BasicLinalg::inner(m_U + i*size, d, size) / invRho_i;
			BasicLinalg::axpy(alpha, m_U + i*size, size, Bd);
		}});
	};

	Base::m_innerIts.clear();
	
	oracle.setCurrentPoint(x);
	oracle.getGradient(m_gk);
	
	Base::m_fx = oracle.getValue();
	Base::m_squaredNormGrad = BasicLinalg::squaredNorm(m_gk, size);
	
	Scalar delta = std::pow(10.0, std::floor(std::log10(std::sqrt(Scalar(size)))));
	
	const Scalar tol2 = Base::m_tol*Base::m_tol*std::max(Scalar(1), Base::m_squaredNormGrad);
	
	Size curr_idx = 0;
	
	if (Base::m_out != nullptr) { fmt::print(Base::m_out, "#L-SR1 Trust region method\n#Iteration f(x) delta usedVectors residual tol\n"); }
	Base::m_info = Info::FAILURE;
	for (Base::m_nIt=0;Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{				
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e} {} {:10.2e} {:10.2e}\n", Base::m_nIt, Base::m_fx, delta, std::ranges::count(isVectorKept, true), Base::m_squaredNormGrad, tol2); }
		if (Base::m_squaredNormGrad < tol2) { Base::m_info = Info::SUCCESS; break; }
		// building the Bk matrix from the last saved vectors
		std::fill(isVectorKept.begin(), isVectorKept.end(), false);
		invRho.foreach([this, &invRho, &isVectorKept, &sr1DropTol, &BkOp, size](CircularBuffer_size k, Scalar& invRho_k)
		{
			// Computes B_k s_k
			BkOp(m_S + k*size, m_Bksk);
			// Computes u_k := y_k - B_k s_k
			#pragma omp simd
			for (Size i=0; i!=size; ++i) { m_U[i + k*size] = m_Y[i + k*size] - m_Bksk[i]; }
			invRho_k = BasicLinalg::inner(m_U + k*size, m_S + k*size, size);
			isVectorKept[k] = (std::abs(invRho_k) > sr1DropTol*BasicLinalg::norm(m_S + k*size, size)*BasicLinalg::norm(m_U + k*size, size));
		});
		// now resume as usual TR method
		m_trsSolver.solve(BkOp, m_gk, size, delta, m_S + curr_idx*size); 
		
		Base::m_innerIts.push_back(m_trsSolver.getIterations());
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_xTrial[i] = x[i] + m_S[i + curr_idx*size]; }
		
		oracle.setCurrentPoint(m_xTrial);
		oracle.getGradient(m_gkp1);
		
		const Scalar fxTrial = oracle.getValue();
		const Scalar pred     = -m_trsSolver.getModelReduction();
		const Scalar ared     = Base::m_fx - fxTrial;
		const Scalar normS    = BasicLinalg::norm(m_S + curr_idx*size, size);
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Y[i + curr_idx*size] = m_gkp1[i] - m_gk[i]; } // u_k = y_k - Bks_k
		invRho.push(0);
		++curr_idx; if (curr_idx == m_memory) { curr_idx = 0; }
		
		const bool isStepFeasible       = std::isfinite(fxTrial);
		const bool isStepSuccessful     = ared > 0 and ared > pred*TRSBase::m_etaSuccessful;
		const bool isStepVerySuccessful = ared > 0 and ared > pred*TRSBase::m_etaVerySuccessful;
		const bool isStepAccepted       = isStepFeasible and ared > 0 and ared > pred*TRSBase::m_etaAccept;
		const bool normS_eq_delta       = (std::abs(normS - delta) <= m_trsSolver.getTolTR()*std::max(normS, delta));
		
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

#endif // LNOT_LSR1_TRUST_REGION_SOLVER_IMPL_HPP
