#ifndef LNOT_LANCZOS_TRS_SOLVER_BASE_IMPL_HPP
#define LNOT_LANCZOS_TRS_SOLVER_BASE_IMPL_HPP

#include <LNOT/TRSSolvers/LanczosTRSSolverBase.hpp>
#include <LNOT/FloatingPoint/FPComparator.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <LNOT/misc/ZipLeadByView.hpp>

#include <span>
#include <ranges>

namespace LNOT
{

template<class Derived> 
LanczosTRSSolverBase<Derived>::LanczosTRSSolverBase(const Size maxIt, const Scalar tol, const Size maxItTr, const Scalar tolTr) 
	: Base(maxIt, tol, tolTr)
	, m_maxItTr(maxItTr) 
{
}

template<class Derived> 
void LanczosTRSSolverBase<Derived>::clearLanczosWorkSpace()
{	
	if (m_v      != nullptr) { delete[] m_v;      m_v      = nullptr; }
	if (m_p      != nullptr) { delete[] m_p;      m_p      = nullptr; }
	if (m_Bp     != nullptr) { delete[] m_Bp;     m_Bp     = nullptr; }
	if (m_Hp     != nullptr) { delete[] m_Hp;     m_Hp     = nullptr; }
	if (m_w      != nullptr) { delete[] m_w;      m_w      = nullptr; }
		
	m_workCapacity = 0;
}

template<class Derived> 
void LanczosTRSSolverBase<Derived>::allocateLanczosWorkSpace()
{
	m_v       = new Scalar[m_workCapacity];
	m_p       = new Scalar[m_workCapacity];
	m_Bp      = new Scalar[m_workCapacity];
	m_Hp      = new Scalar[m_workCapacity];
	m_w       = new Scalar[m_workCapacity];
}

template<class Derived> template<typename HesOp, typename PrecOp, typename ASize> 
auto LanczosTRSSolverBase<Derived>::solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Scalar& delta, Scalar* x) -> Scalar requires(AreHessianOps<HesOp,PrecOp>::value and IsSize<ASize>::value)
{
	using AdlMath::abs;
	using AdlMath::sqrt;
	
	const std::span<const Scalar> spanG(g, size);
	
	const auto minusG = spanG | std::views::transform([](const Scalar gi) -> Scalar
	{
		return -gi;
	});
	
	CRTP::derived().resizeWorkSpace(size);
	
	m_alpha.clear();
	m_beta.clear();
	
	std::fill(x, x + size, 0);
	Scalar precSqNormX{};
	
	m_lambda = 0;
	m_modelReduction = 0;
	
	bool isInterior = true;
	
	resetBv();
	
	std::fill(getBvPrev(), getBvPrev() + size, Scalar{});
	std::ranges::copy(minusG, getBvCurr());
	
	invB(getBvCurr(), m_v);
	
	const Scalar precNormR0 = sqrt(BasicLinalg::inner(getBvCurr(), m_v, size));
	
	m_beta.push_back(Scalar{});
	
	Scalar l_old = {};
	Scalar eta   = -precNormR0;
	
	BasicLinalg::scal(Scalar(1) / eta, size, getBvCurr());
	BasicLinalg::scal(Scalar(1) / eta, size, m_v);
	
	std::fill(m_p,  m_p  + size, Scalar{});
	std::fill(m_Bp, m_Bp + size, Scalar{});
	std::fill(m_Hp, m_Hp + size, Scalar{});
	
	m_precNormR = precNormR0;
	
	const Scalar relTol = m_relTol*m_precNormR;
	const Scalar delta2 = delta*delta;
	
	const FPComparator<Scalar> cmp;
	const FPComparator<Scalar> cmpTr(m_relTolTr, m_absTolTr);
	
	m_info = Info::FAILURE;
	if (m_out) { fmt::println(m_out, "#Preconditioned Lanczos TRS solver : \n#Iteration residual lambda relative_tol absolute_tol"); }
	
	// first try to solve within the trust region
	for (m_nIt=0; m_nIt!=m_maxIt and isInterior; ++m_nIt)
	{
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_precNormR, m_lambda, relTol, m_absTol); std::fflush(m_out); }
		if (m_precNormR < relTol or m_precNormR < m_absTol) { m_info = Info::SUCCESS; return sqrt(precSqNormX); }
		
		const Scalar* Bv      = getBvCurr();
		const Scalar* Bv_prev = getBvPrev();
		
		H(m_v, m_w);
		// m_w = Hv_{k}
		m_alpha.push_back( BasicLinalg::inner(m_v, m_w, size) );
		// solve T_k h_k = -\|r_0\|e_1
		const Scalar d = m_alpha.back() - m_beta.back()*l_old;
		const Scalar invD = Scalar(1) / d;
		
		if (not cmp.isDefPositive(d)) 
		{
			isInterior = false; 
			if (not solveBoundary(precNormR0, delta)) { m_info = Info::BREAKDOWN; return delta; }
		}
		else
		{
			#pragma omp simd
			for (Size i=0; i!=size; ++i) 
			{ 
				m_p[i]  = invD*(m_v[i] - m_beta.back()*m_p[i]); 
				m_Bp[i] = invD*( Bv[i] - m_beta.back()*m_Bp[i]); 
				m_Hp[i] = invD*(m_w[i] - m_beta.back()*m_Hp[i]); 
			} 
			m_modelReduction = BasicLinalg::updateModelReduction(m_modelReduction, x, g, eta, m_p, m_Hp, size);
			
			precSqNormX += 2*eta*BasicLinalg::inner(x, m_Bp, size) + eta*eta*BasicLinalg::inner(m_p, m_Bp, size);
			BasicLinalg::axpy(eta, m_p, size, x);
			
			if (cmpTr.isDefGreaterThan(precSqNormX, delta2))
			{
				isInterior = false;
				if (not solveBoundary(precNormR0, delta)) { m_info = Info::BREAKDOWN; return delta; }
			}
		}
		// resume Lanczos iteration
		BasicLinalg::axpbypz(-m_alpha.back(), Bv, -m_beta.back(), Bv_prev, size, m_w);
		reOrthonormalize(size, m_w);
		// m_w = \beta_{k}Bv_{k+1}
		invB(m_w, m_v);
		m_beta.push_back( sqrt(BasicLinalg::inner(m_w, m_v, size)) );
		
		const Scalar invBeta = Scalar(1) / m_beta.back();
		
		BasicLinalg::scal(invBeta, size, m_v);
		addBvNext(invBeta, m_w, size);
		// prepare next solve
		if (isInterior)
		{
			const Scalar l = invD*m_beta.back();
			
			eta         *= -l;
			l_old        = l;
			m_precNormR  = abs(eta); 
		}
		else
		{
			m_precNormR = abs(m_beta.back()*m_h.back());
		}
	}
	// The solution is located on the boundary of the trust region
	for (; m_nIt!=m_maxIt; ++m_nIt)
	{
		if (m_out) { fmt::print(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}\n", m_nIt, m_precNormR, m_lambda, relTol, m_absTol); std::fflush(m_out); }
		if (m_precNormR < relTol or m_precNormR < m_absTol) { m_info = Info::SUCCESS; break; }
		
		H(m_v, m_w);
		// m_w = Hv
		m_alpha.push_back( BasicLinalg::inner(m_v, m_w, size) );
		
		if (not solveBoundary(precNormR0, delta)) { m_info = Info::BREAKDOWN; return delta; }
		
		// resume Lanczos iteration
		BasicLinalg::axpbypz(-m_alpha.back(), getBvCurr(), -m_beta.back(), getBvPrev(), size, m_w);
		reOrthonormalize(size, m_w);
		
		// m_w = \beta_{k}Bv_{k+1}
		invB(m_w, m_v);
		m_beta.push_back( sqrt(BasicLinalg::inner(m_w, m_v, size)) );
		
		const Scalar invBeta = Scalar(1) / m_beta.back();
		BasicLinalg::scal(invBeta, size, m_v);
		addBvNext(invBeta, m_w, size);
		
		m_precNormR = abs(m_beta.back()*m_h.back());
	}
	
	std::fill(x, x + size, Scalar{});
	
	m_modelReduction = {};
	
	if constexpr (storeLanczosVectors)
	{
		for (Size i=0; i!=Size(m_h.size()); ++i)
		{
			invB(getBv(i+1), m_v);
			H(m_v, m_w);
			
			m_modelReduction = BasicLinalg::updateModelReduction(m_modelReduction, x, g, m_h[i], m_v, m_w, size);
			
			BasicLinalg::axpy(m_h[i], m_v, size, x);
		}
	}
	else // re-run Lanczos iteration to compute x = Vh;
	{
		const auto scaledG = spanG | std::views::transform([scale = Scalar(1) / precNormR0](const Scalar gi) -> Scalar
		{
			return scale*gi;
		});
		
		resetBv();
		
		std::fill(getBvPrev(), getBvPrev() + size, Scalar{});
		std::ranges::copy(scaledG, getBvCurr());
		
		invB(getBvCurr(), m_v);
		
		const std::span<const Scalar> spanH(m_h.data(), size_t(m_h.size() - 1));
		
		for (const auto [h_i, alpha_i, beta_im1, beta_i] : misc::zipLeadBy<0>(spanH, m_alpha, m_beta, m_beta | std::views::drop(1)))
		{
			H(m_v, m_w);
			m_modelReduction = BasicLinalg::updateModelReduction(m_modelReduction, x, g, h_i, m_v, m_w, size);
			
			BasicLinalg::axpy(h_i, m_v, size, x);
			
			BasicLinalg::axpbypz(-alpha_i, getBvCurr(), -beta_im1, getBvPrev(), size, m_w);
			reOrthonormalize(size, m_w);
			
			invB(m_w, m_v);
			
			const Scalar invNextBeta = Scalar(1) / beta_i;
			
			BasicLinalg::scal(invNextBeta, size, m_v);
			addBvNext(invNextBeta, m_w, size);
		}
		{
			H(m_v, m_w);
			m_modelReduction = BasicLinalg::updateModelReduction(m_modelReduction, x, g, m_h.back(), m_v, m_w, size);
		
			BasicLinalg::axpy(m_h.back(), m_v, size, x);
		}
	}
	
	m_nIt += Size(m_h.size());
	
	return delta;
}


template<class Derived> 
bool LanczosTRSSolverBase<Derived>::solveBoundary(const Scalar& gamma, const Scalar& delta)
{	
	using AdlMath::abs;
	using AdlMath::sqrt;
	
	namespace TridiagLDLt = BasicLinalg::Tridiag::LDLt;
	
	const Size   size      = Size(m_alpha.size());
	const Scalar delta2    = delta*delta;
	const Scalar norm1_T   = BasicLinalg::Tridiag::norm1(m_alpha.data(), std::next(m_beta.data()), size);
	
	const FPComparator<Scalar> cmp(m_relTolTr, m_absTolTr);

	m_h.reserve(m_alpha.capacity());
	m_invD.reserve(m_alpha.capacity());
	m_l.reserve(m_alpha.capacity()-1);
	
	m_h.resize(size);
	m_invD.resize(size);
	m_l.resize(size-1);
	
	// estimate lambdaMin and lambdaMax such that
	// lambda in [lambdaMin, lambdaMax]
	// then put our initial guess m_lambda, in [lambdaMin, lambdaMax].
	Scalar lambdaMin = std::max(Scalar{}, abs(gamma) / delta - norm1_T);
	Scalar lambdaMax = abs(gamma) / delta + norm1_T;
	m_lambda = std::clamp(m_lambda, lambdaMin, lambdaMax);
		
	for (size_t it=0;it!=m_maxItTr;++it)
	{	
		const bool ldltSuccess = TridiagLDLt::compute(m_alpha.data(), std::next(m_beta.data()), size, m_lambda, m_invD.data(), m_l.data());
		const bool isSpd       = std::ranges::all_of(m_invD, [&cmp](const Scalar& invDi) -> bool { return cmp.isDefPositive(invDi); });
		const bool cholSuccess = ldltSuccess and isSpd;
		
		if (cholSuccess)
		{
			TridiagLDLt::solveUnit(m_invD.data(), m_l.data(), size, -gamma, m_h.data());
			const Scalar sqNormH = BasicLinalg::squaredNorm(m_h.data(), size);
			const Scalar   normH = sqrt(sqNormH);
			
			if (cmp.isApproxEq(normH, delta))         { return true;  }
			if (cmp.isApproxEq(lambdaMax, lambdaMin)) { return false; } // empty interval
			
			TridiagLDLt::solveInplaceLower(m_l.data(), size, m_h.data());
			
			const Scalar sqnorm_w    = BasicLinalg::weightedSquaredNorm(m_h.data(), m_invD.data(), size);
			const Scalar deltaLambda = ((normH - delta) / delta)*(sqNormH / sqnorm_w);
			
			if (cmp.isDefLessThan(sqNormH, delta2)) { lambdaMax = m_lambda; }
			else                                    { lambdaMin = m_lambda; }
			
			const Scalar gap = (lambdaMax - lambdaMin)*Scalar(0.005);
			m_lambda = std::clamp(m_lambda + deltaLambda, lambdaMin + gap, lambdaMax - gap);
		}
		else
		{
			lambdaMin = m_lambda;
			m_lambda  = Scalar(0.5)*(lambdaMax + lambdaMin); 
		}
	}
	// hard case cannot occur, we did not converge due to numerical errors
	return false;
}

} // namespace LNOT

#endif // LNOT_LANCZOS_TRS_SOLVER_BASE_IMPL_HPP
