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
void LanczosTRSSolverBase<Derived>::allocateLanczosWorkSpace()
{
	m_v  = std::make_unique<Scalar[]>(m_workCapacity);
	m_p  = std::make_unique<Scalar[]>(m_workCapacity);
	m_Bp = std::make_unique<Scalar[]>(m_workCapacity);
	m_Hp = std::make_unique<Scalar[]>(m_workCapacity);
	m_w  = std::make_unique<Scalar[]>(m_workCapacity);
}

template<class Derived> template<typename HesOp, typename PrecOp, typename ASize> 
auto LanczosTRSSolverBase<Derived>::solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Scalar& delta, Scalar* x) -> Scalar requires(areHessianOps<HesOp,PrecOp> and isSize<ASize>)
{
	using AdlMath::abs;
	using AdlMath::sqrt;
	
	const std::span<const Scalar> spanG(g, size);
	
	const auto minusG = spanG | std::views::transform([](const Scalar gi) -> Scalar
	{
		return -gi;
	});
	
	this->derived().resizeWorkSpace(size);
	
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
	
	invB(getBvCurr(), m_v.get());
	
	const Scalar precNormR0 = sqrt(BasicLinalg::inner(getBvCurr(), m_v.get(), size));
	
	m_beta.push_back(Scalar{});
	
	Scalar l_old = {};
	Scalar eta   = -precNormR0;
	
	BasicLinalg::scal(Scalar(1) / eta, size, getBvCurr());
	BasicLinalg::scal(Scalar(1) / eta, size, m_v.get());
	
	std::fill(m_p.get(),  m_p.get()  + size, Scalar{});
	std::fill(m_Bp.get(), m_Bp.get() + size, Scalar{});
	std::fill(m_Hp.get(), m_Hp.get() + size, Scalar{});
	
	m_precNormR = precNormR0;
	
	const Scalar relTol = m_relTol*m_precNormR;
	const Scalar delta2 = delta*delta;
	
	constexpr FPComparator<Scalar> cmp;
	const     FPComparator<Scalar> cmpTr(m_relTolTr, m_absTolTr);
	
	m_info = Info::FAILURE;
	if (m_out) { fmt::println(m_out, "#Preconditioned Lanczos TRS solver : \n#Iteration residual lambda relative_tol absolute_tol"); }
	
	// first try to solve within the trust region
	for (m_nIt=0; m_nIt!=m_maxIt and isInterior; ++m_nIt)
	{
		if (m_out) { fmt::println(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", m_nIt, m_precNormR, m_lambda, relTol, m_absTol); std::fflush(m_out); }
		if (m_precNormR < relTol or m_precNormR < m_absTol) { m_info = Info::SUCCESS; return sqrt(precSqNormX); }
		
		const Scalar* Bv      = getBvCurr();
		const Scalar* Bv_prev = getBvPrev();
		
		H(m_v.get(), m_w.get());
		// m_w = Hv_{k}
		m_alpha.push_back( BasicLinalg::inner(m_v.get(), m_w.get(), size) );
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
			m_modelReduction = BasicLinalg::updateModelReduction(m_modelReduction, x, g, eta, m_p.get(), m_Hp.get(), size);
			
			precSqNormX += 2*eta*BasicLinalg::inner(x, m_Bp.get(), size) + eta*eta*BasicLinalg::inner(m_p.get(), m_Bp.get(), size);
			BasicLinalg::axpy(eta, m_p.get(), size, x);
			
			if (cmpTr.isDefGreaterThan(precSqNormX, delta2))
			{
				isInterior = false;
				if (not solveBoundary(precNormR0, delta)) { m_info = Info::BREAKDOWN; return delta; }
			}
		}
		// resume Lanczos iteration
		BasicLinalg::axpbypz(-m_alpha.back(), Bv, -m_beta.back(), Bv_prev, size, m_w.get());
		reOrthonormalize(size, m_w.get());
		// m_w = \beta_{k}Bv_{k+1}
		invB(m_w.get(), m_v.get());
		m_beta.push_back( sqrt(BasicLinalg::inner(m_w.get(), m_v.get(), size)) );
		
		const Scalar invBeta = Scalar(1) / m_beta.back();
		
		BasicLinalg::scal(invBeta, size, m_v.get());
		addBvNext(invBeta, m_w.get(), size);
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
		
		H(m_v.get(), m_w.get());
		// m_w = Hv
		m_alpha.push_back( BasicLinalg::inner(m_v.get(), m_w.get(), size) );
		
		if (not solveBoundary(precNormR0, delta)) { m_info = Info::BREAKDOWN; return delta; }
		
		// resume Lanczos iteration
		BasicLinalg::axpbypz(-m_alpha.back(), getBvCurr(), -m_beta.back(), getBvPrev(), size, m_w.get());
		reOrthonormalize(size, m_w.get());
		
		// m_w = \beta_{k}Bv_{k+1}
		invB(m_w.get(), m_v.get());
		m_beta.push_back( sqrt(BasicLinalg::inner(m_w.get(), m_v.get(), size)) );
		
		const Scalar invBeta = Scalar(1) / m_beta.back();
		BasicLinalg::scal(invBeta, size, m_v.get());
		addBvNext(invBeta, m_w.get(), size);
		
		m_precNormR = abs(m_beta.back()*m_h.back());
	}
	
	std::fill(x, x + size, Scalar{});
	
	m_modelReduction = {};
	
	if constexpr (storeLanczosVectors)
	{
		for (Size i=0; i!=Size(m_h.size()); ++i)
		{
			invB(getBv(i+1), m_v.get());
			H(m_v.get(), m_w.get());
			
			m_modelReduction = BasicLinalg::updateModelReduction(m_modelReduction, x, g, m_h[i], m_v.get(), m_w.get(), size);
			
			BasicLinalg::axpy(m_h[i], m_v.get(), size, x);
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
		
		invB(getBvCurr(), m_v.get());
		
		const std::span<const Scalar> spanH(m_h.data(), size_t(m_h.size() - 1));
		
		for (const auto [h_i, alpha_i, beta_im1, beta_i] : misc::zipLeadBy<0>(spanH, m_alpha, m_beta, m_beta | std::views::drop(1)))
		{
			H(m_v.get(), m_w.get());
			m_modelReduction = BasicLinalg::updateModelReduction(m_modelReduction, x, g, h_i, m_v.get(), m_w.get(), size);
			
			BasicLinalg::axpy(h_i, m_v.get(), size, x);
			
			BasicLinalg::axpbypz(-alpha_i, getBvCurr(), -beta_im1, getBvPrev(), size, m_w.get());
			reOrthonormalize(size, m_w.get());
			
			invB(m_w.get(), m_v.get());
			
			const Scalar invNextBeta = Scalar(1) / beta_i;
			
			BasicLinalg::scal(invNextBeta, size, m_v.get());
			addBvNext(invNextBeta, m_w.get(), size);
		}
		{
			H(m_v.get(), m_w.get());
			m_modelReduction = BasicLinalg::updateModelReduction(m_modelReduction, x, g, m_h.back(), m_v.get(), m_w.get(), size);
		
			BasicLinalg::axpy(m_h.back(), m_v.get(), size, x);
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
