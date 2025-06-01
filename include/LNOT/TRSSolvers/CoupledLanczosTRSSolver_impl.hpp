#ifndef LNOT_COUPLED_LANCZOS_TRS_SOLVER_IMPL_HPP
#define LNOT_COUPLED_LANCZOS_TRS_SOLVER_IMPL_HPP

#include <span>

namespace LNOT
{
	
template<typename T>
CoupledLanczosTRSSolver<T>::CoupledLanczosTRSSolver(const Size maxIt, const Scalar tol, const Size maxItTr, const Scalar tolTr) 
	: Base(maxIt, tol, tolTr)
	, m_maxItTr(maxItTr) 
{
}

template<typename T>
void CoupledLanczosTRSSolver<T>::clearWorkSpace()
{
	if (m_v     != nullptr) { delete[] m_v;     }
	if (m_q     != nullptr) { delete[] m_q;     }
	if (m_w     != nullptr) { delete[] m_w;     }
	if (m_w_old != nullptr) { delete[] m_w_old; }
	Base::m_workCapacity = 0;
}

template<typename T> template<typename Op> 
void CoupledLanczosTRSSolver<T>::solve(const Op& H, const Scalar* g, const Size size, const Scalar& delta, Scalar* x) requires (IsHessianOp<Op>::value)
{
	if (Base::m_workCapacity < size)
	{
		clearWorkSpace();
		Base::m_workCapacity = size;
		m_v     = new Scalar[Base::m_workCapacity];
		m_q     = new Scalar[Base::m_workCapacity];
		m_w     = new Scalar[Base::m_workCapacity];
		m_w_old = new Scalar[Base::m_workCapacity];
	}
	m_alpha.clear();
	m_beta.clear();
	
	std::fill(x, x+size, 0);
	m_lambda = 0;
	Base::m_modelReduction = 0;
	bool isInterior = true;
	
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { m_v[i] = -g[i]; } 
	
	const Scalar normR0 = BasicLinalg::norm(m_v, size);
	
	m_beta.push_back(0);
	Scalar l_old = 0;
	Scalar eta = -normR0;
	Scalar tau_old = 0;
	
	BasicLinalg::scal(1. / eta, size, m_v);
	std::copy(m_v, m_v + size, m_q);
	
	m_normR = normR0;
	
	const Scalar tol = Base::m_tol*m_normR;
	const Scalar deltaTol2 = (delta + Base::m_tolTr)*(delta + Base::m_tolTr);
	
	Base::m_info = Info::FAILURE;
	if (Base::m_out) { fmt::print(Base::m_out, "Coupled recurrence Lanczos TRS solver : \n#Iteration residual lambda tol\n"); }
	for (Base::m_nIt=0; Base::m_nIt!=Base::m_maxIt and isInterior; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e} {:10.2e}\n", Base::m_nIt, m_normR, m_lambda, tol); }
		if (m_normR < tol) { Base::m_info = Info::SUCCESS; return; }
		H(m_q, m_w);
		
		const Scalar tau = BasicLinalg::inner(m_q, m_w, size);
		m_alpha.push_back(tau + tau_old*l_old*l_old);
		// solve T_k h_k = -\|r_0\|e_1
		const Scalar d = m_alpha.back() - m_beta.back()*l_old;
		const Scalar invD = 1. / d;
		if (d < std::numeric_limits<Scalar>::epsilon()) 
		{ 
			isInterior = false; 
			if (not solveBoundary(normR0, delta)) { Base::m_info = Info::BREAKDOWN; return; }
		}
		else
		{
			#pragma omp simd reduction(+:Base::m_modelReduction)
			for (Size i=0; i!=size; ++i) { Base::m_modelReduction += invD*eta*(x[i]*m_w[i] + 0.5*invD*eta*m_q[i]*m_w[i] + m_q[i]*g[i]); }
			BasicLinalg::axpy(invD*eta, m_q, size, x);
			if (BasicLinalg::squaredNorm(x, size) > deltaTol2)
			{
				isInterior = false;
				if (not solveBoundary(normR0, delta)) { Base::m_info = Info::BREAKDOWN; return; }
			}
		}
		// Resume Lanczos iteration
		BasicLinalg::axpy(-tau, m_v, size, m_w);
		m_beta.push_back( BasicLinalg::norm(m_w, size) );
		// prepare next solve
		const Scalar l = invD*m_beta.back();
		
		eta     *= -l;
		m_normR  = isInterior ? std::abs(eta) : std::abs(m_beta.back()*m_h.back()); 
		l_old    = l;
		tau_old  = tau;
		const Scalar invBeta = 1. / m_beta.back();
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_v[i] = invBeta*m_w[i]; m_q[i] = m_v[i] - l*m_q[i]; } 
	}
	for (; Base::m_nIt!=Base::m_maxIt; ++Base::m_nIt)
	{
		if (Base::m_out) { fmt::print(Base::m_out, "{} {:10.2e} {:10.2e} {:10.2e}\n", Base::m_nIt, m_normR, m_lambda, tol); }
		if (m_normR < tol) { Base::m_info = Info::SUCCESS; break; }
		H(m_q, m_w);
		
		const Scalar tau = BasicLinalg::inner(m_q, m_w, size);
		m_alpha.push_back( tau + tau_old*l_old*l_old );
		const Scalar d = m_alpha.back() - m_beta.back()*l_old;
		
		if (not solveBoundary(normR0, delta)) { Base::m_info = Info::BREAKDOWN; return; }
		
		// Resume Lanczos iteration
		BasicLinalg::axpy(-tau, m_v, size, m_w);
		m_beta.push_back( BasicLinalg::norm(m_w, size) );
		// prepare next solve
		const Scalar l = m_beta.back() / d;
		
		l_old    = l;
		tau_old  = tau;
		m_normR  = std::abs(m_beta.back()*m_h.back());
		const Scalar invBeta = 1. / m_beta.back();
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_v[i] = invBeta*m_w[i]; m_q[i] = m_v[i] - l*m_q[i]; } 
	}
	// re-run Lanczos iteration to compute x = Vh;
	using Iterator = typename std::vector<Scalar>::const_iterator;
	 
	std::fill(x, x + size, 0);
	Base::m_modelReduction = 0;
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { m_v[i] = g[i] / normR0; } 
	std::copy(m_v, m_v + size, m_q);
	
	std::fill(m_w_old, m_w_old + size, 0);
	
	l_old = 0;
	
	Iterator it_alpha = m_alpha.cbegin();
	Iterator it_beta  = m_beta.cbegin();
	for (Iterator it_h = m_h.cbegin(); it_h!=m_h.cend()-1; ++it_h, ++it_alpha, ++it_beta, ++Base::m_nIt)
	{		
		H(m_q, m_w);

		#pragma omp simd reduction(+:Base::m_modelReduction)
		for (Size i=0; i!=size; ++i) 
		{ 
			const Scalar Hv_i = m_w[i] + l_old*m_w_old[i];
			Base::m_modelReduction += (*it_h)*(x[i]*Hv_i + 0.5*(*it_h)*m_v[i]*Hv_i + m_v[i]*g[i]); 
		}
		
		BasicLinalg::axpy(*it_h, m_v, size, x);
		
		const Scalar tau = BasicLinalg::inner(m_q, m_w, size);
		const Scalar d = (*it_alpha) - (*it_beta)*l_old;
		const Scalar l = (*std::next(it_beta)) / d;
		l_old = l;
		const Scalar invNextBeta = 1. / *std::next(it_beta);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_v[i] = invNextBeta*(m_w[i] - tau*m_v[i]); m_q[i] = m_v[i] - l*m_q[i]; }
		std::copy(m_w, m_w + size, m_w_old);
	}
	{		
		H(m_q, m_w);
		
		#pragma omp simd reduction(+:Base::m_modelReduction)
		for (Size i=0; i!=size; ++i) 
		{ 
			const Scalar Hv_i = m_w[i] + l_old*m_w_old[i];
			Base::m_modelReduction += m_h.back()*(x[i]*Hv_i + 0.5*m_h.back()*m_v[i]*Hv_i + m_v[i]*g[i]); 
		}
		
		BasicLinalg::axpy(m_h.back(), m_v, size, x);
	}
}

template<typename T>
bool CoupledLanczosTRSSolver<T>::solveBoundary(const Scalar& __restrict__ gamma, const Scalar& __restrict__ delta)
{	
	namespace TridiagLDLt = BasicLinalg::Tridiag::LDLt;
	
	const Size   size      = Size(m_alpha.size());
	const Scalar delta2    = delta*delta;
	const Scalar tol_delta = std::max(delta*Base::m_tolTr, std::numeric_limits<Scalar>::epsilon());
	const Scalar norm1_T   = BasicLinalg::Tridiag::norm1(m_alpha.data(), std::next(m_beta.data()), size);

	m_h.reserve(m_alpha.capacity());
	m_invD.reserve(m_alpha.capacity());
	m_l.reserve(m_alpha.capacity()-1);
	
	m_h.resize(size);
	m_invD.resize(size);
	m_l.resize(size-1);
	
	Scalar lambdaMin = 0;
	Scalar lambdaMax = std::abs(gamma) / delta + norm1_T;
	m_lambda = std::max(lambdaMin, std::min(m_lambda, lambdaMax));
		
	for (size_t it=0;it!=m_maxItTr;++it)
	{	
		const bool ldltSuccess = TridiagLDLt::compute(m_alpha.data(), std::next(m_beta.data()), size, m_lambda, m_invD.data(), m_l.data());
		const bool isSpd       = std::ranges::all_of(m_invD, [](const Scalar& invDi) -> bool { return invDi > std::numeric_limits<Scalar>::epsilon(); });
		const bool cholSuccess = ldltSuccess and isSpd;
		if (cholSuccess)
		{
			TridiagLDLt::solve_e1(m_invD.data(), m_l.data(), size, -gamma, m_h.data());
			const Scalar sqNormH = BasicLinalg::squaredNorm(m_h.data(), size);
			const Scalar   normH = std::sqrt(sqNormH);
			
			if (std::abs(normH - delta) < tol_delta) { return true; }
			if (lambdaMax - lambdaMin < std::numeric_limits<Scalar>::epsilon()*std::max(lambdaMax, lambdaMin)) { return false; } // empty interval
			
			TridiagLDLt::solve_L_inplace(m_l.data(), size, m_h.data());
			
			const Scalar sqnorm_w    = BasicLinalg::weightedSquaredNorm(m_h.data(), m_invD.data(), size);
			const Scalar deltaLambda = ((normH - delta) / delta)*(sqNormH / sqnorm_w);
			
			if (sqNormH < delta2) { lambdaMax = m_lambda; }
			else                  { lambdaMin = m_lambda; }
			
			const Scalar gap = (lambdaMax - lambdaMin)*0.005;
			m_lambda = std::max(lambdaMin + gap, std::min(m_lambda + deltaLambda, lambdaMax - gap));
		}
		else
		{
			lambdaMin = m_lambda;
			m_lambda  = 0.5*(lambdaMax + lambdaMin); 
		}
	}
	// hard case cannot occur, we did not converge due to numerical errors
	return false;
}
	
} // namespace LNOT

#endif // LNOT_COUPLED_LANCZOS_TRS_SOLVER_IMPL_HPP
