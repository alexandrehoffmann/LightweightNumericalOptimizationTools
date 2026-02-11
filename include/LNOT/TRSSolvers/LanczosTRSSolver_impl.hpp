#ifndef LNOT_LANCZOS_TRS_SOLVER_IMPL_HPP
#define LNOT_LANCZOS_TRS_SOLVER_IMPL_HPP

#include <LNOT/TRSSolvers/LanczosTRSSolver.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <LNOT/BasicLinalg/DiagonalPreconditionerOp.hpp>
#include <LNOT/FloatingPoint/FPComparator.hpp>

#include <algorithm>

#include <fmt/core.h>
#include <fmt/format.h>

namespace LNOT
{

//// explicit template instanciations ////

extern template class LanczosTRSSolver<float>;
extern template class LanczosTRSSolver<double>;

// explicit instanciation for float
extern template void LanczosTRSSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

extern template void LanczosTRSSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<float>::solveImpl(const SymmetricDenseMatrixOp<float, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<float>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

// explicit instanciation for double
extern template void LanczosTRSSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const IdentityPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

extern template void LanczosTRSSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::ROW_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::LOWER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);
extern template void LanczosTRSSolver<double>::solveImpl(const SymmetricDenseMatrixOp<double, StorageOrder::COL_MAJOR, UpLo::UPPER>&  H, const DiagonalPreconditionerOp<double>& invB, const Scalar* g, const Size size, const Scalar& delta, Scalar* x);

//// method implementations ////

template<typename T>
LanczosTRSSolver<T>::LanczosTRSSolver(const Size maxIt, const Scalar tol, const Size maxItTr, const Scalar tolTr) 
	: Base(maxIt, tol, tolTr)
	, m_maxItTr(maxItTr) 
{
}

template<typename T>
void LanczosTRSSolver<T>::clearWorkSpace()
{	
	if (m_Bv_old != nullptr) { delete[] m_Bv_old; m_Bv_old = nullptr; }
	if (m_Bv     != nullptr) { delete[] m_Bv;     m_Bv     = nullptr; }
	if (m_v      != nullptr) { delete[] m_v;      m_v      = nullptr; }
	if (m_p      != nullptr) { delete[] m_p;      m_p      = nullptr; }
	if (m_Bp     != nullptr) { delete[] m_Bp;     m_Bp     = nullptr; }
	if (m_Hp     != nullptr) { delete[] m_Hp;     m_Hp     = nullptr; }
	if (m_w      != nullptr) { delete[] m_w;      m_w      = nullptr; }
	
	m_workCapacity = 0;
}

template<typename T>
void LanczosTRSSolver<T>::resizeWorkSpace(const Size newSize)
{
	if (m_workCapacity < newSize)
	{
		clearWorkSpace();
		m_workCapacity = newSize;
		m_Bv_old  = new Scalar[m_workCapacity];
		m_Bv      = new Scalar[m_workCapacity];
		m_v       = new Scalar[m_workCapacity];
		m_p       = new Scalar[m_workCapacity];
		m_Bp      = new Scalar[m_workCapacity];
		m_Hp      = new Scalar[m_workCapacity];
		m_w       = new Scalar[m_workCapacity];
	}
}

template<typename T> template<typename HesOp, typename PrecOp, typename ASize> 
void LanczosTRSSolver<T>::solveImpl(const HesOp& H, const PrecOp& invB, const Scalar* g, const ASize size, const Scalar& delta, Scalar* x) requires(AreHessianOps<HesOp,PrecOp>::value and IsSize<ASize>::value)
{	
	using std::abs;
	
	resizeWorkSpace(size);
	
	m_alpha.clear();
	m_beta.clear();
	
	std::fill(x, x + size, 0);
	Scalar precSqNormX = 0;
	
	m_lambda = 0;
	m_modelReduction = 0;
	
	bool isInterior = true;
	std::fill(m_Bv_old, m_Bv_old + size, 0);
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { m_Bv[i] = -g[i]; }
	
	invB(m_Bv, m_v);
	
	const Scalar precNormR0 = std::sqrt(BasicLinalg::inner(m_Bv, m_v, size)); 
	
	m_beta.push_back(0);
	
	Scalar l_old = 0;
	Scalar eta   = -precNormR0;
	
	BasicLinalg::scal(Scalar(1) / eta, size, m_Bv);
	BasicLinalg::scal(Scalar(1) / eta, size, m_v);
	
	std::fill(m_p,  m_p  + size, 0);
	std::fill(m_Bp, m_Bp + size, 0);
	std::fill(m_Hp, m_Hp + size, 0);
	
	m_precNormR = precNormR0;
	
	const Scalar relTol = m_relTol*m_precNormR;
	const Scalar delta2 = delta*delta;
	
	const FPComparator<Scalar> cmp;
	const FPComparator<Scalar> cmpTr(m_relTolTr, m_absTolTr);
	
	m_info = Info::FAILURE;
	if (m_out) { fmt::print(m_out, "#Preconditioned Lanczos TRS solver : \n#Iteration residual lambda relative_tol absolute_tol\n"); }
	
	// first try to solve within the trust region
	for (m_nIt=0; m_nIt!=m_maxIt and isInterior; ++m_nIt)
	{
		if (m_out) { fmt::print(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}\n", m_nIt, m_precNormR, m_lambda, relTol, m_absTol); }
		if (m_precNormR < relTol or m_precNormR < m_absTol) { m_info = Info::SUCCESS; return; }
		
		H(m_v, m_w);
		// m_w = Hv_{k}
		m_alpha.push_back( BasicLinalg::inner(m_v, m_w, size) );
		// solve T_k h_k = -\|r_0\|e_1
		const Scalar d = m_alpha.back() - m_beta.back()*l_old;
		const Scalar invD = Scalar(1) / d;
		
		if (not cmp.isDefPositive(d)) 
		{
			isInterior = false; 
			if (not solveBoundary(precNormR0, delta)) { m_info = Info::BREAKDOWN; return; }
		}
		else
		{
			#pragma omp simd
			for (Size i=0; i!=size; ++i) 
			{ 
				m_p[i]  = invD*(m_v[i]  - m_beta.back()*m_p[i]); 
				m_Bp[i] = invD*(m_Bv[i] - m_beta.back()*m_Bp[i]); 
				m_Hp[i] = invD*(m_w[i]  - m_beta.back()*m_Hp[i]); 
			} 
			for (Size i=0; i!=size; ++i) { m_modelReduction += eta*(x[i]*m_Hp[i] + Scalar(0.5)*eta*m_p[i]*m_Hp[i] + m_p[i]*g[i]); }
			
			precSqNormX += 2*eta*BasicLinalg::inner(x, m_Bp, size) + eta*eta*BasicLinalg::inner(m_p, m_Bp, size);
			BasicLinalg::axpy(eta, m_p, size, x);
			
			if (cmpTr.isDefGreaterThan(precSqNormX, delta2))
			{
				isInterior = false;
				if (not solveBoundary(precNormR0, delta)) { m_info = Info::BREAKDOWN; return; }
			}
		}
		// resume Lanczos iteration
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_w[i] += -m_alpha.back()*m_Bv[i] - m_beta.back()*m_Bv_old[i]; }
		// m_w = \beta_{k}Bv_{k+1}
		invB(m_w, m_v);
		m_beta.push_back( std::sqrt(BasicLinalg::inner(m_w, m_v, size)) );
		
		const Scalar invBeta = Scalar(1) / m_beta.back();
		BasicLinalg::scal(invBeta, size, m_v);
		std::copy(m_Bv, m_Bv + size, m_Bv_old);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Bv[i] = invBeta*m_w[i]; }
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
		if (m_out) { fmt::print(m_out, "{} {:10.2e} {:10.2e} {:10.2e} {:10.2e}\n", m_nIt, m_precNormR, m_lambda, relTol, m_absTol); }
		if (m_precNormR < relTol or m_precNormR < m_absTol) { m_info = Info::SUCCESS; break; }
		
		H(m_v, m_w);
		// m_w = Hv
		m_alpha.push_back( BasicLinalg::inner(m_v, m_w, size) );
		
		if (not solveBoundary(precNormR0, delta)) { m_info = Info::BREAKDOWN; return; }
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_w[i] += -m_alpha.back()*m_Bv[i] - m_beta.back()*m_Bv_old[i]; }
		// m_w = \beta_{k}Bv_{k+1}
		invB(m_w, m_v);
		m_beta.push_back( std::sqrt(BasicLinalg::inner(m_w, m_v, size)) );
		
		const Scalar invBeta = Scalar(1) / m_beta.back();
		BasicLinalg::scal(invBeta, size, m_v);
		std::copy(m_Bv, m_Bv + size, m_Bv_old);
		#pragma omp simd
		for (Size i=0; i!=size; ++i) { m_Bv[i] = invBeta*m_w[i]; }
		
		m_precNormR = abs(m_beta.back()*m_h.back());
	}
	// re-run Lanczos iteration to compute x = Vh;
	using Iterator = typename std::vector<Scalar>::const_iterator;
	
	std::fill(x, x + size, 0);
	std::fill(m_Bv_old, m_Bv_old + size, 0);
	m_modelReduction = 0;
	#pragma omp simd
	for (Size i=0; i!=size; ++i) { m_Bv[i] = g[i] / precNormR0; } 
	
	invB(m_Bv, m_v);
	
	Iterator it_alpha = m_alpha.cbegin();
	Iterator it_beta  = m_beta.cbegin();
	for (Iterator it_h = m_h.cbegin(); it_h!=m_h.cend()-1; ++it_h, ++it_alpha, ++it_beta, ++m_nIt)
	{
		H(m_v, m_w);
		for (Size i=0; i!=size; ++i) { m_modelReduction += (*it_h)*(x[i]*m_w[i] + Scalar(0.5)*(*it_h)*m_v[i]*m_w[i] + m_v[i]*g[i]); }
		
		BasicLinalg::axpy(*it_h, m_v, size, x);
		const Scalar invNextBeta = Scalar(1) / *std::next(it_beta);
		
		#pragma omp simd
		for (Size i=0; i!=size; ++i)
		{
			const Scalar new_Bv_i = invNextBeta*(m_w[i] - (*it_alpha)*m_Bv[i] - (*it_beta)*m_Bv_old[i]);
			m_Bv_old[i] = m_Bv[i];
			m_Bv[i] = new_Bv_i;
		}
		invB(m_Bv, m_v);
	}
	{
		H(m_v, m_w);
		for (Size i=0; i!=size; ++i) { m_modelReduction += m_h.back()*(x[i]*m_w[i] + Scalar(0.5)*m_h.back()*m_v[i]*m_w[i] + m_v[i]*g[i]); }
		
		BasicLinalg::axpy(m_h.back(), m_v, size, x);
		
		++m_nIt;
	}	
}

template<typename T>
bool LanczosTRSSolver<T>::solveBoundary(const Scalar& gamma, const Scalar& delta)
{	
	using std::abs;
	using std::sqrt;
	
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
	Scalar lambdaMin = std::max(Scalar(0), abs(gamma) / delta - norm1_T);
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

#endif // LNOT_LANCZOS_TRS_SOLVER_IMPL_HPP
