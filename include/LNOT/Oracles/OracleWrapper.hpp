#ifndef LNOT_ORACLE_WRAPPER_HPP
#define LNOT_ORACLE_WRAPPER_HPP

#include <algorithm>   // for std::copy
#include <type_traits> // for std::is_invocable
#include <cmath>       // for std::isfinite

#include <LNOT/Oracles/OracleBase.hpp>
#include <LNOT/Oracles/concepts.hpp>

namespace LNOT
{

struct VoidFunctor {}; ///<  @brief Marker functor representing absence of a callable object.

template<class Function> concept CVoidFunctor = std::same_as<Function, VoidFunctor>;

template<typename T, typename S, class Function, class Gradient, class HessianOp, class PrecondOp>  class OracleWrapper; 


template<typename T, typename S, class Function, class Gradient, class HessianOp, class PrecondOp> 
struct OracleTraits< OracleWrapper<T, S, Function, Gradient, HessianOp, PrecondOp> > 
{
	using Size   = S; 
	using Scalar = T; 

	static_assert(not CVoidFunctor<Function>  and CFunction<Function,   Scalar>);
	static_assert(    CVoidFunctor<Gradient>   or CGradient<Gradient,   Scalar>);
	static_assert(    CVoidFunctor<HessianOp>  or CHessianOp<HessianOp, Scalar>);
	static_assert(    CVoidFunctor<PrecondOp>  or CHessianOp<PrecondOp, Scalar>);
};

/**
 * @brief Basic Oracle that wraps functors that evaluates a functions, value, its gradient and also the product between its Hessian and a given direction.
 * @tparam T Scalar type used in computations.
 * @tparam Function Functor type for function evaluation. Must be callable as `T f(const T*)`.
 * @tparam Gradient Optional gradient functor. Callable as `void grad(const T*, T*)`.
 * @tparam HessianOp Optional Hessian-vector product functor. Callable as `void hessOp(const T*, const T*, T*)`.
 */
template<typename T, typename S, class Function, class Gradient = VoidFunctor, class HessianOp = VoidFunctor, class PrecondOp=VoidFunctor> 
class OracleWrapper : public OracleBase< OracleWrapper<T, S, Function, Gradient, HessianOp, PrecondOp> >
{
	using Self = OracleWrapper<T, S, Function, Gradient, HessianOp, PrecondOp>;
public:
	LNOT_DEFINE_ORACLE

	static constexpr bool hasGradient     = CGradient<Gradient,   Scalar>;
	static constexpr bool hasHessianProd  = CHessianOp<HessianOp, Scalar>;
	static constexpr bool hasApplyPrecond = CHessianOp<PrecondOp, Scalar>;

	constexpr OracleWrapper() = delete;
	
	constexpr OracleWrapper(const Size nDims, const Function& func) requires (not hasGradient and not hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(func) {}
		
	constexpr OracleWrapper(const Size nDims, Function&& func) requires (not hasGradient and not hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(std::move(func)) {}
		
	constexpr OracleWrapper(const Size nDims, const Function& func, const Gradient& grad) requires (hasGradient and not hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(func)
		, m_gradient(grad) {}
		
	constexpr OracleWrapper(const Size nDims, Function&& func, Gradient&& grad) requires (hasGradient and not hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(std::move(func))
		, m_gradient(std::move(grad)) {}
		
	constexpr OracleWrapper(const Size nDims, const Function& func, const Gradient& grad, const HessianOp& hessOp) requires (hasGradient and hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(func)
		, m_gradient(grad)
		, m_hessianOp(hessOp) {}
		
	constexpr OracleWrapper(const Size nDims, Function&& func, Gradient&& grad, HessianOp&& hessOp) requires (hasGradient and hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(std::move(func))
		, m_gradient(std::move(grad))
		, m_hessianOp(std::move(hessOp)) {}
		
	constexpr OracleWrapper(const Size nDims, const Function& func, const Gradient& grad, const HessianOp& hessOp, const PrecondOp& precOp) requires (hasGradient and hasHessianProd and hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(func)
		, m_gradient(grad)
		, m_hessianOp(hessOp)
		, m_precondOp(precOp) {}
		
	constexpr OracleWrapper(const Size nDims, Function&& func, Gradient&& grad, HessianOp&& hessOp, PrecondOp&& precOp) requires (hasGradient and hasHessianProd and hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(std::move(func))
		, m_gradient(std::move(grad))
		, m_hessianOp(std::move(hessOp))
		, m_precondOp(std::move(precOp)) {}
	
	constexpr Size getNDimsImpl() const { return m_nDims; }

	constexpr void setCurrentPointImpl(const Scalar* x) { m_x = x; }
	
	constexpr bool isFeasibleImpl() const { return std::isfinite(getValueImpl()); }
	
	constexpr Scalar getValueImpl() const { return m_function(m_x); }

	constexpr void getGradientImpl(Scalar* g) const requires (hasGradient) { m_gradient(m_x, g); }

	constexpr void getHessianProdImpl(const Scalar* d, Scalar* Hd) const requires(hasHessianProd) { m_hessianOp(m_x, d, Hd); }
	
	constexpr void initHessianProd() requires (hasHessianProd) {}
	
	constexpr void applyPrecondImpl(const Scalar* d, Scalar* invBd) const requires(hasApplyPrecond) { m_precondOp(m_x, d, invBd); }
	
	constexpr void initApplyPrecond() requires(hasApplyPrecond) {}
private:
	Size m_nDims;
	const Scalar* m_x = nullptr;
	
	Function  m_function;
	Gradient  m_gradient;
	HessianOp m_hessianOp;
	PrecondOp m_precondOp;
};

template<typename T, typename S, class Function, class Gradient = VoidFunctor, class HessianOp = VoidFunctor, class PrecondOp=VoidFunctor>
using OracleWrapperFor = OracleWrapper<T, S, std::remove_cvref_t<Function>, std::remove_cvref_t<Gradient>, std::remove_cvref_t<HessianOp>, std::remove_cvref_t<PrecondOp>>; 

template<typename T, typename S, class Function>
constexpr COracle auto makeOracle(const S size, Function&& function) { return OracleWrapperFor<T, S, Function>(size, std::forward<Function>(function)); }

template<typename T, typename S, class Function, class Gradient>
constexpr CFirstOrderOracle auto makeOracle(const S size, Function&& function, Gradient&& gradient) { return OracleWrapperFor<T, S, Function, Gradient>(size, std::forward<Function>(function), std::forward<Gradient>(gradient)); }

template<typename T, typename S, class Function, class Gradient, class Hessian>
constexpr CSecondOrderOracle auto makeOracle(const S size, Function&& function, Gradient&& gradient, Hessian&& hessian) { return OracleWrapperFor<T, S, Function, Gradient, Hessian>(size, std::forward<Function>(function), std::forward<Gradient>(gradient), std::forward<Hessian>(hessian)); }

template<typename T, typename S, class Function, class Gradient, class Hessian, class PrecondOp>
constexpr CSecondOrderOracle auto makeOracle(const S size, Function&& function, Gradient&& gradient, Hessian&& hessian, PrecondOp&& precond) { return OracleWrapperFor<T, S, Function, Gradient, Hessian, PrecondOp>(size, std::forward<Function>(function), std::forward<Gradient>(gradient), std::forward<Hessian>(hessian), std::forward<PrecondOp>(precond)); }

} // namespace LNOT

#endif // LNOT_ORACLE_WRAPPER_HPP
