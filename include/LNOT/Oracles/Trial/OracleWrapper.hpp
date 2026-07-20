#ifndef LNOT_ORACLE_WRAPPER_HPP
#define LNOT_ORACLE_WRAPPER_HPP

#include <concepts>

#include <LNOT/Oracles/OracleBase.hpp>
#include <LNOT/Oracles/helpers.hpp>
#include <LNOT/Oracles/concepts.hpp>

namespace LNOT
{

struct VoidFunctor {}; ///<  @brief Marker functor representing absence of a callable object.

template<class Function> concept CVoidFunctor = std::same_as<Function, VoidFunctor>;

template<typename T, typename S, class Function, class Gradient, class HessianOp, class PrecondOp>  class OracleWrapper; 

template<typename T, typename S, CFunction<T> Function, CGradient<T> Gradient, CHessianOp<T> HOp, CHessianOp<T> PrecOp>
struct OracleTraits< OracleWrapper<T, S, Function, Gradient, HOp, PrecOp> > 
{
	using Size      = S; 
	using Scalar    = T; 
	using HessianOp = detail::AsLinearOp<HOp, Scalar>;
	using PrecondOp = detail::AsLinearOp<PrecOp, Scalar>;
};

template<typename T, typename S, CFunction<T> Function, CGradient<T> Gradient, CHessianOp<T> HOp>
struct OracleTraits< OracleWrapper<T, S, Function, Gradient, HOp, VoidFunctor> > 
{
	using Size      = S; 
	using Scalar    = T; 
	using HessianOp = detail::AsLinearOp<HOp, Scalar>;
};

template<typename T, typename S, CFunction<T> Function, CGradient<T> Gradient>
struct OracleTraits< OracleWrapper<T, S, Function, Gradient, VoidFunctor, VoidFunctor> > 
{
	using Size      = S; 
	using Scalar    = T; 
};

/**
 * @brief Basic Oracle that wraps functors that evaluates a functions, value, its gradient and also the product between its Hessian and a given direction.
 * @tparam T Scalar type used in computations.
 * @tparam Function Functor type for function evaluation. Must be callable as `T f(const T*)`.
 * @tparam Gradient Optional gradient functor. Callable as `void grad(const T*, T*)`.
 * @tparam HessianOp Optional Hessian-vector product functor. Callable as `void hessOp(const T*, const T*, T*)`.
 */
template<typename T, typename S, class Function, class Gradient = VoidFunctor, class HOp = VoidFunctor, class PrecOp = VoidFunctor> 
class OracleWrapper : public OracleBase< OracleWrapper<T, S, Function, Gradient, HOp, PrecOp> >
{
	using Self = OracleWrapper<T, S, Function, Gradient, HOp, PrecOp>;
public:
	LNOT_DEFINE_ORACLE

	constexpr OracleWrapper() = delete;
	
	constexpr OracleWrapper(const Size nDims, const Function& func) requires (CVoidFunctor<Gradient> and CVoidFunctor<HOp> and CVoidFunctor<PrecOp>)
		: m_nDims(nDims)
		, m_function(func) {}
		
	constexpr OracleWrapper(const Size nDims, Function&& func) requires (CVoidFunctor<Gradient> and CVoidFunctor<HOp> and CVoidFunctor<PrecOp>) 
		: m_nDims(nDims)
		, m_function(std::move(func)) {}
		
	constexpr OracleWrapper(const Size nDims, const Function& func, const Gradient& grad) requires (hasGradient and CVoidFunctor<HOp> and CVoidFunctor<PrecOp>) 
		: m_nDims(nDims)
		, m_function(func)
		, m_gradient(grad) {}
		
	constexpr OracleWrapper(const Size nDims, Function&& func, Gradient&& grad) requires (hasGradient and CVoidFunctor<HOp> and CVoidFunctor<PrecOp>)
		: m_nDims(nDims)
		, m_function(std::move(func))
		, m_gradient(std::move(grad)) {}
		
	constexpr OracleWrapper(const Size nDims, const Function& func, const Gradient& grad, const HOp& hessOp) requires (hasHessianOp and CVoidFunctor<PrecOp>) 
		: m_nDims(nDims)
		, m_function(func)
		, m_gradient(grad)
		, m_hessianOp(hessOp) {}
		
	constexpr OracleWrapper(const Size nDims, Function&& func, Gradient&& grad, HOp&& hessOp) requires (hasHessianOp and CVoidFunctor<PrecOp>) 
		: m_nDims(nDims)
		, m_function(std::move(func))
		, m_gradient(std::move(grad))
		, m_hessianOp(std::move(hessOp)) {}
		
	constexpr OracleWrapper(const Size nDims, const Function& func, const Gradient& grad, const HOp& hessOp, const PrecOp& precOp) requires (hasHessianOp and hasPrecondOp) 
		: m_nDims(nDims)
		, m_function(func)
		, m_gradient(grad)
		, m_hessianOp(hessOp)
		, m_precondOp(precOp) {}
		
	constexpr OracleWrapper(const Size nDims, Function&& func, Gradient&& grad, HOp&& hessOp, PrecOp&& precOp) requires (hasHessianOp and hasPrecondOp)
		: m_nDims(nDims)
		, m_function(std::move(func))
		, m_gradient(std::move(grad))
		, m_hessianOp(std::move(hessOp))
		, m_precondOp(std::move(precOp)) {}
	
	constexpr Size getNDimsImpl() const { return m_nDims; }
	
	constexpr HessianOp makeHessianOpImpl() requires(not CVoidFunctor<HOp>) { return HessianOp(m_hessianOp); }
	
	constexpr PrecondOp makePrecondOpImpl() requires(not CVoidFunctor<PrecOp>) { return PrecondOp(m_precondOp); }
	
	constexpr Scalar getValue(const Scalar* x) { return m_func(x); }
	
	constexpr void getGradient(const Scalar* x, Scalar* g) { m_grad(x, g); }
	
	constexpr void getHessianOp(const Scalar* x, HessianOp& hessOp) requires(not CVoidFunctor<HOp>) { hessOp.setCurrentPoint(x); }
	
	constexpr void getPrecondOp(const Scalar* x, PrecondOp& precOp) requires(not CVoidFunctor<PrecOp>) { precOp.setCurrentPoint(x); }
private:
	Size m_nDims;
	
	Function  m_function;
	Gradient  m_gradient;
	HOp       m_hessianOp;
	PrecOp    m_precondOp;
};


template<typename T, typename S, class Function>
COracle auto makeOracle(const S size, Function&& function) { return OracleWrapper<T, S, std::remove_cvref_t<Function>>(size, std::forward<Function>(function)); }

template<typename T, typename S, class Function, class Gradient>
CFirstOrderOracle auto makeOracle(const S size, Function&& function, Gradient&& gradient) { return OracleWrapper<T, S, std::remove_cvref_t<Function>, std::remove_cvref_t<Gradient>>(size, std::forward<Function>(function), std::forward<Gradient>(gradient)); }

template<typename T, typename S, class Function, class Gradient, class Hessian>
CSecondOrderOracle auto makeOracle(const S size, Function&& function, Gradient&& gradient, Hessian&& hessian) { return OracleWrapper<T, S, std::remove_cvref_t<Function>, std::remove_cvref_t<Gradient>, std::remove_cvref_t<Hessian>>(size, std::forward<Function>(function), std::forward<Gradient>(gradient), std::forward<Hessian>(hessian)); }

template<typename T, typename S, class Function, class Gradient, class Hessian, class PrecondOp>
CSecondOrderOracle auto makeOracle(const S size, Function&& function, Gradient&& gradient, Hessian&& hessian, PrecondOp&& precond) { return OracleWrapper<T, S, std::remove_cvref_t<Function>, std::remove_cvref_t<Gradient>, std::remove_cvref_t<Hessian>, std::remove_cvref_t<PrecondOp>>(size, std::forward<Function>(function), std::forward<Gradient>(gradient), std::forward<Hessian>(hessian), std::forward<PrecondOp>(precond)); }

} // namespace LNOT

#endif // LNOT_ORACLE_WRAPPER_HPP
