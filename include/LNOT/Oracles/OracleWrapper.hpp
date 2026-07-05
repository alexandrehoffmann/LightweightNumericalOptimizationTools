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

template<class Function> struct IsVoidFunctor : std::bool_constant< std::is_same<Function, VoidFunctor>::value > {}; ///<  @brief Trait to determine whether a type is VoidFunctor.

template<typename T, typename S, class Function, class Gradient, class HessianOp, class PrecondOp>  class OracleWrapper; 


template<typename T, typename S, class Function, class Gradient, class HessianOp, class PrecondOp> 
struct OracleTraits< OracleWrapper<T, S, Function, Gradient, HessianOp, PrecondOp> > 
{
	using Size   = S; 
	using Scalar = T; 

	static_assert(not IsVoidFunctor<Function>::value  and CFunction<Function, Scalar>);
	static_assert(    IsVoidFunctor<Gradient>::value   or CGradient<Gradient, Scalar>);
	static_assert(    IsVoidFunctor<HessianOp>::value  or CHessianOp<HessianOp, Scalar>);
	static_assert(    IsVoidFunctor<PrecondOp>::value  or CHessianOp<PrecondOp, Scalar>);

	static constexpr bool hasGradient     = (not IsVoidFunctor<Gradient>::value);
	static constexpr bool hasHessianProd  = (not IsVoidFunctor<HessianOp>::value);
	static constexpr bool hasApplyPrecond = (not IsVoidFunctor<PrecondOp>::value);
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

	OracleWrapper() = delete;
	
	OracleWrapper(const Size nDims, const Function& func) requires (not hasGradient and not hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(func) {}
		
	OracleWrapper(const Size nDims, Function&& func) requires (not hasGradient and not hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(std::move(func)) {}
		
	OracleWrapper(const Size nDims, const Function& func, const Gradient& grad) requires (hasGradient and not hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(func)
		, m_gradient(grad) {}
		
	OracleWrapper(const Size nDims, Function&& func, Gradient&& grad) requires (hasGradient and not hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(std::move(func))
		, m_gradient(std::move(grad)) {}
		
	OracleWrapper(const Size nDims, const Function& func, const Gradient& grad, const HessianOp& hessOp) requires (hasGradient and hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(func)
		, m_gradient(grad)
		, m_hessianOp(hessOp) {}
		
	OracleWrapper(const Size nDims, Function&& func, Gradient&& grad, HessianOp&& hessOp) requires (hasGradient and hasHessianProd and not hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(std::move(func))
		, m_gradient(std::move(grad))
		, m_hessianOp(std::move(hessOp)) {}
		
	OracleWrapper(const Size nDims, const Function& func, const Gradient& grad, const HessianOp& hessOp, const PrecondOp& precOp) requires (hasGradient and hasHessianProd and hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(func)
		, m_gradient(grad)
		, m_hessianOp(hessOp)
		, m_precondOp(precOp) {}
		
	OracleWrapper(const Size nDims, Function&& func, Gradient&& grad, HessianOp&& hessOp, PrecondOp&& precOp) requires (hasGradient and hasHessianProd and hasApplyPrecond) 
		: m_nDims(nDims)
		, m_function(std::move(func))
		, m_gradient(std::move(grad))
		, m_hessianOp(std::move(hessOp))
		, m_precondOp(std::move(precOp)) {}
	
	Size getNDimsImpl() const { return m_nDims; }

	void setCurrentPointImpl(const Scalar* x) { m_x = x; }
	
	bool isFeasibleImpl() const { return std::isfinite(getValueImpl()); }
	
	Scalar getValueImpl() const { return m_function(m_x); }

	void getGradientImpl(Scalar* g) const requires (hasGradient) { m_gradient(m_x, g); }

	void getHessianProdImpl(const Scalar* d, Scalar* Hd) const requires (hasHessianProd) { m_hessianOp(m_x, d, Hd); }
	
	void applyPrecondImpl(const Scalar* d, Scalar* invBd) const requires (hasApplyPrecond) { m_precondOp(m_x, d, invBd); }
private:
	Size m_nDims;
	const Scalar* m_x = nullptr;
	
	Function  m_function;
	Gradient  m_gradient;
	HessianOp m_hessianOp;
	PrecondOp m_precondOp;
};


template<typename T, typename S, class Function>
auto  makeOracle(const S size, Function&& function) { return OracleWrapper<T, S, std::decay_t<Function>>(size, std::forward<Function>(function)); }

template<typename T, typename S, class Function, class Gradient>
auto makeOracle(const S size, Function&& function, Gradient&& gradient) { return OracleWrapper<T, S, std::decay_t<Function>, std::decay_t<Gradient>>(size, std::forward<Function>(function), std::forward<Gradient>(gradient)); }

template<typename T, typename S, class Function, class Gradient, class Hessian>
auto makeOracle(const S size, Function&& function, Gradient&& gradient, Hessian&& hessian) { return OracleWrapper<T, S, std::decay_t<Function>, std::decay_t<Gradient>, std::decay_t<Hessian>>(size, std::forward<Function>(function), std::forward<Gradient>(gradient), std::forward<Hessian>(hessian)); }

template<typename T, typename S, class Function, class Gradient, class Hessian, class PrecondOp>
auto makeOracle(const S size, Function&& function, Gradient&& gradient, Hessian&& hessian, PrecondOp&& precond) { return OracleWrapper<T, S, std::decay_t<Function>, std::decay_t<Gradient>, std::decay_t<Hessian>, std::decay_t<PrecondOp>>(size, std::forward<Function>(function), std::forward<Gradient>(gradient), std::forward<Hessian>(hessian), std::forward<PrecondOp>(precond)); }

} // namespace LNOT

#endif // LNOT_ORACLE_WRAPPER_HPP
