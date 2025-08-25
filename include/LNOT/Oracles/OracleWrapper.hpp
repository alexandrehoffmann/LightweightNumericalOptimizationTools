#ifndef LNOT_ORACLE_WRAPPER_HPP
#define LNOT_ORACLE_WRAPPER_HPP

#include <algorithm>   // for std::copy
#include <type_traits> // for std::is_invocable
#include <cmath>       // for std::isfinite

#include <LNOT/Oracles/OracleBase.hpp>

namespace LNOT
{

struct VoidFunctor {}; ///<  @brief Marker functor representing absence of a callable object.

template<class Function> struct IsVoidFunctor : std::bool_constant< std::is_same<Function, VoidFunctor>::value > {}; ///<  @brief Trait to determine whether a type is VoidFunctor.

template<typename T, class Function, class Gradient, class HessianOp, class PrecondOp>  class OracleWrapper; 


template<typename T, class Function, class Gradient, class HessianOp, class PrecondOp> 
struct OracleTraits< OracleWrapper<T, Function, Gradient, HessianOp, PrecondOp> > 
{
	using Size   = unsigned int; 
	using Scalar = T; 

	static_assert(not IsVoidFunctor<Function>::value  and std::is_invocable<Function, const Scalar*>::value and std::is_same<std::invoke_result_t<Function, Scalar*>, Scalar>::value);
	static_assert(    IsVoidFunctor<Gradient>::value   or std::is_invocable<Gradient, const Scalar*, Scalar*>::value);
	static_assert(    IsVoidFunctor<HessianOp>::value  or std::is_invocable<HessianOp, const Scalar*, const Scalar*, Scalar*>::value);
	static_assert(    IsVoidFunctor<PrecondOp>::value  or std::is_invocable<PrecondOp, const Scalar*, const Scalar*, Scalar*>::value);

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
template<typename T, class Function, class Gradient = VoidFunctor, class HessianOp = VoidFunctor, class PrecondOp=VoidFunctor> 
class OracleWrapper : public OracleBase< OracleWrapper<T, Function, Gradient, HessianOp, PrecondOp> >
{
	using Base = OracleBase< OracleWrapper<T, Function, Gradient, HessianOp, PrecondOp> >;
public:
	using Size   = typename Base::Size;
	using Scalar = typename Base::Scalar;
	
	static constexpr bool hasGradient     = Base::hasGradient;
	static constexpr bool hasHessianProd  = Base::hasHessianProd;
	static constexpr bool hasApplyPrecond = Base::hasApplyPrecond;

	OracleWrapper() = delete;
	OracleWrapper(const Size nDims, Function func)                                                    requires (not hasGradient and not hasHessianProd and not hasApplyPrecond) : m_nDims(nDims), m_function(func) {}
	OracleWrapper(const Size nDims, Function func, Gradient grad)                                     requires (    hasGradient and not hasHessianProd and not hasApplyPrecond) : m_nDims(nDims), m_function(func), m_gradient(grad) {}
	OracleWrapper(const Size nDims, Function func, Gradient grad, HessianOp hessOp)                   requires (    hasGradient and     hasHessianProd and not hasApplyPrecond) : m_nDims(nDims), m_function(func), m_gradient(grad), m_hessianOp(hessOp) {}
	OracleWrapper(const Size nDims, Function func, Gradient grad, HessianOp hessOp, PrecondOp precOp) requires (    hasGradient and     hasHessianProd and     hasApplyPrecond) : m_nDims(nDims), m_function(func), m_gradient(grad), m_hessianOp(hessOp), m_precondOp(precOp) {}
	
	Size getNDims() const { return m_nDims; }

	void setCurrentPoint(const Scalar* x) { m_x = x; }
	
	bool isFeasible() const { return std::isfinite(getValue()); }
	
	Scalar getValue() const { return m_function(m_x); }

	void getGradient(Scalar* g) const requires (hasGradient) { m_gradient(m_x, g); }

	void getHessianProd(const Scalar* d, Scalar* Hd) const requires (hasHessianProd) { m_hessianOp(m_x, d, Hd); }
	
	void applyPrecond(const Scalar* d, Scalar* invBd) const requires (hasApplyPrecond) { m_precondOp(m_x, d, invBd); }
private:
	Size    m_nDims;
	const Scalar* m_x = nullptr;
	
	Function  m_function;
	Gradient  m_gradient;
	HessianOp m_hessianOp;
	PrecondOp m_precondOp;
};


} // namespace LNOT

#endif // LNOT_ORACLE_WRAPPER_HPP
