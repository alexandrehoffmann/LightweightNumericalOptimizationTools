#ifndef LNOT_ORACLE_HELPERS_HPP
#define LNOT_ORACLE_HELPERS_HPP

#include <LNOT/BasicLinalg/IdentityOp.hpp>
#include <LNOT/BasicLinalg/NoOp.hpp>

namespace LNOT::detail
{
	template<class OracleTrait>
	struct HessianOpForTrait
	{
		using Type = NoOp;
	};

	template<class OracleTrait> requires requires { typename OracleTrait::HessianOp; }
	struct HessianOpForTrait<OracleTrait>
	{
		using Type = typename OracleTrait::HessianOp;
	};

	template<class OracleTrait>
	struct PrecondOpForTrait
	{
		using Type = IdentityOp;
	};

	template<class OracleTrait> requires requires { typename OracleTrait::PrecondOp; }
	struct PrecondOpForTrait<OracleTrait>
	{
		using Type = typename OracleTrait::PrecondOp;
	};
	
	template<class Op, typename Scalar>
	class AsLinearOp
	{
	public:	
		constexpr AsLinearOp(Op& op) : m_op(op) {}
		
		constexpr void setCurrentPoint(const Scalar* x) { m_x = x; }
		
		constexpr void eval(const Scalar* d, Scalar* Opd) const { m_op.eval(m_x, d, Opd); }
	private:
		Op&           m_op;
		const Scalar* m_x = nullptr;
	};

} // namespace LNOT::detail

#endif // LNOT_ORACLE_HELPERS_HPP
