#ifndef LNOT_WEIGHTED_L2_NORM_CONVERGENCE_CRITERION_HPP
#define LNOT_WEIGHTED_L2_NORM_CONVERGENCE_CRITERION_HPP

#include <LNOT/ConvergenceCriterions/CConvergenceCriterionFor.hpp>

#include <LNOT/Assert.hpp>
#include <LNOT/BasicLinalg/BasicLinalg.hpp>

#include <span>

namespace LNOT
{

class WeightedL2NormConvergenceCriterion
{
public:
	template<typename Scalar, typename Size>
	class Scope
	{
	public:
		Scope(const Size size, const Scalar& value = Scalar(1));
		
		Scope(const Scalar* values, const Size size);
		
		virtual ~Scope();
		
		Scope(const Scope&) = delete;
	
		Scope(Scope&&) = delete;
	
		Scope& operator=(const Scope&) = delete;
	
		Scope& operator=(Scope&&) = delete;
		
		[[nodiscard]] 
		inline constexpr static Size getActiveSize() { LNOT_ASSERT(s_activeSize == 0, "No active weights — wrap your call in a WeightedL2NormConvergenceCriterion::Scope"); return s_activeSize; }
		
		[[nodiscard]] 
		inline constexpr static Scalar* getActiveWeightsPtr() { LNOT_ASSERT(s_activeWeights.get(), "No active weights — wrap your call in a WeightedL2NormConvergenceCriterion::Scope"); return s_activeWeights.get(); }
		
		[[nodiscard]] 
		inline constexpr static const Scalar* getActiveConstWeightsPtr() { LNOT_ASSERT(s_activeWeights.get(), "No active weights — wrap your call in a WeightedL2NormConvergenceCriterion::Scope"); return s_activeWeights.get(); }
				
		[[nodiscard]] 
		inline constexpr static std::span<Scalar> getActiveWeights() { LNOT_ASSERT(s_activeWeights.get(), "No active weights — wrap your call in a WeightedL2NormConvergenceCriterion::Scope"); return std::span(s_activeWeights.get(), s_activeSize); }
		
		[[nodiscard]] 
		inline constexpr static std::span<const Scalar> getActiveConstWeights() { LNOT_ASSERT(s_activeWeights.get(), "No active weights — wrap your call in a WeightedL2NormConvergenceCriterion::Scope"); return std::span(s_activeWeights.get(), s_activeSize); }
	private:
		std::unique_ptr<Scalar[]> m_previousWeights;
		Size                      m_previousSize;
		
		inline static thread_local std::unique_ptr<Scalar[]> s_activeWeights;
		inline static thread_local Size                      s_activeSize;
	};

	template<typename Scalar> 
	constexpr inline Scalar getAbsTol(const Scalar& tol) const { return tol*tol; }
	
	template<typename Scalar> 
	constexpr inline Scalar getRelTol(const Scalar& tol, const Scalar& residual) const { return tol*tol*residual; }
	
	template<typename Scalar, typename Size>
	constexpr inline Scalar getResidual(const Scalar* gradient, const Size size) const { return BasicLinalg::weightedSquaredNorm(gradient, Scope<Scalar, Size>::getActiveConstWeightsPtr(), size); }
};

} // namespace LNOT

#include <LNOT/ConvergenceCriterions/WeightedL2NormConvergenceCriterion_impl.hpp>

#endif // LNOT_L2_NORM_CONVERGENCE_CRITERION_HPP
