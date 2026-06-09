#ifndef LNOT_CCONVERGENCE_CRITERION_FOR_HPP
#define LNOT_CCONVERGENCE_CRITERION_FOR_HPP

#include <concepts>

namespace LNOT
{

/**
 * @concept CConvergenceCriterionFor
 * @brief Defines the requirements for a convergence policy.
 *
 * A convergence policy specifies how convergence is measured by an
 * optimization algorithm. It defines:
 *
 * - how a residual is computed from a gradient vector,
 * - how an absolute convergence threshold is obtained,
 * - how a relative convergence threshold is obtained.
 *
 * The optimization algorithm uses these quantities to determine whether
 * an iterative process has converged.
 *
 * A model of this concept must be default constructible and provide
 * the following operations:
 *
 * @code
 * policy.getAbsTol(tol);
 * policy.getRelTol(tol, initialResidual);
 * policy.getResidual(gradient, size);
 * @endcode
 *
 * where:
 *
 * - `gradient` points to a gradient vector of length `size`,
 * - `initialResidual` is the reference residual, typically computed
 *   at the beginning of the optimization process,
 * - `tol` is a user-supplied tolerance value.
 *
 * The return values of `getAbsTol()` and `getRelTol()` define the
 * convergence thresholds, while `getResidual()` computes the current
 * residual according to the policy-specific metric.
 *
 * Typical implementations use the L1, L2, or L∞ norm of the gradient,
 * although custom residual metrics may also be provided.
 *
 * @tparam Policy Convergence policy type.
 * @tparam Scalar Scalar type used for residuals and tolerances.
 * @tparam Size Size type used.
 */
template<typename Policy, typename Scalar>
concept CConvergenceCriterionFor = std::default_initializable<Policy> and requires(const Policy policy, const Scalar* gradient, const size_t size, const Scalar initialResidual, const Scalar tol)
{
	{ policy.getAbsTol(tol) } -> std::convertible_to<Scalar>;
	{ policy.getRelTol(tol, initialResidual) } -> std::convertible_to<Scalar>;
	{ policy.getResidual(gradient, size) } -> std::convertible_to<Scalar>;
};

} // namespace LNOT

#endif // LNOT_CCONVERGENCE_CRITERION_FOR_HPP
