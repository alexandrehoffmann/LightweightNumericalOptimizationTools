#include <LNOT/GeneralSolvers/NewtonSolver.hpp>

namespace LNOT
{
	
// explicit instanciation for float
template class NewtonSolver<ConjugateGradient<float>,    BisectionLineSearch<float> >;
template class NewtonSolver<LanczosSolver<float>,        BisectionLineSearch<float> >;
template class NewtonSolver<CoupledLanczosSolver<float>, BisectionLineSearch<float> >;

template class NewtonSolver<ConjugateGradient<float>,    BacktrackingLineSearch<float> >;
template class NewtonSolver<LanczosSolver<float>,        BacktrackingLineSearch<float> >;
template class NewtonSolver<CoupledLanczosSolver<float>, BacktrackingLineSearch<float> >;
// explicit instanciation for double
template class NewtonSolver<ConjugateGradient<double>,    BisectionLineSearch<double> >;
template class NewtonSolver<LanczosSolver<double>,        BisectionLineSearch<double> >;
template class NewtonSolver<CoupledLanczosSolver<double>, BisectionLineSearch<double> >;

template class NewtonSolver<ConjugateGradient<double>,    BacktrackingLineSearch<double> >;
template class NewtonSolver<LanczosSolver<double>,        BacktrackingLineSearch<double> >;
template class NewtonSolver<CoupledLanczosSolver<double>, BacktrackingLineSearch<double> >;

} // namespace LNOT
