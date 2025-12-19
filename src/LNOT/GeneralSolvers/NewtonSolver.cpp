#include <LNOT/GeneralSolvers/NewtonSolver.hpp>

namespace LNOT
{
	
// explicit instanciation for float
template class NewtonSolver<ConjugateGradient<float>,    BisectionLineSearch<float> >;
template class NewtonSolver<LanczosSolver<float>,        BisectionLineSearch<float> >;

template class NewtonSolver<ConjugateGradient<float>,    BacktrackingLineSearch<float> >;
template class NewtonSolver<LanczosSolver<float>,        BacktrackingLineSearch<float> >;
// explicit instanciation for double
template class NewtonSolver<ConjugateGradient<double>,    BisectionLineSearch<double> >;
template class NewtonSolver<LanczosSolver<double>,        BisectionLineSearch<double> >;

template class NewtonSolver<ConjugateGradient<double>,    BacktrackingLineSearch<double> >;
template class NewtonSolver<LanczosSolver<double>,        BacktrackingLineSearch<double> >;
// explicit instanciation for long double
template class NewtonSolver<ConjugateGradient<long double>,    BisectionLineSearch<long double> >;
template class NewtonSolver<LanczosSolver<long double>,        BisectionLineSearch<long double> >;

template class NewtonSolver<ConjugateGradient<long double>,    BacktrackingLineSearch<long double> >;
template class NewtonSolver<LanczosSolver<long double>,        BacktrackingLineSearch<long double> >;

} // namespace LNOT
