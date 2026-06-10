#include <LNOT/GeneralSolvers/GradientDescent.hpp>

namespace LNOT
{

template class GradientDescent< BisectionLineSearch<float> >;
template class GradientDescent< BacktrackingLineSearch<float> >;
template class GradientDescent< NoLineSearch<float> >;

template class GradientDescent< BisectionLineSearch<double> >;
template class GradientDescent< BacktrackingLineSearch<double> >;
template class GradientDescent< NoLineSearch<double> >;

template class GradientDescent< BisectionLineSearch<long double> >;
template class GradientDescent< BacktrackingLineSearch<long double> >;
template class GradientDescent< NoLineSearch<long double> >;

} // namespace LNOT
