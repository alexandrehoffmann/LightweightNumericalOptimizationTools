#include <LNOT/GeneralSolvers/LBFGS.hpp>

namespace LNOT
{

template class LBFGS< BisectionLineSearch<float> >;
template class LBFGS< BacktrackingLineSearch<float> >;

template class LBFGS< BisectionLineSearch<double> >;
template class LBFGS< BacktrackingLineSearch<double> >;

template class LBFGS< BisectionLineSearch<long double> >;
template class LBFGS< BacktrackingLineSearch<long double> >;

} // namespace LNOT
