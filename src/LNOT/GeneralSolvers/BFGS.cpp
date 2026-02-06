#include <LNOT/GeneralSolvers/BFGS.hpp>

namespace LNOT
{

template class BFGS< BisectionLineSearch<float> >;
template class BFGS< BacktrackingLineSearch<float> >;
template class BFGS< NoLineSearch<float> >;

template class BFGS< BisectionLineSearch<double> >;
template class BFGS< BacktrackingLineSearch<double> >;
template class BFGS< NoLineSearch<double> >;

template class BFGS< BisectionLineSearch<long double> >;
template class BFGS< BacktrackingLineSearch<long double> >;
template class BFGS< NoLineSearch<long double> >;

} // namespace LNOT
