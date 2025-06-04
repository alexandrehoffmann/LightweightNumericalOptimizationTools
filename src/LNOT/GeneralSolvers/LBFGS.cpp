#include <LNOT/GeneralSolvers/LBFGS.hpp>

namespace LNOT
{

template class LBFGS< BisectionLineSearch<float> >;
template class LBFGS< BacktrackingLineSearch<float> >;

template class LBFGS< BisectionLineSearch<double> >;
template class LBFGS< BacktrackingLineSearch<double> >;

} // namespace LNOT
