#include <LNOT/GeneralSolvers/BFGS.hpp>

namespace LNOT
{

template class BFGS< BisectionLineSearch<float> >;
template class BFGS< BacktrackingLineSearch<float> >;

template class BFGS< BisectionLineSearch<double> >;
template class BFGS< BacktrackingLineSearch<double> >;

} // namespace LNOT
