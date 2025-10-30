#include <LNOT/GeneralSolvers/LSR1TrustRegionSolver.hpp>

namespace LNOT
{

template class LSR1TrustRegionSolver< TruncatedConjugateGradient<float> >;
template class LSR1TrustRegionSolver< LanczosTRSSolver<float> >;

template class LSR1TrustRegionSolver< TruncatedConjugateGradient<double> >;
template class LSR1TrustRegionSolver< LanczosTRSSolver<double> >;

} // namespace LNOT
