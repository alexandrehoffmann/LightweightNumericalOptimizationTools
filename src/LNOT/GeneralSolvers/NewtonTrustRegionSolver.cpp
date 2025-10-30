#include <LNOT/GeneralSolvers/NewtonTrustRegionSolver.hpp>

namespace LNOT
{
	
template class NewtonTrustRegionSolver< TruncatedConjugateGradient<float> >;
template class NewtonTrustRegionSolver< LanczosTRSSolver<float> >;

template class NewtonTrustRegionSolver< TruncatedConjugateGradient<double> >;
template class NewtonTrustRegionSolver< LanczosTRSSolver<double> >;
	
} // namespace LNOT
