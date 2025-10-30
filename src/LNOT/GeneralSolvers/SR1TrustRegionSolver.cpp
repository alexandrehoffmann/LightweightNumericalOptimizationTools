#include <LNOT/GeneralSolvers/SR1TrustRegionSolver.hpp>

namespace LNOT
{
	
template class SR1TrustRegionSolver< TruncatedConjugateGradient<float> >;
template class SR1TrustRegionSolver< LanczosTRSSolver<float> >;

template class SR1TrustRegionSolver< TruncatedConjugateGradient<double> >;
template class SR1TrustRegionSolver< LanczosTRSSolver<double> >;
	
} // namespace LNOT
