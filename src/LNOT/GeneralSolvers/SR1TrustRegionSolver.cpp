#include <LNOT/GeneralSolvers/SR1TrustRegionSolver.hpp>

namespace LNOT
{
	
template class SR1TrustRegionSolver< TruncatedConjugateGradient<float> >;
template class SR1TrustRegionSolver< LanczosTRSSolver<float> >;

template class SR1TrustRegionSolver< TruncatedConjugateGradient<double> >;
template class SR1TrustRegionSolver< LanczosTRSSolver<double> >;

template class SR1TrustRegionSolver< TruncatedConjugateGradient<long double> >;
template class SR1TrustRegionSolver< LanczosTRSSolver<long double> >;
	
} // namespace LNOT
