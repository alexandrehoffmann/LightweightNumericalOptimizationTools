#include <LNOT/TRSSolvers/LanczosFullOrthTRSSolver.hpp>

namespace LNOT
{
	
template class LanczosFullOrthTRSSolver<float, 1>;
template class LanczosFullOrthTRSSolver<double, 1>;
	
template class LanczosFullOrthTRSSolver<float, 2>;
template class LanczosFullOrthTRSSolver<double, 2>;

} // namespace LNOT
