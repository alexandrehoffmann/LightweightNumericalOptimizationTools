#include <LNOT/GeneralSolvers/NonLinearConjugateGradient.hpp>

namespace LNOT
{
	
//// for floats
template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::FLETCHER_REEVES>;
template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::POLAK_RIBIERE>;
template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::FLETCHER>;
template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::LIU_STOREY>;
template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::DAI_YUAN>;
template class NonLinearConjugateGradient<BisectionLineSearch<float>, NLCGUpdateStrategy::HAGER_ZHANG>;

template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::FLETCHER_REEVES>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::POLAK_RIBIERE>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::FLETCHER>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::LIU_STOREY>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::DAI_YUAN>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<float>, NLCGUpdateStrategy::HAGER_ZHANG>;
//// for doubles
template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::FLETCHER>;
template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::LIU_STOREY>;
template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::DAI_YUAN>;
template class NonLinearConjugateGradient<BisectionLineSearch<double>, NLCGUpdateStrategy::HAGER_ZHANG>;

template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::HESTENES_STIEFEL>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::FLETCHER_REEVES>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::POLAK_RIBIERE>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::FLETCHER>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::LIU_STOREY>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::DAI_YUAN>;
template class NonLinearConjugateGradient<BacktrackingLineSearch<double>, NLCGUpdateStrategy::HAGER_ZHANG>;
	
} // namespace LNOT
