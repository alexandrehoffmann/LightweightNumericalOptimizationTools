#ifndef LNOT_NLCG_UPDATE_STRATEGY_HPP
#define LNOT_NLCG_UPDATE_STRATEGY_HPP

namespace LNOT
{
// from doi.org/10.1051/ro/2022213 
enum class NLCGUpdateStrategy {	HESTENES_STIEFEL, FLETCHER_REEVES, POLAK_RIBIERE, FLETCHER, LIU_STOREY, DAI_YUAN, HAGER_ZHANG }; 

template<NLCGUpdateStrategy UpdateStrategy>
using NLCGUpdateStrategy_constant = std::integral_constant<NLCGUpdateStrategy, UpdateStrategy>;

using HestenesStiefelStrategy = NLCGUpdateStrategy_constant<NLCGUpdateStrategy::HESTENES_STIEFEL>;
using FletcherReevesStrategy  = NLCGUpdateStrategy_constant<NLCGUpdateStrategy::FLETCHER_REEVES>;
using PolakRibiereStrategy    = NLCGUpdateStrategy_constant<NLCGUpdateStrategy::POLAK_RIBIERE>;
using FletcherStrategy        = NLCGUpdateStrategy_constant<NLCGUpdateStrategy::FLETCHER>;
using LiuStoreyStrategy       = NLCGUpdateStrategy_constant<NLCGUpdateStrategy::LIU_STOREY>;
using DaiYuanStrategy         = NLCGUpdateStrategy_constant<NLCGUpdateStrategy::DAI_YUAN>;
using HagerZhangStrategy      = NLCGUpdateStrategy_constant<NLCGUpdateStrategy::HAGER_ZHANG>;
	
} // namespace LNOT

#endif // LNOT_NLCG_UPDATE_STRATEGY_HPP
