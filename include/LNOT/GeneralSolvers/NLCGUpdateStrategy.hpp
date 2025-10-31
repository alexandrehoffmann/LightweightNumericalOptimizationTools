#ifndef LNOT_NLCG_UPDATE_STRATEGY_HPP
#define LNOT_NLCG_UPDATE_STRATEGY_HPP

namespace LNOT
{
	
/**
 * @brief Enumeration indicating the update strategy of the Non Linear Conjugate Gradient (NLCG)
 * 
 * The NLCG computes it's next step as follows : \f$d_{k+1} = -\nabla f(x_{k+1}) + \beta_k d_k\f$
 * Multiple update parameters \f$\beta_k\f$ have been proposed in the literature. See Stanimirović et al. for more details
 * 
 * Stanimirović, P. S., Ivanov, B., Ma, H., & Mosić, D. (2020). A survey of gradient methods for solving nonlinear optimization. Electronic Research Archive, 28(4), 1573–1624. doi:10.3934/era.2020115
 */
enum class NLCGUpdateStrategy 
{	
	HESTENES_STIEFEL, ///<  @brief \f$\beta_k = \frac{g_{k+1}^T y_k}{d_k^T y_k}\f$
	FLETCHER_REEVES,  ///<  @brief \f$\beta_k = \frac{\|g_{k+1}\|^2}{\|g_{k}\|^2}\f$
	POLAK_RIBIERE,    ///<  @brief \f$\beta_k = \frac{g_{k+1}^T y_k}{\|g_{k}\|^2}\f$
	FLETCHER,         ///<  @brief \f$\beta_k = \frac{\|g_{k+1}\|^2}{-d_k^Tg_k}\f$
	LIU_STOREY,       ///<  @brief \f$\beta_k = \frac{g_{k+1}^Ty_k}{-d_k^Tg_k}\f$
	DAI_YUAN,         ///<  @brief \f$\beta_k = \frac{\|g_{k+1}\|^2}{d_k^T y_k}\f$
	HAGER_ZHANG       ///<  @brief \f$\beta_k = \left(y_k - 2d_k\frac{\|y_k\|^2}{d_k^Ty_k}\right)^T\frac{g_{k+1}}{d_k^T y_k}\f$
}; 
	
} // namespace LNOT

#endif // LNOT_NLCG_UPDATE_STRATEGY_HPP
