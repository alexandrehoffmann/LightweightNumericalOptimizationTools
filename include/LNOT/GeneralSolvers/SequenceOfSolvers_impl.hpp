#ifndef LNOT_SEQUENCE_OF_SOLVERS_IMPL_HPP
#define LNOT_SEQUENCE_OF_SOLVERS_IMPL_HPP

#include <LNOT/GeneralSolvers/SequenceOfSolvers.hpp>

#include <ranges>

namespace LNOT
{

template<CSolver... Solvers>
SequenceOfSolvers<Solvers...>::SequenceOfSolvers(const Size maxIt, const Scalar relTol, const Scalar absTol)
	: Base(maxIt, relTol, absTol)
{
	BIC::foreach(BIC::fixed<size_t, 0>, BIC::fixed<size_t, nSolvers>, [&maxIt, &relTol, &absTol, &solvers = m_solvers](const auto i)
	{
		std::get<i>(solvers).setMaxIt(maxIt / nSolvers);
		std::get<i>(solvers).setRelTol(relTol);
		std::get<i>(solvers).setAbsTol(absTol);
	});
}

template<CSolver... Solvers>
SequenceOfSolvers<Solvers...>::SequenceOfSolvers(const std::array<Size, nSolvers>& maxIt, const std::array<Scalar, nSolvers>& relTol, const std::array<Scalar, nSolvers>& absTol)
	: Base(std::reduce(maxIt.begin(), maxIt.end()), std::ranges::min(relTol), std::ranges::min(absTol))
{
	BIC::foreach(BIC::fixed<size_t, 0>, BIC::fixed<size_t, nSolvers>, [&maxIt, &relTol, &absTol, &solvers = m_solvers](const auto i)
	{
		std::get<i>(solvers).setMaxIt(maxIt[i]);
		std::get<i>(solvers).setRelTol(relTol[i]);
		std::get<i>(solvers).setAbsTol(absTol[i]);
	});
}

template<CSolver... Solvers>
void SequenceOfSolvers<Solvers...>::setMaxIt(const std::array<Size, nSolvers>& maxIt)
{
	BIC::foreach(BIC::fixed<size_t, 0>, BIC::fixed<size_t, nSolvers>, [&maxIt, &solvers = m_solvers](const auto i)
	{
		std::get<i>(solvers).setMaxIt(maxIt[i]);
	});
}

template<CSolver... Solvers>
void SequenceOfSolvers<Solvers...>::setRelTol(const std::array<Scalar, nSolvers>& tol)
{
	BIC::foreach(BIC::fixed<size_t, 0>, BIC::fixed<size_t, nSolvers>, [&tol, &solvers = m_solvers](const auto i)
	{
		std::get<i>(solvers).setRelTol(tol[i]);
	});
}

template<CSolver... Solvers>
void SequenceOfSolvers<Solvers...>::setAbsTol(const std::array<Scalar, nSolvers>& tol)
{
	BIC::foreach(BIC::fixed<size_t, 0>, BIC::fixed<size_t, nSolvers>, [&tol, &solvers = m_solvers](const auto i)
	{
		std::get<i>(solvers).setAbsTol(tol[i]);
	});
}

template<CSolver... Solvers> template<size_t I, CFirstOrderOracle Oracle, typename ABool> 
void SequenceOfSolvers<Solvers...>::solveImplRec(BIC::Fixed<size_t, I> i, Oracle& oracle, const ABool solveInPlace, Scalar* x) requires(IsBool<ABool>::value)
{
	using IthSolverInfo = typename IthSolver<I>::Info;
	
	IthSolver<I>& solver = std::get<i>(m_solvers);
	
	if constexpr (i == 0) { m_innerIts.clear(); }
	
	solver.solve(oracle, solveInPlace, x);
	m_innerIts.push_back(solver.getIterations());
	
	if (i == 0 and m_out) { fmt::println(m_out, "#SequenceOfSolvers\n#Solver #Iterations f(x) residual relative_tol absolute_tol"); }
	if (m_out)            { fmt::println(m_out, "{} {} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", i+1, solver.getIterations(), solver.getValue(), solver.getError(), solver.getRelTol(), solver.getAbsTol()); std::fflush(m_out); }
	
	if (solver.getInfo() != IthSolverInfo::SUCCESS)
	{
		if constexpr (i+1 != nSolvers) { solveImplRec(BIC::next(i), oracle, BIC::fixed<bool, false>, x); }
		else
		{
			m_info     = Info::FAILURE;
			m_nIt      = i;
			m_fx       = solver.getValue(); 
			m_residual = solver.getError();
		}
	}
	else
	{
		m_info     = Info::SUCCESS;
		m_nIt      = i;
		m_fx       = solver.getValue(); 
		m_residual = solver.getError(); 
	}
}
	
} // namespace LNOT

#endif // LNOT_SEQUENCE_OF_SOLVERS_IMPL_HPP
