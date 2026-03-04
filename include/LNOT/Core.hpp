#include <LNOT/FloatingPoint/NumTraits.hpp>

#include <LNOT/BasicLinalg/BasicLinalg.hpp>
#include <LNOT/BasicLinalg/SymmetricDenseMatrixOp.hpp>
#include <LNOT/BasicLinalg/DiagonalPreconditionerOp.hpp>

#include <LNOT/LinearSolvers/ConjugateGradient.hpp>
#include <LNOT/LinearSolvers/LanczosSolver.hpp>

#include <LNOT/TRSSolvers/LanczosTRSSolver.hpp>
#include <LNOT/TRSSolvers/TruncatedConjugateGradient.hpp>

#include <LNOT/LineSearches/BacktrackingLineSearch.hpp>
#include <LNOT/LineSearches/BisectionLineSearch.hpp>

#include <LNOT/GeneralSolvers/NewtonSolver.hpp>
#include <LNOT/GeneralSolvers/NewtonTrustRegionSolver.hpp>
#include <LNOT/GeneralSolvers/BFGS.hpp>
#include <LNOT/GeneralSolvers/LBFGS.hpp>
#include <LNOT/GeneralSolvers/NonLinearConjugateGradient.hpp>
#include <LNOT/GeneralSolvers/SR1TrustRegionSolver.hpp>
#include <LNOT/GeneralSolvers/LSR1TrustRegionSolver.hpp>
