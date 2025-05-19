#include "GISolver.h"

using namespace pba;

GISolver pba::CreateBackwardEulerSolver(GISolver& A, GISolver& B)
{
  return std::make_shared<BackwardEulerSolver>(A, B);
}
GISolver pba::CreateForwardEulerSolver(GISolver& A, GISolver& B)
{
  return std::make_shared<ForwardEulerSolver>(A, B);
}
GISolver pba::CreateLeapFrogSolver(GISolver& A, GISolver& B)
{
  return std::make_shared<LeapFrogSolver>(A, B);
}
GISolver pba::CreateGISolverSubstep(GISolver& s, int nbsteps )
{
  return std::make_shared<GISolverSubstep>(s, nbsteps);
}GISolver pba::CreateGISolverFourthOrder(GISolver& A)
{
  return std::make_shared<GISolverFourthOrder>(A);
}GISolver pba::CreateGISolverSixthOrder(GISolver& A)
{
  return std::make_shared<GISolverSixthOrder>(A);
}

