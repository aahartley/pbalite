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
