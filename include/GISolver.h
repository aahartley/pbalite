#ifndef __PBA_GISOLVER_H__
#define __PBA_GISOLVER_H__

#include <memory>

namespace pba
{

class GISolverBase
{
  public:
    GISolverBase(){}
    virtual void init() = 0;
    virtual void solve(const double dt) = 0;
    virtual ~GISolverBase(){};
};
typedef std::shared_ptr<pba::GISolverBase> GISolver;
// A == P', B == V'

//semi implicit
class BackwardEulerSolver : public GISolverBase
{
  public:
    BackwardEulerSolver(GISolver& A, GISolver& B) :
      a(A),
      b(B)
      {}
    ~BackwardEulerSolver(){}
    void init(){a->init(); b->init();}
    void solve(const double dt)
    {
        b->solve(dt);
        a->solve(dt);
    }
  private:
    GISolver a;
    GISolver b;
};

//semi-implicit
class ForwardEulerSolver : public GISolverBase
{
  public:
    ForwardEulerSolver(GISolver& A, GISolver& B) :
      a(A),
      b(B)
    {}
    ~ForwardEulerSolver(){}
    void init(){a->init(); b->init();}
    void solve(const double dt)
    {
        a->solve(dt);
        b->solve(dt);
    }
  private:
    GISolver a;
    GISolver b;
};

//vel,pos,vel
class LeapFrogSolver : public GISolverBase
{
  public:
    LeapFrogSolver(GISolver& A, GISolver& B) :
      a(A),
      b(B)
    {}
    ~LeapFrogSolver(){}
    void init(){a->init(); b->init();}
    void solve(const double dt)
    {
        const double dtd2 = 0.5*dt;
        a->solve(dtd2);
        b->solve(dt);
        a->solve(dtd2);
    }
  private:
    GISolver a;
    GISolver b;
};

GISolver CreateLeapFrogSolver(GISolver& A, GISolver& B);
GISolver CreateForwardEulerSolver(GISolver& A, GISolver& B);
GISolver CreateBackwardEulerSolver(GISolver& A, GISolver& B);


}


#endif