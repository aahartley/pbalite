//*******************************************************************
//
//   GISolver.h
//
//  Bounding Volume Hierarchy for continous collison detection
//
//
//
//*******************************************************************

#ifndef __PBA_GISOLVER_H__
#define __PBA_GISOLVER_H__

#include <memory>
#include <cmath>

namespace pba
{

class GISolverBase
{
  public:
    GISolverBase(){}
    virtual void Init() = 0;
    virtual void Solve(double dt) = 0;
    virtual ~GISolverBase(){};
};
using GISolverPtr = std::unique_ptr<pba::GISolverBase>;

// a == P', b == V'

//semi implicit
class BackwardEulerSolver : public GISolverBase
{
  public:
    BackwardEulerSolver(GISolverBase* a, GISolverBase* b)
      : a_(a),
        b_(b) {}
    ~BackwardEulerSolver() {}
    void Init() { a_->Init(); b_->Init(); }
    void Solve(double dt)
    {
      if (dt >= 0)
      {
        b_->Solve(dt);
        a_->Solve(dt);
      }
      else
      {
        //not use negative future vel for postion update, update pos first for -dt
        a_->Solve(dt);
        b_->Solve(dt);
      }
    }
  private:
    GISolverBase* a_;
    GISolverBase* b_;
};

//semi-implicit
class ForwardEulerSolver : public GISolverBase
{
  public:
    ForwardEulerSolver(GISolverBase* a, GISolverBase* b) 
     : a_(a),
       b_(b) {}
    ~ForwardEulerSolver() {}
    void Init() { a_->Init(); b_->Init(); }
    void Solve(double dt)
    {
        a_->Solve(dt);
        b_->Solve(dt);
    }
  private:
    GISolverBase* a_;
    GISolverBase* b_;
};

//vel,pos,vel
class LeapFrogSolver : public GISolverBase
{
  public:
    LeapFrogSolver(GISolverBase* a, GISolverBase* b)
      : a_(a),
        b_(b) {}
    ~LeapFrogSolver() {}
    void Init(){ a_->Init(); b_->Init(); }
    void Solve(double dt)
    {
        const double dtd2 = 0.5*dt;
        b_->Solve(dtd2);
        a_->Solve(dt);
        b_->Solve(dtd2);
    }
  private:
    GISolverBase* a_;
    GISolverBase* b_;
};
//! Implements multiple substeps of specified solver
class GISolverSubstep : public GISolverBase
{
  public:
    GISolverSubstep(GISolverPtr s, int nbsteps) 
      : solver_(std::move(s)),
        steps_(nbsteps) {}

    ~GISolverSubstep() {}
    void Init(){ solver_->Init(); }
    void Solve(double dt)
    {
      const double dta = dt/steps_;
	    for (int i = 0; i <steps_; ++i)
      { 
        solver_->Solve(dta);
      }
    }
  private:
    GISolverPtr solver_;
    double steps_;
};


//! Fourth order accurate
class GISolverFourthOrder : public GISolverBase
{
  public:
    GISolverFourthOrder(GISolverPtr s)
      : solver_(std::move(s)) 
	  {
	    a_ = 1.0/( 2.0 - std::pow(2.0, 1.0/3.0) );
      b_ = 1.0 - 2.0*a_;
	  }
    ~GISolverFourthOrder() {}
    void Init(){ solver_->Init(); }
    void Solve(double dt)
    {
      const double dta = a_ * dt;
	    const double dtb = b_ * dt;
	    solver_->Solve(dta);
	    solver_->Solve(dtb);
	    solver_->Solve(dta);
    }

  private:
    GISolverPtr solver_;
    double a_, b_;
};


//! Sixth order accurate (goes back in time)
class GISolverSixthOrder : public GISolverBase
{
  public:
    GISolverSixthOrder(GISolverPtr s) 
      : solver_(std::move(s)) 
	  {
	    a_ = 1.0/( 4.0 - std::pow(4.0, 1.0/3.0) );
      b_ = 1.0 - 4.0*a_;
	  }
    ~GISolverSixthOrder(){}

    void Init(){ solver_->Init(); }
    void Solve(double dt)
    {
      const double dta = a_ * dt;
      const double dtb = b_ * dt;
      solver_->Solve(dta);
      solver_->Solve(dta);
      solver_->Solve(dtb);
      solver_->Solve(dta);
      solver_->Solve(dta);
    }
  private:
    GISolverPtr solver_;
    double a_, b_;
};

inline GISolverPtr CreateLeapFrogSolver(GISolverBase* a, GISolverBase* b)
{
  return std::make_unique<LeapFrogSolver>(a, b);
}
inline GISolverPtr CreateForwardEulerSolver(GISolverBase* a, GISolverBase* b)
{
  return std::make_unique<ForwardEulerSolver>(a, b);
}
inline GISolverPtr CreateBackwardEulerSolver(GISolverBase* a, GISolverBase* b)
{
  return std::make_unique<BackwardEulerSolver>(a, b);
}
inline GISolverPtr CreateGISolverSubstep(GISolverPtr s, int nb_steps)
{
  return std::make_unique<GISolverSubstep>(std::move(s), nb_steps);
}
inline GISolverPtr CreateGISolverFourthOrder(GISolverPtr s)
{
  return std::make_unique<GISolverFourthOrder>(std::move(s));
}
inline GISolverPtr CreateGISolverSixthOrder(GISolverPtr s)
{
  return std::make_unique<GISolverSixthOrder>(std::move(s));
}

}


#endif