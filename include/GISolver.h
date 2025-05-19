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
        b->solve(dtd2);
        a->solve(dt);
        b->solve(dtd2);
    }
  private:
    GISolver a;
    GISolver b;
};
//! Implements multiple substeps of specified solver
class GISolverSubstep : public GISolverBase
{
  public:
    GISolverSubstep( GISolver& s, int nbsteps ) : 
        _solver (s), _steps (nbsteps)
	{}

    ~GISolverSubstep(){}

    void init(){ _solver->init(); }

    void solve( const double dt )
    {
        const double dta = dt/_steps;
	for( int i=0;i<_steps;i++){ _solver->solve(dta); }
    }

  private:

    GISolver _solver;
    double _steps;
};



//! Fourth order accurate
class GISolverFourthOrder : public GISolverBase
{
  public:
    GISolverFourthOrder( GISolver& s ) : 
        _solver (s) 
	{
	    _a = 1.0/( 2.0 - std::pow(2.0, 1.0/3.0) );
            _b = 1.0 - 2.0*_a;
	}

    ~GISolverFourthOrder(){}

    void init(){ _solver->init(); }

    void solve( const double dt )
    {
        const double dta = _a * dt;
	const double dtb = _b * dt;
	_solver->solve(dta);
	_solver->solve(dtb);
	_solver->solve(dta);
    }


  private:

    GISolver _solver;
    double _a,_b;
};


//! Sixth order accurate
class GISolverSixthOrder : public GISolverBase
{
  public:
    GISolverSixthOrder( GISolver& s ) : 
        _solver (s) 
	{
	    _a = 1.0/( 4.0 - std::pow(4.0, 1.0/3.0) );
            _b = 1.0 - 4.0*_a;
	}

    ~GISolverSixthOrder(){}

    void init(){ _solver->init(); }

    void solve( const double dt )
    {
        const double dta = _a * dt;
	const double dtb = _b * dt;
	_solver->solve(dta);
	_solver->solve(dta);
	_solver->solve(dtb);
	_solver->solve(dta);
	_solver->solve(dta);
    }


  private:

    GISolver _solver;
    double _a, _b;
};

GISolver CreateLeapFrogSolver(GISolver& A, GISolver& B);
GISolver CreateForwardEulerSolver(GISolver& A, GISolver& B);
GISolver CreateBackwardEulerSolver(GISolver& A, GISolver& B);
GISolver CreateGISolverSubstep( GISolver& s, int nbsteps );
GISolver CreateGISolverFourthOrder( GISolver& s );
GISolver CreateGISolverSixthOrder( GISolver& s );

}


#endif