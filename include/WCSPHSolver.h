//*******************************************************************
//
//   WCSPHSolver.h
//
//  Weakly Incompressible SPH solver
//
//
//
//*******************************************************************


#ifndef ____PBA_WCSPHSOLVER_H____
#define ____PBA_WCSPHSOLVER_H____

#include "SPHState.h"
#include "Force.h"
#include "GISolver.h"
#include "CollisionHandler.h"
#include "ExplicitDynamics.h"


namespace pba
{

/*!
  Weakly Inncompressible SPH solver
 */
class WCSPHSolver : public GISolverBase
{
  public:
    WCSPHSolver(SPHStateData& pq, GISolverPtr solver)
      : pq_(pq),
        solver_(std::move(solver)),
        user_dt_(0.001),
        dt_(0.001) {}

    ~WCSPHSolver(){}
    
    void Init();
    void Solve(double dt);

    void get_timestep();

  private:
    SPHStateData& pq_;
    GISolverPtr solver_;
    float user_dt_;
    float dt_;



};

GISolverPtr CreateWCSPHSolver(SPHStateData& pq, GISolverPtr solver);


} //end of pba namespace

#endif //____PBA_WCSPHSOLVER_H____