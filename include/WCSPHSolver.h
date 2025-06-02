#ifndef ____PBA_WCSPHSOLVER_H____
#define ____PBA_WCSPHSOLVER_H____

#include "SPHState.h"
#include "Force.h"
#include "GISolver.h"
#include "CollisionHandler.h"
#include "ExplicitDynamics.h"
#include <chrono>


namespace pba
{


class WCSPHSolver : public GISolverBase
{
  public:
    WCSPHSolver(SPHStateData& pq, ForceBase<SPHStateData>& f, GISolverPtr solver)
      : pq_(pq),
        force_(f),
        solver_(std::move(solver)),
        user_dt_(0.001),
        dt_(0.001) {}

    ~WCSPHSolver(){}
    
    void Init();
    void Solve(double dt);

    void get_timestep();

  private:
    SPHStateData& pq_;
    ForceBase<SPHStateData>& force_;
    GISolverPtr solver_;
    float user_dt_;
    float dt_;



};

GISolverPtr CreateWCSPHSolver(SPHStateData& pq, ForceBase<SPHStateData>& f, GISolverPtr solver);





}

#endif