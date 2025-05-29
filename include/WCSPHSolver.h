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
    WCSPHSolver(SPHState& pq, Force& f, GISolver& sol);
    ~WCSPHSolver(){}
    
    void init();
    void solve(const double dt);

    void get_timestep();

  private:
    SPHState PQ;
    Force force;
    GISolver solver;
    float user_dt;
    float dt;


};

GISolver CreateWCSPHSolver( SPHState& pq, Force& f, GISolver& sol );


class WCSPHSolverWithCollisions : public GISolverBase
{
  public:
    WCSPHSolverWithCollisions(SPHState& pq, Force& f, ElasticCollisionHandler& coll, GISolver& sol);
    ~WCSPHSolverWithCollisions(){}
    
    void init();
    void solve(const double dt);

    void get_timestep();




  private:
    SPHState PQ;
    Force force;
    ElasticCollisionHandler& CS;
    GISolver solver;
    float user_dt;
    float dt;


};

GISolver CreateWCSPHSolver( SPHState& pq, Force& f, ElasticCollisionHandler& cs, GISolver& sol );



}

#endif