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


// class WCSPHSolver : public GISolverBase
// {
//   public:
//     WCSPHSolver(SPHState& pq, Force& f, double vclamp, double aclamp);
//     ~WCSPHSolver(){}
    
//     void init();
//     void solve(const double dt);

//     void advance_velocity();
//     void advance_position();
//     void correct_density_error();
//     void correct_divergence_error();
//     void get_timestep();

//   private:
//     SPHState PQ;
//     Force force;
//     float velocity_clamp;
//     float acceleration_clamp;
//     float dt;


// };

// GISolver CreateWCSPHSolver( SPHState& pq, Force& f, float vel_clamp, float accel_clamp );


class WCSPHSolverWithCollisions : public GISolverBase
{
  public:
    WCSPHSolverWithCollisions(SPHState& pq, Force& f, double vclamp, double aclamp, ElasticCollisionHandler& coll, GISolver& sol);
    ~WCSPHSolverWithCollisions(){}
    
    void init();
    void solve(const double dt);

    void get_timestep();

    const float get_velocity_clamp() const { return velocity_clamp; }
    void set_velocity_clamp(const float& v ) { velocity_clamp = v; }

    const float get_acceleration_clamp() const { return acceleration_clamp; }
    void set_acceleration_clamp(const float& v ) { acceleration_clamp = v; }


  private:
    SPHState PQ;
    Force force;
    ElasticCollisionHandler& CS;
    GISolver solver;
    float velocity_clamp;
    float acceleration_clamp;
    float user_dt;
    float dt;


};

GISolver CreateWCSPHSolver( SPHState& pq, Force& f, float vel_clamp, float accel_clamp, ElasticCollisionHandler& cs, GISolver& sol );



}

#endif