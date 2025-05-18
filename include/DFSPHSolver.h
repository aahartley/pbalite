//  DFSPH functions are based off the 2015 paper: https://animation.rwth-aachen.de/media/papers/2015-SCA-DFSPH.pdf
//  2017 paper, where the equations are written differently (still equivalent): https://animation.rwth-aachen.de/media/papers/2017-TVCG-ViscousDFSPH.pdf


#ifndef ____PBA_DFSPHSOLVER_H____
#define ____PBA_DFSPHSOLVER_H____

#include "SPHState.h"
#include "Force.h"
#include "Viscosity.h"
#include "GISolver.h"
#include "CollisionHandler.h"
#include <chrono>


namespace pba
{


// class DFSPHSolver : public GISolverBase
// {
//   public:
//     DFSPHSolver(SPHState& pq, Force& f, double vclamp, double aclamp, DFSPHForce& sphf);
//     ~DFSPHSolver(){}
    
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
//     DFSPHForce sphforce;
//     float velocity_clamp;
//     float acceleration_clamp;
//     float dt;


// };

// GISolver CreateDFSPHSolver( SPHState& pq, Force& f, float vel_clamp, float accel_clamp, DFSPHForce& sphforce );


class DFSPHSolverWithCollisions : public GISolverBase
{
  public:
    DFSPHSolverWithCollisions(SPHState& pq, Force& f, double vclamp, double aclamp, ElasticCollisionHandler& coll);
    ~DFSPHSolverWithCollisions(){}
    
    void init();
    void solve(const double dt);

    void advance_velocity();
    void advance_position();
    void correct_density_error();
    void density_solve_iteration(float& avg_density_error);
    void correct_divergence_error();
    void divergence_solve_iteration(float& avg_density_error);
    void compute_pressure_acc(size_t p, const std::string& type);
    float compute_error_force(size_t p, const std::string& type);
    void get_timestep();
    void fakecs();

    const float get_velocity_clamp() const { return velocity_clamp; }
    void set_velocity_clamp(const float& v ) { velocity_clamp = v; }

    const float get_acceleration_clamp() const { return acceleration_clamp; }
    void set_acceleration_clamp(const float& v ) { acceleration_clamp = v; }


  private:
    SPHState PQ;
    Force force;
    ElasticCollisionHandler& CS;
    float velocity_clamp;
    float acceleration_clamp;
    float user_dt;
    float dt;


};

GISolver CreateDFSPHSolver( SPHState& pq, Force& f, float vel_clamp, float accel_clamp, ElasticCollisionHandler& cs );



}

#endif