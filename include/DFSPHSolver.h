//*******************************************************************
//
//   DFSPHSolver.h
//
//   DFSPH Solver with CFL condition
//   Two Pressure solvers
//
//
//
//*******************************************************************

//  DFSPH functions are based off the 2015 paper: https://animation.rwth-aachen.de/media/papers/2015-SCA-DFSPH.pdf
//  2017 paper, where the equations are written differently (still equivalent): https://animation.rwth-aachen.de/media/papers/2017-TVCG-ViscousDFSPH.pdf

#ifndef ____PBA_DFSPHSOLVER_H____
#define ____PBA_DFSPHSOLVER_H____

#include "SPHState.h"
#include "Force.h"
#include "Viscosity.h"
#include "GISolver.h"
#include "CollisionHandler.h"


namespace pba
{
/*!
  DFSPH solver 
  Two pressure solves: constant density, divergence free
 */  
class DFSPHSolver : public GISolverBase
{
  public:
    DFSPHSolver(
      SPHStateData& pq, ForceBase<SPHStateData>& f, ElasticCollisionHandler* coll = nullptr, 
      float aclamp = 1.0e9, float vclamp = 1.0e9)
      : pq_(pq),
        force_(f),
        coll_handler_(coll),
        use_coll_(coll_handler_ != nullptr),
        acceleration_clamp_(aclamp),
        velocity_clamp_(vclamp),
        user_dt_(0.001),
        dt_(0.001) {}

    ~DFSPHSolver(){}
    
    void Init();
    void Solve(double dt);
    //! Integrate velocities
    void advance_velocity();
    //! Integrate positions
    void advance_position();
    //! Enforce constant density: PPE 1 (Alg. 3): ρ∗i - ρ0 = 0
    void correct_density_error();
    //! Jacobi like iterativ solve until end of iterations or correct density error (ki)
    void density_solve_iteration(float& avg_density_error);
    //! Enforce divergence free velocity field: Dρi/Dt = 0
    void correct_divergence_error();
    //! Jacobi like iterativ solve until end of iterations or correct density error (kiv)
    void divergence_solve_iteration(float& avg_density_error);
    //! Compute pressure accelerations using ki/v
    void compute_pressure_acc(size_t p, const std::string& type);
    //! Compute error in density due to pressure
    float compute_error_force(size_t p, const std::string& type);
    //! CFL condition or user input
    void get_timestep();
    //! Temporary collision detection test
    void fakecs();

    float get_velocity_clamp() const { return velocity_clamp_; }
    void set_velocity_clamp(float v) { velocity_clamp_ = v; }

    float get_acceleration_clamp() const { return acceleration_clamp_; }
    void set_acceleration_clamp(float a) { acceleration_clamp_ = a; }

  private:
    SPHStateData& pq_;
    ForceBase<SPHStateData>& force_;
    ElasticCollisionHandler* coll_handler_;
    bool use_coll_;
    float acceleration_clamp_;
    float velocity_clamp_;
    float user_dt_;
    float dt_;
};

inline GISolverPtr CreateDFSPHSolver(
  SPHStateData& pq, ForceBase<SPHStateData>& f, ElasticCollisionHandler* coll = nullptr,
  float accel_clamp = 1.0e9, float vel_clamp = 1.0e9)
{
  return std::make_unique<DFSPHSolver>(pq, f, coll, accel_clamp, vel_clamp);
}


}// end of pba namespace

#endif // ____PBA_DFSPHSOLVER_H____