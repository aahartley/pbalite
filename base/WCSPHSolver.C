//-------------------------------------------------------
//
//  WCSPHSolver.C
//
//  Solvers for DFSPH Dynamics
//
//  DFSPH functions are based off the 2015 paper: https://animation.rwth-aachen.de/media/papers/2015-SCA-DFSPH.pdf
//  2017 paper, where the equations are written differently (still equivalent): https://animation.rwth-aachen.de/media/papers/2017-TVCG-ViscousDFSPH.pdf
//
//  Copyright (c) 2024 Jerry Tessendorf
//
//
//--------------------------------------------------------

#include "WCSPHSolver.h"
#include "ForceLibrary.h"

namespace pba
{

//sim loop does this already --- no need call
//requirements for sim to start
void WCSPHSolver::Init()
{
    pq_.Populate();
    pq_.ComputeDensity();
}


void WCSPHSolver::Solve(const double userdt)
{
    user_dt_ = userdt;

    //Occupancy Grid
    pq_.Populate(); 

    pq_.ComputeDensity();

    solver_->Solve(dt_);
    //CFL Condition
    get_timestep();

    //auto start2 = std::chrono::high_resolution_clock::now();

   //  auto end2 = std::chrono::high_resolution_clock::now();
   //  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
   //  std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;

    // DynamicalState ss =  std::dynamic_pointer_cast<DynamicalStateData, SPHStateData>(pq_);
    // CS.handle_collisions( dt, ss );
    // std::cout << "collissions good\n";
}

//CFL condition
void WCSPHSolver::get_timestep()
{
    if(pq_.get_use_user_dt())
    {
      dt_ = user_dt_;
      std::cout << "Dt: " << dt_ << '\n';
      return;
    }
    double pdiameter = pq_.get_radius() / 2.0; //particle radius is kernel_radius /4
    float max_vel = pq_.MaxVelocity();
    float lambda = 0.4f;
    double max_dt;
    if(max_vel != 0)
       max_dt = lambda * (pdiameter/max_vel);
    else max_dt = user_dt_;   

    dt_ = max_dt;
    //if(dt > user_dt_ || dt==0) dt = user_dt_; //user_dt_ dictactes max dt from cfl, change?
    //return (dt < 0.0001) ? 0.0001 : dt;

    //hard coding to clamp between 0.005 and 0.0001, can remove
    if(dt_ > 0.005) dt_ = 0.005;
    if(dt_ < 0.0001) dt_ = 0.0001;

    std::cout << "Dt: " << dt_ << '\n';

}

GISolverPtr CreateWCSPHSolver(SPHStateData& pq, ForceBase<SPHStateData>& f, GISolverPtr solver)
{
  return std::make_unique<WCSPHSolver>(pq, f, std::move(solver));
}

}
