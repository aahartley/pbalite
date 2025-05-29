//-------------------------------------------------------
//
//  DFSPHSPHSolver.C
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

using namespace pba;
using namespace std;

//colls controlled by solver passed to it***, but to keep design pattern
WCSPHSolver::WCSPHSolver(SPHState& pq, Force& f, GISolver& sol) :
  PQ(pq),
  force(f),
  solver(sol),
  dt(0.001)
{  
}

//sim loop does this already --- no need call
//requirements for sim to start
void WCSPHSolver::init()
{
    PQ->populate();
    PQ->compute_density();
}


void WCSPHSolver::solve(const double userdt)
{
    user_dt = userdt;

    //Occupancy Grid
    PQ->populate(); 

    PQ->compute_density();

    solver->solve(dt);
    //CFL Condition
    get_timestep();

    //auto start2 = std::chrono::high_resolution_clock::now();

   //  auto end2 = std::chrono::high_resolution_clock::now();
   //  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
   //  std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;

    // DynamicalState ss =  std::dynamic_pointer_cast<DynamicalStateData, SPHStateData>(PQ);
    // CS.handle_collisions( dt, ss );
    // std::cout << "collissions good\n";
}

//CFL condition
void WCSPHSolver::get_timestep()
{
    if(PQ->get_useUserDT())
    {
      dt = user_dt;
      std::cout << "Dt: " << dt << '\n';
      return;
    }
    double pdiameter = PQ->get_radius() / 2.0; //particle radius is kernel_radius /4
    float max_vel = PQ->max_velocity();
    float lambda = 0.4f;
    double max_dt;
    if(max_vel != 0)
       max_dt = lambda * (pdiameter/max_vel);
    else max_dt = user_dt;   

    dt = max_dt;
    //if(dt > user_dt || dt==0) dt = user_dt; //user_dt dictactes max dt from cfl, change?
    //return (dt < 0.0001) ? 0.0001 : dt;

    //hard coding to clamp between 0.005 and 0.0001, can remove
    if(dt > 0.005) dt = 0.005;
    if(dt < 0.0001) dt = 0.0001;

    std::cout << "Dt: " << dt << '\n';

}

pba::GISolver pba::CreateWCSPHSolver( SPHState& pq, Force& f, GISolver& sol )
{
   return GISolver( new WCSPHSolver( pq, f, sol ) );
}



WCSPHSolverWithCollisions::WCSPHSolverWithCollisions(SPHState& pq, Force& f, ElasticCollisionHandler& coll, GISolver& sol) :
  PQ(pq),
  force(f),
  CS(coll),
  solver(sol),
  dt(0.001)
{  
}

//sim loop does this already --- no need call
//requirements for sim to start
void WCSPHSolverWithCollisions::init()
{
    PQ->populate();
    PQ->compute_density();
}


void WCSPHSolverWithCollisions::solve(const double userdt)
{
    user_dt = userdt;

    //Occupancy Grid
    PQ->populate(); 

    PQ->compute_density();

    solver->solve(dt);
    //CFL Condition
    get_timestep();

    //auto start2 = std::chrono::high_resolution_clock::now();

   //  auto end2 = std::chrono::high_resolution_clock::now();
   //  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
   //  std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;

    // DynamicalState ss =  std::dynamic_pointer_cast<DynamicalStateData, SPHStateData>(PQ);
    // CS.handle_collisions( dt, ss );
    // std::cout << "collissions good\n";
}

//CFL condition
void WCSPHSolverWithCollisions::get_timestep()
{
    if(PQ->get_useUserDT())
    {
      dt = user_dt;
      std::cout << "Dt: " << dt << '\n';
      return;
    }
    double pdiameter = PQ->get_radius() / 2.0; //particle radius is kernel_radius /4
    float max_vel = PQ->max_velocity();
    float lambda = 0.4f;
    double max_dt;
    if(max_vel != 0)
       max_dt = lambda * (pdiameter/max_vel);
    else max_dt = user_dt;   

    dt = max_dt;
    //if(dt > user_dt || dt==0) dt = user_dt; //user_dt dictactes max dt from cfl, change?
    //return (dt < 0.0001) ? 0.0001 : dt;

    //hard coding to clamp between 0.005 and 0.0001, can remove
    if(dt > 0.005) dt = 0.005;
    if(dt < 0.0001) dt = 0.0001;

    std::cout << "Dt: " << dt << '\n';

}

pba::GISolver pba::CreateWCSPHSolver( SPHState& pq, Force& f, ElasticCollisionHandler& coll, GISolver& sol )
{
   return GISolver( new WCSPHSolverWithCollisions( pq, f, coll, sol ) );
}
