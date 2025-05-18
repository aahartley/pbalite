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

// DFSPHSolver::DFSPHSolver(SPHState& pq, Force& f, double vclamp, double aclamp, DFSPHForce& sphf) :
//   PQ(pq),
//   force(f),
//   sphforce(sphf),
//   velocity_clamp(vclamp),
//   acceleration_clamp(aclamp),
//   dt(1.0/48.0)
// {}

// void DFSPHSolver::init() 
// {
//     PQ->populate();
//     PQ->compute_density();
//     PQ->compute_factor();
// }

// void DFSPHSolver::solve(const double deltaTime)
// {
//     force->compute(PQ, dt);
//     get_timestep();
//     advance_velocity();
//     correct_density_error();
//     advance_position();
//     PQ->populate();
//     PQ->compute_density();
//     PQ->compute_factor();
//     correct_divergence_error();

// }

// void DFSPHSolver::correct_density_error()
// {
// //    float densAvg = PQ->average_density();
// //    float densRest = 1000;
// //    int iter = 0;
// //    float threshold = 100.f;
// //    while((((densAvg - densRest) > threshold) || iter < 2)&& iter < 100)
// //    {
// //       //predict dens
// //       PQ->compute_predicted_density(dt);
// //       #pragma omp parallel for
// //       for(size_t p = 0; p < PQ->nb(); p++)
// //       {
// //          const Vector P = PQ->pos(p);
// //          const Vector V = PQ->vel(p);
// //          float density = PQ->get_float_attr("density", p);
// //          size_t pindex = PQ->index(P);
// //          size_t i,j,k;
// //          PQ->anti_index(pindex, i,j,k);
// //          float ki = ((PQ->get_float_attr("predicted_density", p) - densRest) / (dt * dt)) * PQ->get_float_attr("factor", p);

// //          Vector sum{};
// // //cell
// //          const std::vector<size_t>& cells = PQ->cell_contents(i,j,k); 
// //          for(size_t a=0;a<cells.size();a++)
// //          {
// //             size_t pid = cells[a]; 
// //             float kj = ((PQ->get_float_attr("predicted_density", pid) - densRest) / (dt * dt)) * PQ->get_float_attr("factor", pid);
// //             float density_j = PQ->get_float_attr("density", pid);
// //             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
// //          }
// // //neighbors
// //          std::vector<size_t> neighbors;
// //          PQ->neighbor_cells(i,j,k,neighbors);
// //          for(size_t a=0;a<neighbors.size();a++)
// //          {
// //             size_t pid = cells[a]; 
// //             float kj = ((PQ->get_float_attr("predicted_density", pid) - densRest) / (dt * dt)) * PQ->get_float_attr("factor", pid);
// //             float density_j = PQ->get_float_attr("density", pid);
// //             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
// //          }

// //          PQ->set_vel(p, V -  (dt * sum));

// //       }
// //       //what to do about new avg??
// //       densAvg = PQ->average_predicted_density();
// //       iter++;
      
// //    }
// }

// void DFSPHSolver::correct_divergence_error()
// {
// //    PQ->compute_density_derivative(); //use from prev?
// //    float densDAvg = PQ->average_density_derivative();
// //    int iter = 0;
// //    float threshold = 10.f;
// //    while((densDAvg > threshold || iter < 1) && iter < 100)   {
// //       //predict dens
// //       PQ->compute_density_derivative();
// //       #pragma omp parallel for
// //       for(size_t p = 0; p < PQ->nb(); p++)
// //       {
// //          const Vector P = PQ->pos(p);
// //          const Vector V = PQ->vel(p);
// //          float density = PQ->get_float_attr("density", p);
// //          size_t pindex = PQ->index(P);
// //          size_t i,j,k;
// //          PQ->anti_index(pindex, i,j,k);
// //          float ki = (1.0/dt) * PQ->get_float_attr("density_derivative",p) * PQ->get_float_attr("factor", p);

// //          Vector sum{};
// // //cell
// //          const std::vector<size_t>& cells = PQ->cell_contents(i,j,k); 
// //          for(size_t a=0;a<cells.size();a++)
// //          {
// //             size_t pid = cells[a]; 
// //             float kj = (1.0/dt) * PQ->get_float_attr("density_derivative",pid) * PQ->get_float_attr("factor", pid);
// //             float density_j = PQ->get_float_attr("density", pid);
// //             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
// //          }
// // //neighbors
// //          std::vector<size_t> neighbors;
// //          PQ->neighbor_cells(i,j,k,neighbors);
// //          for(size_t a=0;a<neighbors.size();a++)
// //          {
// //             size_t pid = cells[a]; 
// //             float kj = (1.0/dt) * PQ->get_float_attr("density_derivative",pid) * PQ->get_float_attr("factor", pid);
// //             float density_j = PQ->get_float_attr("density", pid);
// //             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
// //          }

// //          PQ->set_vel(p, V -  (dt * sum));

// //       }
// //       //what to do about new avg??
// //       //PQ->compute_density_derivative();
// //       densDAvg = PQ->average_density_derivative();
// //       iter++;
      
// //    }
// }

// void DFSPHSolver::get_timestep()
// {
//     double pdiameter = PQ->get_radius() / 2.0;
//     float max_vel = PQ->max_velocity();
//     float lambda = 0.4f;
//     double max_dt;
//     if(max_vel != 0)
//        max_dt = lambda * (pdiameter/max_vel);
//     else max_dt = 1.0/48.0;   
//     if(dt > max_dt) dt = max_dt;
//     if(dt > 1.0/48.0) dt = 1.0/48.0;
//     //dt=1.0/48.0;

// }

// void DFSPHSolver::advance_velocity()
// {
// #pragma omp parallel for
//     for( size_t i=0;i<PQ->nb();i++ )
//     {    
//       Vector A = PQ->accel(i);
//       float Amag = A.magnitude();
//       if(Amag > acceleration_clamp)
//       {
//          A *= acceleration_clamp/Amag;
//       }
//       Vector V = PQ->vel(i) + A*dt;
//       float Vmag = V.magnitude();
//       if(Vmag > velocity_clamp)
//       {
//          V *= velocity_clamp/Vmag;
//       }
//       PQ->set_vel( i, V );
//    }
// }

// void DFSPHSolver::advance_position()
// {
// #pragma omp parallel for
//    for( size_t i=0;i<PQ->nb();i++ )
//    {
//       PQ->set_pos( i, PQ->pos(i) + PQ->vel(i)*dt );
//    }
// }

// pba::GISolver pba::CreateDFSPHSolver( SPHState& pq, Force& f, float vclamp, float aclamp, DFSPHForce& sphforce )
// {
//    return GISolver( new DFSPHSolver( pq, f, vclamp, aclamp, sphforce ) );
// }



WCSPHSolverWithCollisions::WCSPHSolverWithCollisions(SPHState& pq, Force& f, double vclamp, double aclamp, ElasticCollisionHandler& coll, GISolver& sol) :
  PQ(pq),
  force(f),
  CS(coll),
  solver(sol),
  velocity_clamp(vclamp),
  acceleration_clamp(aclamp),
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





pba::GISolver pba::CreateWCSPHSolver( SPHState& pq, Force& f, float vclamp, float aclamp, ElasticCollisionHandler& coll, GISolver& sol )
{
   return GISolver( new WCSPHSolverWithCollisions( pq, f, vclamp, aclamp, coll, sol ) );
}
