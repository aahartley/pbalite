//*******************************************************************
//
//   DFSPHSolver.C
//
//   DFSPH Solver with CFL condition
//   Two Pressure solvers
//
//
//
//*******************************************************************
//  DFSPH functions are based off the 2015 paper: https://animation.rwth-aachen.de/media/papers/2015-SCA-DFSPH.pdf
//  2017 paper, where the equations are written differently (still equivalent): https://animation.rwth-aachen.de/media/papers/2017-TVCG-ViscousDFSPH.pdf

#include "DFSPHSolver.h"

#include <chrono>

namespace pba
{

//sim loop does this already --- no need call
//requirements for sim to start
void DFSPHSolver::Init()
{
   pq_.Populate();
   pq_.ComputeDensity();
   pq_.ComputeFactor();
}

//naive collisions to keep particles in box
void DFSPHSolver::fakecs()
{
   float length = 19 * (0.025*2);
   float botx = -(length/2)*3;
   float WIDTHx = std::abs(botx) ;

   float boty = -length;
   float WIDTHy = std::abs(boty) ;

   float botz = -(length/2);
   float WIDTHz = std::abs(botz) ;

   float radius = pq_.get_particle_radius();
   float coef = 0.9;
   #pragma omp parallel for
   for (size_t i = 0; i < pq_.nb(); ++i)
   {
      if (pq_.pos(i).X() >= WIDTHx-radius)
      { 
         pq_.set_pos(i, Vector((WIDTHx-radius), pq_.pos(i).Y(), pq_.pos(i).Z()));
         pq_.set_vel(i, Vector((pq_.vel(i).X() * -coef),pq_.vel(i).Y(),pq_.vel(i).Z()));
         //pq_.set_vel(i, Vector(0,pq_.vel(i).Y(),pq_.vel(i).Z()));

      }
      if (pq_.pos(i).X() <= botx+radius)
      {
         pq_.set_pos(i, Vector((botx+radius), pq_.pos(i).Y(), pq_.pos(i).Z()));
         pq_.set_vel(i, Vector((pq_.vel(i).X() * -coef),pq_.vel(i).Y(),pq_.vel(i).Z()));
         //pq_.set_vel(i, Vector(0,pq_.vel(i).Y(),pq_.vel(i).Z()));

      }
      if (pq_.pos(i).Y() >= WIDTHy-radius)
      { 
         pq_.set_pos(i, Vector(pq_.pos(i).X(), (WIDTHy-radius), pq_.pos(i).Z()));
         pq_.set_vel(i, Vector(pq_.vel(i).X(),(pq_.vel(i).Y() * -coef),pq_.vel(i).Z()));
         //pq_.set_vel(i, Vector(pq_.vel(i).X(),0,pq_.vel(i).Z()));

      }
      if (pq_.pos(i).Y() <= boty+radius)
      {
         //std::cout << "botplane hit\n";
         pq_.set_pos(i, Vector(pq_.pos(i).X(), (boty+radius), pq_.pos(i).Z()));
         pq_.set_vel(i, Vector(pq_.vel(i).X(),(pq_.vel(i).Y() * -coef),pq_.vel(i).Z()));
         //pq_.set_vel(i, Vector(pq_.vel(i).X(),0,pq_.vel(i).Z()));

      }
      if (pq_.pos(i).Z() >= WIDTHz-radius)
      { 
         pq_.set_pos(i, Vector(pq_.pos(i).X(), pq_.pos(i).Y(),(WIDTHz-radius)));
         pq_.set_vel(i, Vector(pq_.vel(i).X(),pq_.vel(i).Y(),(pq_.vel(i).Z() * -coef)));
         //pq_.set_vel(i, Vector(pq_.vel(i).X(),pq_.vel(i).Y(),0));

      }
      if(pq_.pos(i).Z() <= botz+radius)
      {
         pq_.set_pos(i, Vector(pq_.pos(i).X(), pq_.pos(i).Y(),(botz+radius)));
         pq_.set_vel(i, Vector(pq_.vel(i).X(),pq_.vel(i).Y(),(pq_.vel(i).Z() * -coef)));
         //pq_.set_vel(i, Vector(pq_.vel(i).X(),pq_.vel(i).Y(),0));
      }
   }
}

void DFSPHSolver::Solve(const double userdt)
{
   user_dt_ = userdt;

   //Occupancy Grid
   pq_.Populate(); 

   pq_.ComputeDensity();

   //Equ (8)
   pq_.ComputeFactor();

   //PPE 2 (Alg. 2): Dρi/Dt = 0
   correct_divergence_error();

   //Non-pressure forces (viscosity, gravity)
   //TODO: surface tension?
   force_.Compute(pq_, dt_);

   //CFL Condition
   get_timestep();

   // ___ Euler integration
   advance_velocity();

   //correct_divergence_error();
   // PPE 1 (Alg. 3): ρ∗i - ρ0 = 0
   correct_density_error();

   // ___ Euler integration
   advance_position();

   //Boundary Handling
   //fakecs();

   //auto start2 = std::chrono::high_resolution_clock::now();

   // auto end2 = std::chrono::high_resolution_clock::now();
   // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
   // std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;

   coll_handler_->HandleCollisions(dt_, pq_);
   //std::cout << "collissions good\n";

}

// PPE 1 (Alg. 3): ρ∗i - ρ0 = 0
void DFSPHSolver::correct_density_error()
{
   #pragma omp parallel for
   for (size_t p = 0; p < pq_.nb(); ++p)
   {
      pq_.ComputePredictedDensity(p, dt_);

      //ρ∗_i - ρ0
      const float source_i =  pq_.get_float_attr("predicted_density", p) - pq_.get_density0(); 

      //clamp source -- if source is neg == expansion, only solve for compression
      const float residuum = std::max(source_i, 0.0f); 

      //Under Equation (9)
      //pressure stifness parameter ki = 1/∆t^2 (ρ∗i - ρ0) αi 
      pq_.set_attr("k_i", p, residuum * (pq_.get_float_attr("factor",p)* (1.0f/(dt_*dt_)))); 

   }

   int iter = 0;
   float average_density_error = 0.0f;
   bool check = false;

   while ((!check || iter < 2) && iter < pq_.get_max_iter() && pq_.nb() != 0)   
   {
      check = true;
      average_density_error = 0.0f;
      //Jacobi Iteration to solve for ki
      density_solve_iteration(average_density_error);

      //0.01 == m_maxError,  m_maxError is given as a percent. p0*maxerror = 0.01% * density0
      const float eta = pq_.get_density0() * pq_.get_m_max_error() * 0.01; 
      check = check && (average_density_error <= eta);
      iter++;
   }  

   std::cout << "DE dens error: " <<average_density_error << "iter: " << iter << '\n';

   #pragma omp parallel for
   for (size_t p = 0; p < pq_.nb(); ++p)
   {
      //End of section 3.2 is Fpi,total = -mi ∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij
      //this is acc =  -∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij
      compute_pressure_acc(p, "density"); 
      const Vector& vel = pq_.vel(p);
      //pressure_acc has negative sign accounted for, unlike psuedo-code in paper
      pq_.set_vel(p, vel + dt_ * pq_.get_vector_attr("pressure_acc", p));
      //#pragma omp critical
      //{
         // if(std::isnan(pq_.vel(p).X()) || std::isinf(pq_.vel(p).X())) std::cout << "dense pq_.vel(p).X() x bad" << "P: "<< p <<'\n';
         // if(std::isnan(pq_.vel(p).Y()) || std::isinf(pq_.vel(p).Y())) std::cout << "dense pq_.vel(p).Y() y bad" << "P: "<< p <<'\n';
         // if(std::isnan(pq_.vel(p).Z()) || std::isinf(pq_.vel(p).Z())) std::cout << "dense pq_.vel(p).Y() z bad" << "P: "<< p <<'\n';
      //}
   }

}

//Jacobi Iteration to solve for ki
void DFSPHSolver::density_solve_iteration(float& average_density_error)
{
   float density_error = 0.0f;

   #pragma omp parallel for
   for (size_t p = 0; p < pq_.nb(); ++p)
   {
      compute_pressure_acc(p, "density"); 
   }   
   #pragma omp  parallel for reduction(+:density_error)
   for (size_t p = 0; p < pq_.nb(); ++p)
   {
      float error_force = compute_error_force(p, "density");
   
      //const float predicted_density = pq_.get_float_attr("predicted_density", p);
 
      const float s_i =  pq_.get_float_attr("predicted_density", p) - pq_.get_density0(); 

      float k_i = pq_.get_float_attr("k_i", p);
  
      //Equation (9)
      //(ρ∗i - ρ0) = 1/∆t^2  ∑_j m_j( f_Pi / m_i - f_Pj / m_i)∇Wij (in my derivations it should be -(1/∆t^2), typo in paper?)
      //(ρ∗i - ρ0) = -1/∆t^2  ∑_j m_j( f_Pi / m_i - f_Pj / m_i)∇Wij
      //(ρ∗i - ρ0) + 1/∆t^2  ∑_j m_j( f_Pi / m_i - f_Pj / m_i)∇Wij = error
      const float residuum = std::max(s_i + error_force, 0.0f); 

      //adjust ki by half of the error during jacobi, clamp to 0 for no negative pressure
      //error is not clamped here, to fix neighborhood pressure even if particle itself is not compressed
      k_i = std::max(k_i + 0.5f * (s_i + error_force) * (pq_.get_float_attr("factor",p) * (1.0f/(dt_*dt_))), 0.0f);

      pq_.set_attr("k_i", p, k_i);

      density_error +=  residuum;

   }

   average_density_error = density_error / pq_.nb();
   // std::cout << "dens error: " <<average_density_error << '\n';
}

//Dρi/Dt = 0
void DFSPHSolver::correct_divergence_error()
{

   #pragma omp parallel for
   for (size_t p = 0; p < pq_.nb(); ++p)
   {
      pq_.ComputeDensityDerivative(p); 
      float density_derivative = pq_.get_float_attr("density_derivative", p);
      density_derivative = std::max(density_derivative, 0.0f); //divergence free so we want density to never be negative? 
      if (pq_.get_dd_clamp())
         pq_.set_attr("density_derivative", p, density_derivative);

      int num_neighbors = 0;
      const Vector& pos = pq_.pos(p);
      std::vector<size_t> neighbors;
      pq_.NeighborsList(neighbors, pos, pq_.get_neighbor_parallel());
      num_neighbors += static_cast<int>(neighbors.size()); //i dont think they are counting the particle itself -1 or not to -1
      //do not solve for particles without enough influence
      if (num_neighbors < 20)
      {
         density_derivative = 0.0f;
         if (pq_.get_dd_clamp())
            pq_.set_attr("density_derivative", p, density_derivative);
      }
      float kv_i = density_derivative * (pq_.get_float_attr("factor", p) * (1.0f/dt_)); //pressure value used to compute acc
      pq_.set_attr("kv_i", p, kv_i);

   }

   int iter = 0;

   float average_density_error = 0.0f;
   bool check = false;

   while ((!check || iter < 1) && iter < pq_.get_max_iter() && pq_.nb() != 0)   
   {
      check = true;
      
      average_density_error = 0.0f;
      divergence_solve_iteration(average_density_error);
      //this jsut equals 1/dt?? 1/0.005 is 200
      //we want density error dp/dt - forces*dt = 0 to be less than the rate of change of time * a fraction of the rest denstiy
      //0.1% * 1000=1
      float eta = (1.0f / dt_) * pq_.get_max_error() * 0.01f * pq_.get_density0();  // maxError is given in percent, 0.1 == max error
      check = check && (average_density_error <= eta);
      
      iter++;
   }
   std::cout << "DivE dens error: " <<average_density_error << "iter: " << iter << '\n';

   #pragma omp parallel for 
   for(size_t p = 0; p < pq_.nb(); p++)
   {
      //End of section 3.2 is Fpi,total = -mi ∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij
      //this is acc =  -∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij
      compute_pressure_acc(p, "divergence"); //has negative accoutned for
      const Vector& vel = pq_.vel(p);
      if(std::isnan(vel.X()) || std::isnan(vel.Y()) || std::isnan(vel.Z())) std::cout <<"divergence bad vel nan\n";

      pq_.set_vel(p, vel + dt_ * pq_.get_vector_attr("pressure_acc", p));
      #pragma omp critical
      {
         if(std::isnan(pq_.vel(p).X()) || std::isinf(pq_.vel(p).X())) std::cout << "dive pq_.vel(p).X() x bad\n";
         if(std::isnan(pq_.vel(p).Y()) || std::isinf(pq_.vel(p).Y())) std::cout << "dive pq_.vel(p).Y() x bad\n";
         if(std::isnan(pq_.vel(p).Z()) || std::isinf(pq_.vel(p).Z())) std::cout << "dive pq_.vel(p).Z() x bad\n";
      }
   }

}

void DFSPHSolver::divergence_solve_iteration(float& average_density_error)
{
   float density_error = 0.0f;

   #pragma omp parallel for
   for (size_t p = 0; p < pq_.nb(); ++p)
   {
      compute_pressure_acc(p, "divergence"); 
   }
   #pragma omp parallel for reduction(+:density_error)
   for (size_t p = 0; p < pq_.nb(); ++p)
   {
      float error_force = compute_error_force(p, "divergence");
   
      const float density_derivative = pq_.get_float_attr("density_derivative", p);
      const float s_i = density_derivative;//source

      float kv_i = pq_.get_float_attr("kv_i", p);

      //(dp/dt) + 1/∆t  ∑_j m_j( f_Pi / m_i - f_Pj / m_i)∇Wij = error
      //should dd be clamped?
      float residuum = std::max(s_i + error_force, 0.0f); 

      int num_neighbors = 0;
      const Vector& pos = pq_.pos(p);
      std::vector<size_t> neighbors;
      pq_.NeighborsList(neighbors, pos, pq_.get_neighbor_parallel());
      num_neighbors += static_cast<int>(neighbors.size()); //i dont think they are counting the particle itself
      // #pragma omp critical
      // std::cout << "neighbors " << num_neighbors <<'\n';
      if (num_neighbors < 20)
         residuum = 0.0f;
      
      //adjust ki by half of the error during jacobi, clamp to 0 for no negative pressure
      //error is not clamped here, to fix neighborhood pressure even if particle itself is not compressed
      kv_i = std::max(kv_i + 0.5f*(s_i + error_force) * (pq_.get_float_attr("factor",p)*(1.0f/dt_)), 0.0f); 

      pq_.set_attr("kv_i", p, kv_i);

      density_error += residuum; 
   }

   average_density_error = density_error / pq_.nb();
}

//End of section 3.2 is Fpi,total = -mi ∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij
//this is acc =  -∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij
void DFSPHSolver::compute_pressure_acc(const size_t p, const std::string& type)
{
   Vector pressure_acci_i(0.f, 0.f, 0.f);
   float ki;
   if (type == "divergence")
      ki = pq_.get_float_attr("kv_i", p);
   else if (type == "density")
      ki = pq_.get_float_attr("k_i", p);
   else
   {
      std::cout <<"ERROR_compute_pressure_acc\n";
      ki = 0;
   }

   const Vector& pos = pq_.pos(p);
   std::vector<size_t> neighbors;
   pq_.NeighborsList(neighbors, pos, pq_.get_neighbor_parallel());

   for (size_t a = 0; a < neighbors.size(); ++a)
   {
      size_t pid = neighbors[a]; 

      float kj;
      if (type == "divergence")
         kj = pq_.get_float_attr("kv_i", pid);
      else if (type == "density")
         kj = pq_.get_float_attr("k_i", pid);
      else kj = 0;
      float psum = (ki) + (kj);

      if(fabs(psum) > pq_.get_meps())//m_eps = (1.0e-5); b/c we set pressure to 0 for certain scenarios
      {
         Vector grad_pj = pq_.get_float_attr("volume", pid) * pq_.GradWeight(pid, pos)  * pq_.get_density0();
         pressure_acci_i += psum * grad_pj;
      }
   }

   pq_.set_attr("pressure_acc", p, -pressure_acci_i );
   #pragma omp critical
   {
      if(std::isnan(pressure_acci_i.X()) || std::isinf(pressure_acci_i.X())) std::cout << "pressure_acci_i x bad\n";
      if(std::isnan(pressure_acci_i.Y()) || std::isinf(pressure_acci_i.Y())) std::cout << "pressure_acci_i y bad\n";
      if(std::isnan(pressure_acci_i.Z()) || std::isinf(pressure_acci_i.Z())) std::cout << "pressure_acci_i z bad\n";
      // if(pressure_acci_i.X() < 0)  std::cout << "pressure_acci_i x <=0 " << pressure_acci_i.X() << '\n';
      // if(pressure_acci_i.Y() < 0) std::cout << "pressure_acci_i y <=0 " << pressure_acci_i.Y() << '\n';
      // if(pressure_acci_i.Z() < 0) std::cout << "pressure_acci_i z <=0 " << pressure_acci_i.Z() << '\n';
   }
}

//RHS Equation (9) -- forces to correct density error (ρ∗i - ρ0) (change in density due to pressure acc)
//∑_j m_j( f_Pi / m_i - f_Pj / m_i)∇Wij
//which is ∑_j m_j( Acc_Pi - Acc_Pj)∇Wij
float DFSPHSolver::compute_error_force(const size_t p, const std::string& type)
{
   float force = 0.f;
   const Vector& pos = pq_.pos(p);
   const Vector& pressure_acc_i = pq_.get_vector_attr("pressure_acc", p);

   std::vector<size_t> neighbors;
   pq_.NeighborsList(neighbors, pos, pq_.get_neighbor_parallel());

   for (size_t a = 0; a < neighbors.size(); ++a)
   {
      size_t pid = neighbors[a]; 

      const Vector& pressure_acc_j = pq_.get_vector_attr("pressure_acc", pid);
      force += pq_.mass(pid) * ((pressure_acc_i) - (pressure_acc_j)) * pq_.GradWeight(pid, pos); 
   }

   if (std::isnan(force) || std::isinf(force)) std::cout << "force bad\n";
   // if(force > 0 || force < 0)
   // {
   //    #pragma omp critical
   //    std::cout << "force: " << force << '\n';
   // }
   if (type == "divergence") force *= dt_;
   else if (type == "density") force *= dt_*dt_;
   else std::cout << "error in error_force\n";
   return force;
}

//CFL condition
void DFSPHSolver::get_timestep()
{
   if (pq_.get_use_user_dt())
   {
      dt_ = user_dt_;
      std::cout << "Dt: " << dt_ << '\n';
      return;
   }
   double pdiameter = pq_.get_radius() / 2.0; //particle radius is kernel_radius /4
   float max_vel = pq_.MaxVelocity();
   float lambda = 0.4f;
   double max_dt;
   if (max_vel != 0)
      max_dt = lambda * (pdiameter/max_vel);
   else max_dt = user_dt_;   

   dt_ = max_dt;
   if (dt_ > user_dt_ || dt_==0) dt_ = user_dt_; //user_dt dictactes max dt from cfl, change?
   //return (dt < 0.0001) ? 0.0001 : dt;

   //hard coding to clamp between 0.005 and 0.0001, can remove
   //if(dt > 0.005) dt = 0.005;
   if (dt_ < 0.0001) dt_ = 0.0001;

   std::cout << "Dt: " << dt_ << '\n';

}

void DFSPHSolver::advance_velocity()
{
   #pragma omp parallel for
   for (size_t i = 0; i <pq_.nb(); ++i)
   {    
      Vector acc = pq_.accel(i);
      float a_mag = acc.magnitude();
      if(a_mag > acceleration_clamp_)
      {
         acc *= acceleration_clamp_/a_mag;
      }
      Vector vel = pq_.vel(i) + acc*dt_;
      float v_mag = vel.magnitude();
      if(v_mag > velocity_clamp_)
      {
         vel *= velocity_clamp_/v_mag;
      }
      pq_.set_vel(i, vel);
      //if(std::isnan(V.X()) || std::isnan(V.Y()) || std::isnan(V.Z())) std::cout <<"vel bad nan\n";

      // #pragma omp critical
      // std::cout << "particle:"<< i <<" vel: \t" << pq_.vel(i).X() << ' ' << pq_.vel(i).Y() << ' ' << pq_.vel(i).Z() << '\n';
   }
}

void DFSPHSolver::advance_position()
{
   #pragma omp parallel for
   for (size_t i = 0; i < pq_.nb(); ++i)
   {
      pq_.set_pos(i, pq_.pos(i) + pq_.vel(i)*dt_);
      // #pragma omp critical
      // std::cout << "particle:"<< i <<" pos: \t" << pq_.pos(i).X() << ' ' << pq_.pos(i).Y() << ' ' << pq_.pos(i).Z() << '\n';
   }
}


} //end of pba namespace