#include "Viscosity.h"

using namespace pba;

void ExplicitViscosity::compute(DynamicalState& s, const double dt){}

void ExplicitViscosity::compute(SPHState& s, const double dt)
{
    #pragma omp parallel for
    for( size_t p = 0; p<s->nb(); p++)
    {
        float density = s->get_float_attr( "density", p );
       float mass = s->mass(p);
       if(mass==0 ) std::cout << "force bad mass\n";
       if( density ==0) std::cout << "force bad density\n";
       float kinematic_viscosity = dynamic_viscosity / density;
       float h = s->get_radius();
       const int dimensions = 3;
       const Vector& P = s->pos(p);
       const Vector& V = s->vel(p);
       std::vector<size_t> neighbors;
       s->neighbors_list(neighbors, P, s->get_neighborParallel());
 
       Vector laplacian{};
       Vector viscosity_force{}; 
 
 // PARTICLES IN SAME CELL 
       for(size_t a = 0; a < neighbors.size(); a++)
       {
         size_t pid = neighbors[a];
         Vector v_ij = V - s->vel(pid);
         Vector x_ij = P - s->pos(pid);
         float pid_density = s->get_float_attr("density", pid);
         float pid_mass = s->mass(pid);
 
         laplacian += ( (pid_mass / pid_density) * 
                         ( (v_ij * x_ij) / ((x_ij.magnitude() * x_ij.magnitude()) + 0.01 * (h*h)) ) ) * s->grad_weight(pid, P); 
       
       }
 
       laplacian = (2*(dimensions+2) * laplacian);
       viscosity_force = (mass * kinematic_viscosity) * laplacian;
       s->set_accel(p, s->accel(p) + (viscosity_force)/s->mass(p)); // / by mass?
       #pragma omp critical
       if(std::isnan(s->accel(p).X()) || std::isnan(s->accel(p).Y()) || std::isnan(s->accel(p).Z())) std::cout <<"force bad vel nan\n";
 
    } 
}

Force pba::CreateExplicitViscosity(const float kv)
{
    return std::make_shared<ExplicitViscosity>(kv);
}
