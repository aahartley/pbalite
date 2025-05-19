#include "ForceLibrary.h"

using namespace pba;

void GravityForce::compute(DynamicalState& s, const double dt)
{
  #pragma omp parallel for
  for(size_t p=0; p<s->nb(); p++)
  {
      s->set_accel(p, s->accel(p) + gravity);
  }
}
void GravityForce::compute(SPHState& s, const double dt)
{
  #pragma omp parallel for
  for(size_t p=0; p<s->nb(); p++)
  {
      s->set_accel(p, s->accel(p) + gravity);
  }
}
void GravityForce::compute(SoftBodyState& s, const double dt)
{
  #pragma omp parallel for
  for(size_t p=0; p<s->nb(); p++)
  {
      s->set_accel(p, s->accel(p) + gravity);
  }
}

void TaitPressureForce::compute(DynamicalState& s, const double dt)
{
  std::cout << "not fluid\n";
}
void TaitPressureForce::compute(SPHState& s, const double dt)
{
  #pragma omp parallel for
  for(size_t p=0; p<s->nb(); p++)
  {
    Vector pressure;
    const Vector P = s->pos(p);
    std::vector<size_t> neighbors;
    s->neighbors_list(neighbors, P, s->get_neighborParallel());
    float dens_a = s->get_float_attr("density",p);
    dens_a = std::fmax(dens_a, s->get_density0());
    for(size_t a = 0; a < neighbors.size(); a++)
    {
      size_t pid = neighbors[a]; 
      float dens_b = s->get_float_attr("density",pid);
      dens_b = std::fmax(dens_b, s->get_density0());

      float pa = strength * (std::pow(dens_a/rho_0,gamma) - 1.0);
      float pb = strength * (std::pow(dens_b/rho_0,gamma) - 1.0);
      pressure += s->grad_weight(pid, P) * s->mass(pid)
                * ((pa/(dens_a*dens_a)) + (pb/(dens_b*dens_b))) ;    
    }

    s->set_accel(p, s->accel(p) - pressure);
  }
}
void TaitPressureForce::compute(SoftBodyState& s, const double dt)
{
  std::cout << "not fluid\n";
}

void HarmonicOscillatorForce::compute(DynamicalState& s, const double dt)
{
    #pragma omp parallel for
    for(size_t p=0; p<s->nb(); p++)
    {
        s->set_accel(p, s->accel(p) - Kd * s->pos(p) / s->mass(p));
    }
}
void HarmonicOscillatorForce::compute(SPHState& s, const double dt)
{
    #pragma omp parallel for
    for(size_t p=0; p<s->nb(); p++)
    {
        s->set_accel(p, s->accel(p) - Kd * s->pos(p) / s->mass(p));
    }
}
void HarmonicOscillatorForce::compute(SoftBodyState& s, const double dt)
{
    #pragma omp parallel for
    for(size_t p=0; p<s->nb(); p++)
    {
        s->set_accel(p, s->accel(p) - Kd * s->pos(p) / s->mass(p));
    }
}
void AccumulatingForce::add_force(Force& f){ forces.push_back(f); }
void AccumulatingForce::compute(DynamicalState& s, const double dt)
{
    #pragma omp parallel for
    for(size_t p=0; p<s->nb(); p++)
    {
        s->set_accel(p, Vector(0,0,0)); // initialize accelerations to zero before we start accumulating
    }
    for(size_t p=0; p<forces.size(); p++)
    {
        forces[p]->compute(s,dt);
    }
}
void AccumulatingForce::compute(SPHState& s, const double dt)
{
    #pragma omp parallel for
    for(size_t p=0; p<s->nb(); p++)
    {
        s->set_accel(p, Vector(0,0,0)); // initialize accelerations to zero before we start accumulating
    }
    for(size_t p=0; p<forces.size(); p++)
    {
        forces[p]->compute(s,dt);
    }
}
void AccumulatingForce::compute(SoftBodyState& s, const double dt)
{
    #pragma omp parallel for
    for(size_t p=0; p<s->nb(); p++)
    {
        s->set_accel(p, Vector(0,0,0)); // initialize accelerations to zero before we start accumulating
    }
    for(size_t p=0; p<forces.size(); p++)
    {
        forces[p]->compute(s,dt);
    }
}

void AccumulatingStrutForce::compute( DynamicalState& s, const double dt ){}
void AccumulatingStrutForce::compute( SPHState& s, const double dt ){}
void AccumulatingStrutForce::compute( SoftBodyState& s, const double dt )
{
  #pragma omp parallel for
  for( size_t i = 0; i < s->nb_pairs(); i++ )
  {
    const SoftEdge& se = s->get_connected_pair(i);
    const size_t& inode = se->get_first_node();
    const size_t& jnode = se->get_second_node();
    Vector d_ij = s->pos(jnode) - s->pos(inode);
    Vector v_ij = s->vel(jnode)-s->vel(inode);
    float r_mass = (s->mass(inode) * s->mass(jnode)) / (s->mass(inode) + s->mass(jnode));
    if(crit_damp)
    {
    float critical_damping = 2 * std::sqrt(r_mass * spring);
    friction = critical_damping;
    }
    Vector F;
    //both particles same spot
    //d_ij unnomrlaized is full error
    if(se->get_edge_length() < 0.000001 )
    {
        F = d_ij * spring;
        F += v_ij * friction;
    }
    else
    {
      double separation = d_ij.magnitude() - se->get_edge_length();
      if(d_ij.magnitude() > 0) d_ij.normalize();
      F = d_ij * (separation*spring);
      F += d_ij * ( d_ij*v_ij ) * friction;
    }
    #pragma omp critical
    {
    s->set_accel(jnode, s->accel(jnode) - F/s->mass(jnode));
    s->set_accel(inode, s->accel(inode) + F/s->mass(inode));
    }
  }
}

Force pba::CreateGravityForce(const Vector& g)
{
  return std::make_shared<GravityForce>(g);
}
Force pba::CreateTaitPressureForce(const float s, const float rest, const float g)
{
  return std::make_shared<TaitPressureForce>(s, rest, g);
}
Force pba::CreateHarmonicOscillatorForce(const double& k)
{
  return std::make_shared<HarmonicOscillatorForce>(k);
}
Force pba::CreateAccumulatingForce()
{
  return std::make_shared<AccumulatingForce>();
}
Force pba::CreateAccumulatingStrutForce(const double g, const double f, const bool c)
{
  return std::make_shared<AccumulatingStrutForce>(g, f, c);
}
