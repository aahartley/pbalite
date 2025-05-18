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
