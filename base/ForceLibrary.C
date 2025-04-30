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

void HarmonicOscillatorForce::compute(DynamicalState& s, const double dt)
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

Force pba::CreateGravityForce(const Vector& g)
{
  return std::make_shared<GravityForce>(g);
}
Force pba::CreateHarmonicOscillatorForce(const double& k)
{
  return std::make_shared<HarmonicOscillatorForce>(k);
}
Force pba::CreateAccumulatingForce()
{
  return std::make_shared<AccumulatingForce>();
}
