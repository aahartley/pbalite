#ifndef __PBA_FORCE_LIBRARY_H__
#define __PBA_FORCE_LIBRARY_H__

#include "DynamicalState.h"
#include "SPHState.h"
#include <memory>
#include "Force.h"
namespace pba
{
class GravityForce : public ForceBase
{
  public:
    GravityForce(const Vector& g):
      gravity (g)
    {}
    ~GravityForce(){}
    void compute(DynamicalState& s, const double dt);
    void compute(SPHState& s, const double dt);
    void set_gravity(const Vector& g){gravity=g;}
    const Vector& get_gravity() const{return gravity;}
  private:
    Vector gravity; 
};

class TaitPressureForce : public ForceBase
{
  public:
    TaitPressureForce(const float s, const float rest, const float e):
      strength (s),
      rho_0 (rest),
      gamma (e)
    {}
    ~TaitPressureForce(){}
    void compute(DynamicalState& s, const double dt);
    void compute(SPHState& s, const double dt);
    void set_strenth(const float s){strength=s;}
    const float get_strength() const{return strength;}
    void set_rho0(const float s){rho_0=s;}
    const float get_rho0() const{return rho_0;}
    void set_gamma(const float s){gamma=s;}
    const float get_gamma() const{return gamma;}
  private:
    float strength;
    float rho_0;
    float gamma;
};


class HarmonicOscillatorForce : public ForceBase
{
  public:
    HarmonicOscillatorForce(const double& kd) :
      Kd(kd)
    {}
    ~HarmonicOscillatorForce(){}
    void compute(DynamicalState& s, const double dt);
    void compute(SPHState& s, const double dt);
    // set and gets to be able to change the force strength
    void set_kd(const double& v){ Kd = v; }
    const double& get_kd() const { return Kd; }
    //EXAMPLE USAGE:
    //Force harmonic_oscillator = make_harmonic_oscillator_force( 1.0 ); // Input parameter is the force strength kd
    //std::shared_ptr<HarmonicOscillatorForce> f = dynamic_pointer_cast<HarmonicOscillatorForce>(harmonic_oscillator);
    //f->set_kd( f->get_kd()/1.1 ); // reduces the magnitude of the strength
  private:
    double Kd;
};



class AccumulatingForce : public ForceBase
{
  public:
    AccumulatingForce(){}
    ~AccumulatingForce(){}
    void compute(DynamicalState& s, const double dt);
    void compute(SPHState& s, const double dt);
    // Build up the collection of forces to accumulate
    void add_force(Force& f);
  private:
    std::vector<Force> forces;
};

//EXAMPLE:
// Create several forces, in this case: gravity and harmonic oscillator
//Force gravity = make_gravity_force( Vector(0,-9.8,0) );
//Force harmonic_oscillator = make_harmonic_oscillator_force( 1.0 );
// Create for accumulator force
//Force accumulator = make_accumulating_force();
// Add our forces to the accumulator, after the needed cast
//std::shared_ptr<AccumulatingForce> f = dynamic_pointer_cast<AccumulatingForce>(accumulator);
//f->add_force( gravity );
//f->add_force( harmonic_oscillator );
pba::Force CreateGravityForce(const Vector& g);
pba::Force CreateTaitPressureForce(const float s, const float rest, const float g);

pba::Force CreateHarmonicOscillatorForce(const double& k);

pba::Force CreateAccumulatingForce();
}

#endif