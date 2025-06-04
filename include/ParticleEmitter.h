//*******************************************************************
//
//   ParticleEmitter.h
//
//   Particle emissions
//
//
//
//*******************************************************************

#ifndef __PBA_PARTICLE_EMITTER_H__
#define __PBA_PARTICLE_EMITTER_H__

#include "Vector.h"
#include "Color.h"
#include "DynamicalState.h"
#include "SPHState.h"

#include <random>

namespace pba
{
/*!
  Particle emission functionality
 */
class ParticleEmitter 
{
  public:
    ParticleEmitter();
    ParticleEmitter(const Vector& loc, const Vector& velocity);
    ~ParticleEmitter(){}

    void EmitCube(DynamicalStateData& state, int numParticlesPerAxis, const Vector& center);
    void EmitCube(SPHStateData& state, int numParticlesPerAxis, const Vector& center);

    float Randf() { return ((float)rand() / RAND_MAX) * 2.0f - 1.0f; } // range: [-1, 1]

    int emission_rate() const { return rate_; }
    void set_emission_rate(int r) { rate_ = r; }
    void set_particle_color(const Color& c) { particle_color_ = c;}

  private:
    Vector location_;
    Vector velocity_;
    int rate_;
    Color particle_color_;
};



} //end of pba namespace


#endif //__PBA_PARTICLE_EMITTER_H__