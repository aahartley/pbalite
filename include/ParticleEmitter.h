#ifndef __PBA_PARTICLE_EMITTER_H__
#define __PBA_PARTICLE_EMITTER_H__
#include "Vector.h"
#include "Color.h"
#include "DynamicalState.h"
#include "SPHState.h"
#include <random>

namespace pba
{

class ParticleEmitter 
{
  public:
    ParticleEmitter();
    ParticleEmitter( const Vector& loc, const Vector& velocity );
    ~ParticleEmitter(){}

    void emitCube(DynamicalState& state, int numParticlesPerAxis, const Vector& center);
    void emitCube(SPHState& state, int numParticlesPerAxis, const Vector& center);


    float randf() { return ((float)rand() / RAND_MAX) * 2.0f - 1.0f; } // range: [-1, 1]


    int emission_rate() const { return rate; }
    void set_emission_rate(const int r){ rate = r; }

    void set_particle_color(const Color& c) {particle_color = c;}

  private:

    Vector location;
    Vector velocity;
    int rate;
    Color particle_color;
};



}


#endif