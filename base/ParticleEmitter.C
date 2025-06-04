//*******************************************************************
//
//   ParticleEmitter.C
//
//   Particle emissions
//
//
//
//*******************************************************************

#include "ParticleEmitter.h"
#include <iostream>

namespace pba
{

ParticleEmitter::ParticleEmitter() 
  : location_(Vector(0,0,0)),
    velocity_(Vector(0,0,0)),
    rate_(0),
    particle_color_(Color(0,0,1,1)) {}


ParticleEmitter::ParticleEmitter(const Vector& location, const Vector& vel) 
  : location_(location),
    velocity_(vel),
    rate_(0),
    particle_color_(Color(0,0,1,1)) {}


void ParticleEmitter::EmitCube(DynamicalStateData& state, const int numParticlesPerAxis, const Vector& center)
{
    int numPoints = numParticlesPerAxis * numParticlesPerAxis * numParticlesPerAxis;
    //int numPoints = rate * rate * rate;
    state.Add(numPoints);
    int i = state.nb()-numPoints;
    std::cout << "Emit: Total Points " << state.nb() << std::endl;

    //numParticlesPerAxis = rate; 
   
    for (int x = 0; x < numParticlesPerAxis; ++x)
    {
        for (int y = 0; y < numParticlesPerAxis; ++y)
        { 
            for (int z = 0; z < numParticlesPerAxis; ++z)
            {
                float spacing = 2.0f * state.rad(i);

                //make cube - center cube
                float px = x * spacing + center.X() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                float py = y * spacing + center.Y() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                float pz = z * spacing + center.Z() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float px = x * spacing + location.X() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float py = y * spacing + location.Y() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float pz = z * spacing + location.Z() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                state.set_pos(i, Vector(px, py, pz));
                state.set_vel(i, Vector(Randf(), Randf(), Randf()));
                state.set_mass(i, 1);
                state.set_ci(i, Color(0,0,1,1));
                state.set_id(i, i);
                i++;
            }
      }
   }


}

void ParticleEmitter::EmitCube(SPHStateData& state, const int numParticlesPerAxis, const Vector& center)
{
    int numPoints = numParticlesPerAxis * numParticlesPerAxis * numParticlesPerAxis;
    //int numPoints = rate * rate * rate;
    state.Add(numPoints);
    int i = state.nb()-numPoints;
    std::cout << "Emit: Total Points " << state.nb() << std::endl;

    //numParticlesPerAxis = rate; 
   
    for (int x = 0; x < numParticlesPerAxis; ++x)
    {
        for (int y = 0; y < numParticlesPerAxis; ++y)
        { 
            for (int z = 0; z < numParticlesPerAxis; ++z)
            {
                float spacing = 2.0f * state.get_particle_radius();

                //make cube - center cube
                float px = x * spacing + center.X() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                float py = y * spacing + center.Y() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                float pz = z * spacing + center.Z() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float px = x * spacing + location.X() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float py = y * spacing + location.Y() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float pz = z * spacing + location.Z() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                state.set_pos(i, Vector(px, py, pz));
                state.set_vel(i, Vector(0, 0, 0));
                state.set_mass(i, state.get_float_attr("volume", i) * state.get_density0());
                state.set_ci(i, Color(0,0,1,1));
                state.set_id(i, i);
                i++;
            }
      }
   }


}

}