#include "ParticleEmitter.h"

using namespace pba;

ParticleEmitter::ParticleEmitter() :
    location      (Vector(0,0,0)),
    velocity (Vector(0,0,0)),
    rate(0),
    particle_color(Color(0,0,1,1))
{
}


ParticleEmitter::ParticleEmitter( const Vector& loc, const Vector& vel) :
    location      (loc),
    velocity (vel),
    rate(0),
    particle_color(Color(0,0,1,1))
{
}


void ParticleEmitter::emitCube(DynamicalState& state, int numParticlesPerAxis, const Vector& center)
{



    int numPoints = numParticlesPerAxis * numParticlesPerAxis * numParticlesPerAxis;
    //int numPoints = rate * rate * rate;
    state->add(numPoints);
    int i = state->nb()-numPoints;
    std::cout << "Emit: Total Points " << state->nb() << std::endl;

    //numParticlesPerAxis = rate; 
   
    for (int x = 0; x < numParticlesPerAxis; x++)
    {
        for (int y = 0; y < numParticlesPerAxis; y++)
        { 
            for (int z = 0; z < numParticlesPerAxis; z++)
            {
                float spacing = 2.0f * state->rad(i);

                //make cube - center cube
                float px = x * spacing + center.X() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                float py = y * spacing + center.Y() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                float pz = z * spacing + center.Z() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float px = x * spacing + location.X() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float py = y * spacing + location.Y() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float pz = z * spacing + location.Z() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                state->set_pos(i, Vector(px,py,pz));
                state->set_vel(i, Vector(randf(),randf(),randf()));
                state->set_mass(i,1);
                state->set_ci(i, Color(0,0,1,1));
                state->set_id(i,i);
                i++;
            }
      }
   }


}

void ParticleEmitter::emitCube(SPHState& state, int numParticlesPerAxis, const Vector& center)
{



    int numPoints = numParticlesPerAxis * numParticlesPerAxis * numParticlesPerAxis;
    //int numPoints = rate * rate * rate;
    state->add(numPoints);
    int i = state->nb()-numPoints;
    std::cout << "Emit: Total Points " << state->nb() << std::endl;

    //numParticlesPerAxis = rate; 
   
    for (int x = 0; x < numParticlesPerAxis; x++)
    {
        for (int y = 0; y < numParticlesPerAxis; y++)
        { 
            for (int z = 0; z < numParticlesPerAxis; z++)
            {
                float spacing = 2.0f * state->get_particle_radius();

                //make cube - center cube
                float px = x * spacing + center.X() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                float py = y * spacing + center.Y() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                float pz = z * spacing + center.Z() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float px = x * spacing + location.X() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float py = y * spacing + location.Y() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                // float pz = z * spacing + location.Z() - (numParticlesPerAxis - 1) * spacing * 0.5f;
                state->set_pos(i, Vector(px,py,pz));
                state->set_vel(i, Vector(0,0,0));
                state->set_mass(i, state->get_float_attr("volume", i) * state->get_density0());
                state->set_ci(i, Color(0,0,1,1));
                state->set_id(i,i);
                i++;
            }
      }
   }


}