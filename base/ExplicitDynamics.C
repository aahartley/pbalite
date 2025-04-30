#include "ExplicitDynamics.h"

using namespace pba;

AdvancePosition::AdvancePosition(DynamicalState& pq) :PQ(pq){}
AdvancePosition::~AdvancePosition(){}
void AdvancePosition::init(){}
void AdvancePosition::solve(const double dt)
{

  #pragma omp parallel for
  for(size_t i=0; i<PQ->nb(); i++)
  {
      PQ->set_pos(i, PQ->pos(i) + PQ->vel(i)*dt);
  }

}

AdvancePositionWithCollisions::AdvancePositionWithCollisions(DynamicalState& pq, CollisionHandler& coll):PQ(pq),CH(coll){}
AdvancePositionWithCollisions::~AdvancePositionWithCollisions(){}
void AdvancePositionWithCollisions::init(){}
void AdvancePositionWithCollisions::solve(const double dt)
{
    #pragma omp parallel for
    for( size_t i=0;i<PQ->nb();i++ )
    {
      PQ->set_pos( i, PQ->pos(i) + PQ->vel(i)*dt );
    }
    CH.handle_collisions(dt, PQ);
}

AdvanceVelocity::AdvanceVelocity(DynamicalState& pq, Force& f): PQ(pq),force(f){}
AdvanceVelocity::~AdvanceVelocity(){}
void AdvanceVelocity::init(){}
void AdvanceVelocity::solve(const double dt)
{

  force->compute(PQ, dt); // computes the force and stores (force/mass) in the state vector acceleration member
  #pragma omp parallel for
  for(size_t i=0; i<PQ->nb(); i++)
  {
      PQ->set_vel(i, PQ->vel(i) + PQ->accel(i)*dt);
  }
}

GISolver pba::CreateAdvanceVelcity(DynamicalState& pq, Force& f)
{
  return std::make_shared<AdvanceVelocity>(pq, f);
}
GISolver pba::CreateAdvancePosition(DynamicalState& pq)
{
    return std::make_shared<AdvancePosition>(pq);
}
GISolver pba::CreateAdvancePositionColl(DynamicalState& pq, CollisionHandler& coll )
{
  return std::make_shared<AdvancePositionWithCollisions>(pq, coll);

}
