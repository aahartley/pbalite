#include "ExplicitDynamics.h"

using namespace pba;
//explicit
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
AdvancePositionSPH::AdvancePositionSPH(SPHState& pq) :PQ(pq){}
AdvancePositionSPH::~AdvancePositionSPH(){}
void AdvancePositionSPH::init(){}
void AdvancePositionSPH::solve(const double dt)
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
AdvancePositionWithCollisionsSPH::AdvancePositionWithCollisionsSPH(SPHState& pq, CollisionHandler& coll):PQ(pq),CH(coll){}
AdvancePositionWithCollisionsSPH::~AdvancePositionWithCollisionsSPH(){}
void AdvancePositionWithCollisionsSPH::init(){}
void AdvancePositionWithCollisionsSPH::solve(const double dt)
{
    #pragma omp parallel for
    for( size_t i=0;i<PQ->nb();i++ )
    {
      PQ->set_pos( i, PQ->pos(i) + PQ->vel(i)*dt );
    }
    DynamicalState pq = PQ;
    CH.handle_collisions(dt, pq);
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
AdvanceVelocitySPH::AdvanceVelocitySPH(SPHState& pq, Force& f): PQ(pq),force(f){}
AdvanceVelocitySPH::~AdvanceVelocitySPH(){}
void AdvanceVelocitySPH::init(){}
void AdvanceVelocitySPH::solve(const double dt)
{

  force->compute(PQ, dt); // computes the force and stores (force/mass) in the state vector acceleration member
  #pragma omp parallel for
  for(size_t i=0; i<PQ->nb(); i++)
  {
      PQ->set_vel(i, PQ->vel(i) + PQ->accel(i)*dt);
  }
}

GISolver pba::CreateAdvanceVelocity(DynamicalState& pq, Force& f)
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
GISolver pba::CreateAdvanceVelocitySPH(SPHState& pq, Force& f)
{
  return std::make_shared<AdvanceVelocitySPH>(pq, f);
}
GISolver pba::CreateAdvancePositionSPH(SPHState& pq)
{
    return std::make_shared<AdvancePositionSPH>(pq);
}
GISolver pba::CreateAdvancePositionCollSPH(SPHState& pq, CollisionHandler& coll )
{
  return std::make_shared<AdvancePositionWithCollisionsSPH>(pq, coll);

}
