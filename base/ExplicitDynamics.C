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
  PQ->set_stage(POSITION_UPDATED);

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
  PQ->set_stage(POSITION_UPDATED);

}
AdvancePositionSoftBody::AdvancePositionSoftBody(SoftBodyState& pq) :PQ(pq){}
AdvancePositionSoftBody::~AdvancePositionSoftBody(){}
void AdvancePositionSoftBody::init(){}
void AdvancePositionSoftBody::solve(const double dt)
{
  #pragma omp parallel for
  for(size_t i=0; i<PQ->nb(); i++)
  {
      PQ->set_pos(i, PQ->pos(i) + PQ->vel(i)*dt);
  }
  PQ->set_stage(POSITION_UPDATED);

}
AdvancePositionRigidBody::AdvancePositionRigidBody(RigidBodyState& pq) :PQ(pq){}
AdvancePositionRigidBody::~AdvancePositionRigidBody(){}
void AdvancePositionRigidBody::init(){}
void AdvancePositionRigidBody::solve(const double dt)
{
  #pragma omp parallel for
  for(size_t i=0; i<PQ->nb(); i++)
  {
      PQ->set_pos(i, PQ->pos(i) + PQ->vel(i)*dt);
  }
  PQ->set_stage(POSITION_UPDATED);

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
    PQ->set_stage(POSITION_UPDATED);

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
    PQ->set_stage(POSITION_UPDATED);

}
AdvancePositionWithCollisionsSoftBody::AdvancePositionWithCollisionsSoftBody(SoftBodyState& pq, CollisionHandler& coll):PQ(pq),CH(coll){}
AdvancePositionWithCollisionsSoftBody::~AdvancePositionWithCollisionsSoftBody(){}
void AdvancePositionWithCollisionsSoftBody::init(){}
void AdvancePositionWithCollisionsSoftBody::solve(const double dt)
{
    #pragma omp parallel for
    for( size_t i=0;i<PQ->nb();i++ )
    {
      PQ->set_pos( i, PQ->pos(i) + PQ->vel(i)*dt );
    }
    DynamicalState pq = PQ;
    CH.handle_collisions(dt, pq);
    PQ->set_stage(POSITION_UPDATED);
}

AdvancePositionWithCollisionsRigidBody::AdvancePositionWithCollisionsRigidBody(RigidBodyState& pq, CollisionHandler& coll):PQ(pq),CH(coll){}
AdvancePositionWithCollisionsRigidBody::~AdvancePositionWithCollisionsRigidBody(){}
void AdvancePositionWithCollisionsRigidBody::init(){}
void AdvancePositionWithCollisionsRigidBody::solve(const double dt)
{
    #pragma omp parallel for
    for( size_t i=0;i<PQ->nb();i++ )
    {
      PQ->set_pos( i, PQ->pos(i) + PQ->vel(i)*dt );
    }
    DynamicalState pq = PQ;
    CH.handle_collisions(dt, pq);
    PQ->set_stage(POSITION_UPDATED);
}

AdvanceVelocity::AdvanceVelocity(DynamicalState& pq, Force& f, float a, float v): PQ(pq),force(f),acceleration_clamp(a), velocity_clamp(v){}
AdvanceVelocity::~AdvanceVelocity(){}
void AdvanceVelocity::init(){}
void AdvanceVelocity::solve(const double dt)
{

  force->compute(PQ, dt); // computes the force and stores (force/mass) in the state vector acceleration member
  #pragma omp parallel for
  for( size_t i=0;i<PQ->nb();i++ )
  {    
    Vector A = PQ->accel(i);
    float Amag = A.magnitude();
    if(Amag > acceleration_clamp)
    {
      A *= acceleration_clamp/Amag;
    }
    Vector V = PQ->vel(i) + A*dt;
    float Vmag = V.magnitude();
    if(Vmag > velocity_clamp)
    {
      V *= velocity_clamp/Vmag;
    }
    PQ->set_vel( i, V );
  }
  PQ->set_stage(VELOCITY_UPDATED);

}
AdvanceVelocitySPH::AdvanceVelocitySPH(SPHState& pq, Force& f, float a, float v): PQ(pq),force(f),acceleration_clamp(a), velocity_clamp(v){}
AdvanceVelocitySPH::~AdvanceVelocitySPH(){}
void AdvanceVelocitySPH::init(){}
void AdvanceVelocitySPH::solve(const double dt)
{

  force->compute(PQ, dt); // computes the force and stores (force/mass) in the state vector acceleration member
  #pragma omp parallel for
  for( size_t i=0;i<PQ->nb();i++ )
  {    
    Vector A = PQ->accel(i);
    float Amag = A.magnitude();
    if(Amag > acceleration_clamp)
    {
      A *= acceleration_clamp/Amag;
    }
    Vector V = PQ->vel(i) + A*dt;
    float Vmag = V.magnitude();
    if(Vmag > velocity_clamp)
    {
      V *= velocity_clamp/Vmag;
    }
    PQ->set_vel( i, V );
  }
  PQ->set_stage(VELOCITY_UPDATED);

}
AdvanceVelocitySoftBody::AdvanceVelocitySoftBody(SoftBodyState& pq, Force& f, float a, float v): PQ(pq),force(f),acceleration_clamp(a), velocity_clamp(v){}
AdvanceVelocitySoftBody::~AdvanceVelocitySoftBody(){}
void AdvanceVelocitySoftBody::init(){}
void AdvanceVelocitySoftBody::solve(const double dt)
{
  force->compute(PQ, dt); // computes the force and stores (force/mass) in the state vector acceleration member
  #pragma omp parallel for
  for( size_t i=0;i<PQ->nb();i++ )
  {    
    Vector A = PQ->accel(i);
    float Amag = A.magnitude();
    if(Amag > acceleration_clamp)
    {
      A *= acceleration_clamp/Amag;
    }
    Vector V = PQ->vel(i) + A*dt;
    float Vmag = V.magnitude();
    if(Vmag > velocity_clamp)
    {
      V *= velocity_clamp/Vmag;
    }
    PQ->set_vel( i, V );
  }
  PQ->set_stage(VELOCITY_UPDATED);
}
AdvanceVelocityRigidBody::AdvanceVelocityRigidBody(RigidBodyState& pq, Force& f, float a, float v): PQ(pq),force(f),acceleration_clamp(a), velocity_clamp(v){}
AdvanceVelocityRigidBody::~AdvanceVelocityRigidBody(){}
void AdvanceVelocityRigidBody::init(){}
void AdvanceVelocityRigidBody::solve(const double dt)
{
  force->compute(PQ, dt); // computes the force and stores (force/mass) in the state vector acceleration member
  #pragma omp parallel for
  for( size_t i=0;i<PQ->nb();i++ )
  {    
    Vector A = PQ->accel(i);
    float Amag = A.magnitude();
    if(Amag > acceleration_clamp)
    {
      A *= acceleration_clamp/Amag;
    }
    Vector V = PQ->vel(i) + A*dt;
    float Vmag = V.magnitude();
    if(Vmag > velocity_clamp)
    {
      V *= velocity_clamp/Vmag;
    }
    PQ->set_vel( i, V );
  }
  PQ->set_stage(VELOCITY_UPDATED);
}

GISolver pba::CreateAdvanceVelocity(DynamicalState& pq, Force& f, float a, float v)
{
  return std::make_shared<AdvanceVelocity>(pq, f, a, v);
}
GISolver pba::CreateAdvancePosition(DynamicalState& pq)
{
    return std::make_shared<AdvancePosition>(pq);
}
GISolver pba::CreateAdvancePositionColl(DynamicalState& pq, CollisionHandler& coll )
{
  return std::make_shared<AdvancePositionWithCollisions>(pq, coll);
}

GISolver pba::CreateAdvanceVelocitySPH(SPHState& pq, Force& f, float a, float v)
{
  return std::make_shared<AdvanceVelocitySPH>(pq, f, a, v);
}
GISolver pba::CreateAdvancePositionSPH(SPHState& pq)
{
    return std::make_shared<AdvancePositionSPH>(pq);
}
GISolver pba::CreateAdvancePositionCollSPH(SPHState& pq, CollisionHandler& coll )
{
  return std::make_shared<AdvancePositionWithCollisionsSPH>(pq, coll);

}

GISolver pba::CreateAdvanceVelocitySoftBody(SoftBodyState& pq, Force& f, float a, float v)
{
  return std::make_shared<AdvanceVelocitySoftBody>(pq, f, a, v);
}
GISolver pba::CreateAdvancePositionSoftBody(SoftBodyState& pq)
{
    return std::make_shared<AdvancePositionSoftBody>(pq);
}
GISolver pba::CreateAdvancePositionCollSoftBody(SoftBodyState& pq, CollisionHandler& coll )
{
  return std::make_shared<AdvancePositionWithCollisionsSoftBody>(pq, coll);
}

GISolver pba::CreateAdvanceVelocityRigidBody(RigidBodyState& pq, Force& f, float a, float v)
{
  return std::make_shared<AdvanceVelocityRigidBody>(pq, f, a, v);
}
GISolver pba::CreateAdvancePositionRigidBody(RigidBodyState& pq)
{
    return std::make_shared<AdvancePositionRigidBody>(pq);
}
GISolver pba::CreateAdvancePositionCollRigidBody(RigidBodyState& pq, CollisionHandler& coll )
{
  return std::make_shared<AdvancePositionWithCollisionsRigidBody>(pq, coll);
}