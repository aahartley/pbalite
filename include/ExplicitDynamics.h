#ifndef __PBA_INTEGRATORS_H__
#define __PBA_INTEGRATORS_H__
#include "DynamicalState.h"
#include "SPHState.h"
#include "SoftBodyState.h"
#include "RigidBodyState.h"
#include "ForceLibrary.h"
#include "CollisionHandler.h"
#include "GISolver.h"

namespace pba
{
class AdvancePosition : public GISolverBase
{
  public:
    AdvancePosition(DynamicalState& pq);
    ~AdvancePosition();
    void init();
    void solve(const double dt);
  private:
    DynamicalState PQ; 
};

class AdvancePositionSPH : public GISolverBase
{
  public:
    AdvancePositionSPH(SPHState& pq);
    ~AdvancePositionSPH();
    void init();
    void solve(const double dt);
  private:
    SPHState PQ; 
};
class AdvancePositionSoftBody : public GISolverBase
{
  public:
    AdvancePositionSoftBody(SoftBodyState& pq);
    ~AdvancePositionSoftBody();
    void init();
    void solve(const double dt);
  private:
    SoftBodyState PQ; 
};
class AdvancePositionRigidBody : public GISolverBase
{
  public:
    AdvancePositionRigidBody(RigidBodyState& pq);
    ~AdvancePositionRigidBody();
    void init();
    void solve(const double dt);
  private:
    RigidBodyState PQ; 
};

class AdvancePositionWithCollisions : public GISolverBase 
{
  public:
 
    AdvancePositionWithCollisions(DynamicalState& pq, CollisionHandler& coll);
    ~AdvancePositionWithCollisions();
 
    void init();
    void solve(const double dt);
    
  private:
    DynamicalState PQ;
    CollisionHandler& CH;
 
};
class AdvancePositionWithCollisionsSPH : public GISolverBase 
{
  public:
 
    AdvancePositionWithCollisionsSPH(SPHState& pq, CollisionHandler& coll);
    ~AdvancePositionWithCollisionsSPH();
    
    void init();
    void solve(const double dt);
    
  private:
    SPHState PQ;
    CollisionHandler& CH;
 
};
class AdvancePositionWithCollisionsSoftBody : public GISolverBase 
{
  public:
 
    AdvancePositionWithCollisionsSoftBody(SoftBodyState& pq, CollisionHandler& coll);
    ~AdvancePositionWithCollisionsSoftBody();
    
    void init();
    void solve(const double dt);
    
  private:
    SoftBodyState PQ;
    CollisionHandler& CH;
 
};
class AdvancePositionWithCollisionsRigidBody : public GISolverBase 
{
  public:
 
    AdvancePositionWithCollisionsRigidBody(RigidBodyState& pq, CollisionHandler& coll);
    ~AdvancePositionWithCollisionsRigidBody();
    
    void init();
    void solve(const double dt);
    
  private:
    RigidBodyState PQ;
    CollisionHandler& CH;
 
};



class AdvanceVelocity : public GISolverBase
{
  public:
    AdvanceVelocity(DynamicalState& pq, Force& f, float a, float v);
    ~AdvanceVelocity();
    void init();
    void solve(const double dt);
    const float get_velocity_clamp() const { return velocity_clamp; }
    void set_velocity_clamp(const float& v ) { velocity_clamp = v; }

    const float get_acceleration_clamp() const { return acceleration_clamp; }
    void set_acceleration_clamp(const float& v ) { acceleration_clamp = v; }
  private:
    DynamicalState PQ; 
    Force force; 
    float acceleration_clamp;
    float velocity_clamp;

};
class AdvanceVelocitySPH : public GISolverBase
{
  public:
    AdvanceVelocitySPH(SPHState& pq, Force& f, float a, float v);
    ~AdvanceVelocitySPH();
    void init();
    void solve(const double dt);
    const float get_velocity_clamp() const { return velocity_clamp; }
    void set_velocity_clamp(const float& v ) { velocity_clamp = v; }

    const float get_acceleration_clamp() const { return acceleration_clamp; }
    void set_acceleration_clamp(const float& v ) { acceleration_clamp = v; }
  private:
    SPHState PQ; 
    Force force; 
    float acceleration_clamp;
    float velocity_clamp;
};
class AdvanceVelocitySoftBody : public GISolverBase
{
  public:
    AdvanceVelocitySoftBody(SoftBodyState& pq, Force& f, float a, float v);
    ~AdvanceVelocitySoftBody();
    void init();
    void solve(const double dt);
    const float get_velocity_clamp() const { return velocity_clamp; }
    void set_velocity_clamp(const float& v ) { velocity_clamp = v; }

    const float get_acceleration_clamp() const { return acceleration_clamp; }
    void set_acceleration_clamp(const float& v ) { acceleration_clamp = v; }
  private:
    SoftBodyState PQ; 
    Force force; 
    float acceleration_clamp;
    float velocity_clamp;
};
class AdvanceVelocityRigidBody : public GISolverBase
{
  public:
    AdvanceVelocityRigidBody(RigidBodyState& pq, Force& f, float a, float v);

    ~AdvanceVelocityRigidBody();
    void init();
    void solve(const double dt);
    const float get_velocity_clamp() const { return velocity_clamp; }
    void set_velocity_clamp(const float& v ) { velocity_clamp = v; }

    const float get_acceleration_clamp() const { return acceleration_clamp; }
    void set_acceleration_clamp(const float& v ) { acceleration_clamp = v; }
  private:
    RigidBodyState PQ; 
    Force force; 
    float acceleration_clamp;
    float velocity_clamp;
};



GISolver CreateAdvanceVelocity(DynamicalState& pq, Force& f, float a, float v);
GISolver CreateAdvanceVelocitySPH(SPHState& pq, Force& f, float a, float v);
GISolver CreateAdvanceVelocitySoftBody(SoftBodyState& pq, Force& f, float a, float v);
GISolver CreateAdvanceVelocityRigidBody(RigidBodyState& pq, Force& f, float a, float v);

GISolver CreateAdvancePositionColl(DynamicalState& pq, pba::CollisionHandler& coll );
GISolver CreateAdvancePositionCollSPH(SPHState& pq, pba::CollisionHandler& coll );
GISolver CreateAdvancePositionCollSoftBody(SoftBodyState& pq, pba::CollisionHandler& coll );
GISolver CreateAdvancePositionCollRigidBody(RigidBodyState& pq, pba::CollisionHandler& coll );

GISolver CreateAdvancePosition(DynamicalState& pq);
GISolver CreateAdvancePositionSPH(SPHState& pq);
GISolver CreateAdvancePositionSoftBody(SoftBodyState& pq);
GISolver CreateAdvancePositionRigidBody(RigidBodyState& pq);

}


#endif
