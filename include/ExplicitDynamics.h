#ifndef __PBA_INTEGRATORS_H__
#define __PBA_INTEGRATORS_H__
#include "DynamicalState.h"
#include "SPHState.h"
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



class AdvanceVelocity : public GISolverBase
{
  public:
    AdvanceVelocity(DynamicalState& pq, Force& f);

    ~AdvanceVelocity();
    void init();
    void solve(const double dt);

  private:
    DynamicalState PQ; 
    Force force; 
};
class AdvanceVelocitySPH : public GISolverBase
{
  public:
    AdvanceVelocitySPH(SPHState& pq, Force& f);

    ~AdvanceVelocitySPH();
    void init();
    void solve(const double dt);

  private:
    SPHState PQ; 
    Force force; 
};





GISolver CreateAdvanceVelocity(DynamicalState& pq, Force& f);
GISolver CreateAdvanceVelocitySPH(SPHState& pq, Force& f);

GISolver CreateAdvancePositionColl(DynamicalState& pq, pba::CollisionHandler& coll );
GISolver CreateAdvancePositionCollSPH(SPHState& pq, pba::CollisionHandler& coll );

GISolver CreateAdvancePosition(DynamicalState& pq);
GISolver CreateAdvancePositionSPH(SPHState& pq);

}


#endif
