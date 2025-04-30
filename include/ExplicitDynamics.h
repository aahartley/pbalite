#ifndef __PBA_INTEGRATORS_H__
#define __PBA_INTEGRATORS_H__
#include "DynamicalState.h"
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





GISolver CreateAdvanceVelcity(DynamicalState& pq, Force& f);
GISolver CreateAdvancePositionColl(DynamicalState& pq, pba::CollisionHandler& coll );
GISolver CreateAdvancePosition(DynamicalState& pq);
}


#endif
