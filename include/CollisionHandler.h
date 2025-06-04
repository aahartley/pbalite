//*******************************************************************
//
//   CollisionHandler.h
//
//   Continous collison detection, resolve all hits during one time step
//
//
//
//*******************************************************************

#ifndef __PBA_COLLISIONHANDLER_H__
#define __PBA_COLLISIONHANDLER_H__
  
#include "CollisionSurface.h"
#include "DynamicalState.h"
#include "BVH.h"

#include <iostream>

namespace pba
{
/*!
  Collision Handler loops until every collision is handled in a time step
 */
class CollisionHandlerBase
{
  public:
    CollisionHandlerBase() 
      : surf_(nullptr),
        bvh_(nullptr), 
        use_bvh_(true) {}
    virtual ~CollisionHandlerBase(){}

    virtual void HandleCollisions(double dt, DynamicalStateData& pq) 
    { 
      std::cout << "CollisionHandler::handle_collisions(double,DynamicalState) called\n"; 
    }

    void set_collision_surface(CollisionSurface* c);
    void set_bvh(const BVH* b);
    void set_use_bvh(bool b) { use_bvh_ = b; }
  
  protected:
    CollisionSurface* surf_;
    const BVH* bvh_;
    bool use_bvh_;
 };

/*!
  Collision Handler loops until every collision is handled in a time step
  Uses sticky and restituion 
 */
class ElasticCollisionHandler : public CollisionHandlerBase
{
  public:
    ElasticCollisionHandler();
    ~ElasticCollisionHandler();

    void HandleCollisions(double dt, DynamicalStateData& pq);
};

}//end of pba namespace

#endif // __PBA_COLLISIONHANDLER_H__
