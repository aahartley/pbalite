#ifndef __PBA_COLLISIONHANDLER_H__
#define __PBA_COLLISIONHANDLER_H__
  
#include "CollisionSurface.h"
#include "DynamicalState.h"
#include "BVH.h"

#include <iostream>

namespace pba
{
  
class CollisionHandler
{
  public:
  
    CollisionHandler() : surf(0), bvh(0), useBVH(true){}
    virtual ~CollisionHandler(){}

    virtual void handle_collisions(const double dt, DynamicalState& s){ std::cout << "CollisionHandler::handle_collisions(double,DynamicalState) called\n"; }
    void set_collision_surface(CollisionSurface& c);
    void set_bvh(BVH& b);

    void use_bvh() { useBVH = true; }
    void dont_use_bvh() { useBVH = false; }
  
  protected:
    CollisionSurface surf;
    BVH bvh;
    bool useBVH;
  
 };


class ElasticCollisionHandler : public CollisionHandler
{
  public:

    ElasticCollisionHandler();
    ~ElasticCollisionHandler();

    void handle_collisions(const double dt, DynamicalState& PQ);

};


}

#endif
