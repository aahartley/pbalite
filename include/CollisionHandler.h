#ifndef __PBA_COLLISIONHANDLER_H__
#define __PBA_COLLISIONHANDLER_H__
  
#include "CollisionSurface.h"
#include "DynamicalState.h"

#include <iostream>

namespace pba
{
  
class CollisionHandler
{
  public:
  
    CollisionHandler() : surf(0){}
    virtual ~CollisionHandler(){}

    virtual void handle_collisions(const double dt, DynamicalState& s){ std::cout << "CollisionHandler::handle_collisions(double,DynamicalState) called\n"; }
    void set_collision_surface(CollisionSurface& c);
  
  protected:
    CollisionSurface surf;
  
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
