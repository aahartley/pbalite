#ifndef __PBA_COLLISIONPLANE_H__
#define __PBA_COLLISIONPLANE_H__
  
#include "Vector.h"

#include <iostream>
#include <algorithm>
  
namespace pba
{
  
class CollisionInfinitePlane
{
  public:
    CollisionInfinitePlane(){ P0=Vector(0,0,0); normal = Vector(0,0,0); }
    CollisionInfinitePlane(const Vector& normal, const Vector& p0);
    ~CollisionInfinitePlane(){}
    
    bool hit(const Vector& X0, const Vector& Xu, const Vector& V, const double dt, Vector& XH_cand, double& dtH_cand, float radius) const;
    // bool hit( const Vector& XS, const Vector VS, const double& dt, Vector& XH, double& dtH ) const;
    
    //It returns true, then XH and dtH are filled with the hit point and hit time.
    // Takes in hit data and returns reflected position and velocity
    void handle(const Vector& XS, const Vector& VS, const double& dt, const Vector& XH, const double& dtH, Vector& XR, Vector& VR, float cs, float cr) const;  

    const Vector& getP0(){ return P0; } 
    const Vector& getNormal(){ return normal; } 

  private:
    Vector normal;
    Vector P0;
  
};
  
}
  
#endif