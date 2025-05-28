#ifndef __PBA_COLLISION_TRI_H__
#define __PBA_COLLISION_TRI_H__
  
#include "Vector.h"
#include "AABB.h"
#include <iostream>
#include <algorithm>
#include <memory>
  
namespace pba
{
  
class CollisionTriangleRaw
{
  public:
    CollisionTriangleRaw() :  v0(Vector(0,0,0)), v1(Vector(0,0,0)), v2(Vector(0,0,0))
    {
        normal = Vector(0,0,0),
        un_normal = normal;
        e1 = v1 - v0;
        e2 = v2 - v0;
        e3 = v2 - v1; //e3 = e2-e1;
        aabb = AABB(v0, v1);
    }
    CollisionTriangleRaw(const Vector& v0, const Vector& v1, const Vector& v2);
    ~CollisionTriangleRaw(){}
    
    bool hit(const Vector& X0, const Vector& Xu, const Vector& V, const double dt, Vector& XH_cand, double& dtH_cand, float radius) const;
    // bool hit( const Vector& XS, const Vector VS, const double& dt, Vector& XH, double& dtH ) const;
    
    //It returns true, then XH and dtH are filled with the hit point and hit time.
    // Takes in hit data and returns reflected position and velocity
    void handle(const Vector& XS, const Vector& VS, const double& dt, const Vector& XH, const double& dtH, Vector& XR, Vector& VR, float cs, float cr) const;  

    const Vector& getV0(){ return v0; } 
    const Vector& getV1(){ return v1; } 
    const Vector& getV2(){ return v2; }

    const Vector& getE1(){ return e1; } 
    const Vector& getE2(){ return e2; } 
    const Vector& getE3(){ return e3; } 

    const Vector& getNormal(){ return normal; } 
    const AABB& getAABB(){ return aabb; }
    void compute_normal();

  private:
    Vector normal;
    Vector un_normal;
    Vector v0;
    Vector v1;
    Vector v2;

    Vector e1;
    Vector e2;
    Vector e3;
    AABB aabb;

  
};
typedef std::shared_ptr<CollisionTriangleRaw> CollisionTriangle;
CollisionTriangle makeCollisionTriangle(  const Vector& p0, const Vector& p1, const Vector& p2 );

}
  
#endif