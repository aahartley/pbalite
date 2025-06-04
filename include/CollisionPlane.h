//*******************************************************************
//
//   CollisionPlane.h
//
//   Collision Infinite Plane determines if hit plane
//   handles reflection
//
//
//
//*******************************************************************

#ifndef __PBA_COLLISIONPLANE_H__
#define __PBA_COLLISIONPLANE_H__
  
#include "Vector.h"


namespace pba
{
/*!
  Collision Triangle does barycentric coords check
  handles refelction
 */   
class CollisionInfinitePlane
{
  public:
    CollisionInfinitePlane(){ p0_ = Vector(0,0,0); normal_ = Vector(0,0,0); }
    CollisionInfinitePlane(const Vector& normal, const Vector& p0);
    ~CollisionInfinitePlane(){}
    
    //! If in plane, returns true, then xh_cand and dt_cand are filled with the hit point and hit time.
    bool Hit(
      const Vector& x_0, const Vector& x_u, const Vector& v, double dt, Vector& xh_cand, double& dt_cand, float rad) const;
    
    //! Takes in hit data and returns reflected position and velocity
    void Handle(
      const Vector& v_s, double dt, const Vector& x_h, double dt_h, Vector& x_r, Vector& v_r, float cs, float cr) const;   

    const Vector& getP0(){ return p0_; } 
    const Vector& getNormal(){ return normal_; } 

  private:
    Vector normal_;
    Vector p0_;
  
};
  
}
  
#endif