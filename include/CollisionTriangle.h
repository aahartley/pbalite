//*******************************************************************
//
//   CollisionTriangle.h
//
//   Collision Triangle determines if inside triangle area
//   handles reflection
//
//
//
//*******************************************************************

#ifndef __PBA_COLLISION_TRI_H__
#define __PBA_COLLISION_TRI_H__
  
#include "Vector.h"
#include "AABB.h"
  
namespace pba
{
/*!
  Collision Triangle does barycentric coords check
  handles refelction
 */   
class CollisionTriangle
{
  public:
    CollisionTriangle() 
      : v0_(Vector(0,0,0)), 
        v1_(Vector(0,0,0)), 
        v2_(Vector(0,0,0))
    {
        normal_ = Vector(0,0,0),
        un_normal_ = normal_;
        e1_ = v1_ - v0_;
        e2_ = v2_ - v0_;
        e3_ = v2_ - v1_; //e3 = e2-e1;
        aabb_ = AABB(v0_, v1_);
    }
    CollisionTriangle(const Vector& v0, const Vector& v1, const Vector& v2);
    ~CollisionTriangle(){}
    
    void ComputeNormal();

    //! If in triangle, returns true, then xh_cand and dt_cand are filled with the hit point and hit time.
    bool Hit(
      const Vector& x_0, const Vector& x_u, const Vector& v, double dt, Vector& xh_cand, double& dt_cand, float rad) const;
    
    //! Takes in hit data and returns reflected position and velocity
    void Handle(
      const Vector& v_s, double dt, const Vector& x_h, double dt_h, Vector& x_r, Vector& v_r, float cs, float cr) const;  

    const Vector& v0() const { return v0_; } 
    const Vector& v1() const { return v1_; } 
    const Vector& v2() const { return v2_; }

    const Vector& e1() const { return e1_; } 
    const Vector& e2() const { return e2_; } 
    const Vector& e3() const { return e3_; } 

    const Vector& normal() const { return normal_; } 
    const AABB& aabb() { return aabb_; }

  private:
    Vector normal_;
    Vector un_normal_;
    Vector v0_;
    Vector v1_;
    Vector v2_;

    Vector e1_;
    Vector e2_;
    Vector e3_;
    AABB aabb_;

  
};


}
  
#endif