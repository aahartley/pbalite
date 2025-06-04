//*******************************************************************
//
//   CollisionTriangle.C
//
//   Collision Triangle determines if inside triangle area
//   handles reflection
//
//
//
//*******************************************************************

#include "CollisionTriangle.h"

#include <iostream>
#include <algorithm>

namespace pba
{

CollisionTriangle::CollisionTriangle(const Vector& v0, const Vector& v1, const Vector& v2)
  : v0_(v0),
    v1_(v1),
    v2_(v2)
{
    e1_ = v1_ - v0_;
    e2_ = v2_ - v0_;
    e3_ = v2_ - v1_; //e3_ = e2_-e1_;
    ComputeNormal();

    double xllc = v0_.X();
    double yllc = v0_.Y();
    double zllc = v0_.Z();
    if (v1_.X() < xllc) { xllc = v1_.X(); }
    if (v2_.X() < xllc) { xllc = v2_.X(); }
    if (v1_.Y() < yllc) { yllc = v1_.Y(); }
    if (v2_.Y() < yllc) { yllc = v2_.Y(); }
    if (v1_.Z() < zllc) { zllc = v1_.Z(); }
    if (v2_.Z() < zllc) { zllc = v2_.Z(); }

    double xurc = v0_.X();
    double yurc = v0_.Y();
    double zurc = v0_.Z();
    if (v1_.X() > xurc) { xurc = v1_.X(); }
    if (v2_.X() > xurc) { xurc = v2_.X(); }
    if (v1_.Y() > yurc) { yurc = v1_.Y(); }
    if (v2_.Y() > yurc) { yurc = v2_.Y(); }
    if (v1_.Z() > zurc) { zurc = v1_.Z(); }
    if (v2_.Z() > zurc) { zurc = v2_.Z(); }

    aabb_ = AABB(Vector(xllc,yllc,zllc), Vector(xurc,yurc,zurc));
}

void CollisionTriangle::ComputeNormal()
{
    Vector n = e1_ ^ e2_;
    un_normal_ = n;
    n.normalize();
    normal_ = n;
}
 
bool CollisionTriangle::Hit(
    const Vector& x_0, const Vector& x_u, const Vector& vel, const double dt,
    Vector& xh_cand, double& dt_cand, const float rad) const
{
    float fx0 = normal_ * (x_0 - v0_);
    float fxu = normal_ * (x_u - v0_);
    float nv = normal_ * vel;
    // moving parallel to normal 
    if (std::abs(nv) == 1.0e-6) 
    {
        return false; 
    }

    bool hit = false;
    // end pos in plane, opposite sides, start pos on/behind plan while moving opposite of normal
    if ((fxu == 0) || (fx0 * fxu < 0) || (fx0 <= 0 && nv < 0)) 
    {
        hit = true;
    }

    if (hit) 
    {
        //if starting slightly below or above plane and moving opposite of normal
        if (fx0 >= -1e-6 && fx0 <= 1e-6 && nv < 0) 
        {
            int sign = (dt > 0) ? 1 : -1;
            if (fx0 > 0) dt_cand = 1.0e-6 * sign;
            else dt_cand = 0;
            if (fx0 < 0) xh_cand = x_0 + (-fx0 * normal_); //project to surface?
            else xh_cand = x_0;
        } 
        else 
        {
            dt_cand = (normal_ * (v0_ - x_0)) / nv;
            if (std::abs(dt_cand) < 1e-6) return false; // avoid jitter
            xh_cand = x_0 + vel * dt_cand;
        }

        // Barycentric check with tolerance
        const double tol = 1e-6;
        double u = (un_normal_ * ((xh_cand - v0_) ^ e2_)) / (un_normal_ * un_normal_);
        double v = (un_normal_ * (e1_ ^ (xh_cand - v0_))) / (un_normal_ * un_normal_);
        
        if (!(u >= -tol && v >= -tol && u + v <= 1.0 + tol)) 
        {
            return false; // Not inside triangle
        }

        // Validate collision time
        bool valid = (dt >= 0) ? (dt_cand >= 0 && dt_cand <= dt) : (dt_cand <= 0 && dt_cand >= dt);
        bool valid2 = std::abs(dt_cand) <= std::abs(dt);
        if (!valid) 
        {
            std::cout << "Invalid dt_cand: " << dt_cand << " dt: " << dt << "\n";
            return false;
        }
         if (!valid2) 
        {
            std::cout << "Invalid dt_cand222222: " << dt_cand << " dt: " << dt << "\n";
        }

        return true;
    }
    return false;
}

//1 sticky means unchanged
void CollisionTriangle::Handle(
    const Vector& v_s, const double dt, const Vector& x_h, 
    const double dt_h, Vector& x_r, Vector& v_r, const float cs, const float cr) const
{
    //std::cout << "handling\n";
    v_r = (cs * v_s) - ((cs + cr) * normal_ * (normal_ * v_s));
    double remaining_dt =  (dt - dt_h);
    x_r = x_h + v_r * remaining_dt;
}

}//end of pba namespace