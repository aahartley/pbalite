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
    const Vector& x_0, const Vector& x_u, const Vector& v, const double dt,
    Vector& xh_cand, double& dt_cand, const float rad) const
{
    bool hit = false;
    float fx0 = normal_ * ((x_0) - v0_);
    float fxu = normal_ * ((x_u) - v0_);
    float vn  = normal_ * v ;

    if (v*v < 1.0e-12) return hit; // no motion, no collision

    // Check collision for forward/backward motion
    bool approaching = (dt >= 0) ? (vn < 0) : (vn > 0);
    //init penetration , prevent tunneling
    if (fx0 <= 0.0f && approaching) 
    {
        dt_cand = 0;
        xh_cand  = x_0;
        hit = true;
        return hit;
    }
    else if ((fxu == 0 || fx0 * fxu < 0 )) hit = true;
    if (hit)// plane of tri hit, eval barycentric coords
    {
        float nv = normal_ * v;  
        if (std::abs(nv) < 1e-6) nv = (v * normal_ < 0) ? -1e-6 : 1e-6;
        
        xh_cand = x_0 + v * ( (normal_ * (v0_ - x_0)) / nv );
        dt_cand = (normal_ * (v0_-x_0)) / nv; // (positive dt results: should be <= dt )
      
        //eval barycentric coords
        double u = (un_normal_ * ((xh_cand - v0_) ^ e2_)) / (un_normal_ * un_normal_);
        double v = (un_normal_ * (e1_ ^ (xh_cand - v0_))) / (un_normal_ * un_normal_);
        if (!((u >=0 && u <=1) && (v >=0 && v <=1) && (u+v >=0 && u+v <=1))) hit = false; // not in tri

        bool valid = (dt >= 0) ? (dt_cand >= 0 && dt_cand <= dt) : (dt_cand <= 0 && dt_cand >= dt);

        if (!valid) 
        {
            std::cout << "Invalid dt_cand: " << dt_cand << " dt: " << dt << "\n";
        }
    }

    return hit;
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