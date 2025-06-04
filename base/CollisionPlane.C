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

#include "CollisionPlane.h"

#include <iostream>
#include <algorithm>

namespace pba
{

CollisionInfinitePlane::CollisionInfinitePlane(const Vector& norm, const Vector& p0) 
  : normal_(norm),
    p0_(p0) {}



   
bool CollisionInfinitePlane::Hit(
    const Vector& x_0, const Vector& x_u, const Vector& vel, double dt, Vector& xh_cand, double& dt_cand, float rad) const
{
    bool hit = false;
    float fx0 = normal_ * ((x_0) - p0_);
    float fxu = normal_ * ((x_u) - p0_);
    std::cout << "fx0: " << fx0 << " fxu: " << fxu << '\n';
    const float epsilon = 1.0e-2;
    if((fxu <= epsilon && fxu >= -epsilon) || fx0 * fxu < 0 /*|| (fxu <=0 && fx0 <=0)*/) hit = true;
    if(hit)
    {
        float nv = normal_ * vel;
        if (std::abs(nv) < 1e-6) 
        {
            std::cout << nv << "\n";
            nv = (dt >= 0) ? 1e-6 : -1e-6;
            return false;
        }
        xh_cand = x_0 + vel * ( (normal_ * (p0_ - x_0)) / nv );
        dt_cand = (normal_ * (p0_-x_0)) / nv; // (positive dt results: should be <= dt )
      
        //std::cout << "Dth: " << dtH_cand << '\n';
        //if(dtH_cand > dt) std::cout << "dth bigger than dt!!" << dtH_cand << ' ' << dt << '\n';

        bool valid = (dt >= 0) ? (dt_cand >= 0 && dt_cand <= dt) : (dt_cand <= 0 && dt_cand >= dt);

        if (!valid) 
        {
            std::cout << "Invalid dtH_cand: " << dt_cand << " dt: " << dt << "\n";
        }
    }

    return hit;
}

//1 sticky means unchanged
void CollisionInfinitePlane::Handle(
    const Vector& v_s, const double dt, const Vector& x_h, 
    const double dt_h, Vector& x_r, Vector& v_r, const float cs, const float cr) const
{
    //std::cout << "handling\n";
    v_r = (cs * v_s) - ((cs + cr) * normal_ * (normal_ * v_s));
    double remaining_dt =  (dt - dt_h);
    x_r = x_h + v_r * remaining_dt;
}

}