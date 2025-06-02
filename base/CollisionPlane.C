#include "CollisionPlane.h"
#include <iostream>
#include <algorithm>
using namespace pba;

CollisionInfinitePlane::CollisionInfinitePlane(const Vector& norm, const Vector& p0) :
  normal(norm),
  P0(p0)
  {}



   
bool CollisionInfinitePlane::hit( const Vector& X0, const Vector& Xu, const Vector& V, const double dt, Vector& XH_cand, double& dtH_cand, float rad) const
{
    Vector effectiveV = (dt >= 0) ? V : -V;

    bool hit = false;
    float fx0 = normal * ((X0) - P0);
    float fxu = normal * ((Xu) - P0);
    float vn  = normal * effectiveV ;

    if (V*V < 1.0e-12) return hit; // no motion, no collision

    // Check collision for forward/backward motion
    bool approaching = (vn < 0 && dt >= 0) || (vn > 0 && dt < 0);
    if (fx0 <= 0.0f && approaching) 
    {
        dtH_cand = 0;
        XH_cand  = X0;
        hit = true;
        return hit;
    }
    else if((fxu == 0 || fx0 * fxu < 0 ))/*|| (fxu <=0 && fx0 <=0)*/ hit = true;
    if(hit)
    {
        float nv = normal * effectiveV;
        if (std::abs(nv) < 1e-6) nv = (dt >= 0) ? 1e-6 : -1e-6;
        XH_cand = X0 + effectiveV * ( (normal * (P0 - X0)) / nv );
        dtH_cand = (normal * (P0-X0)) / nv; // (positive dt results: should be <= dt )
      
        //std::cout << "Dth: " << dtH_cand << '\n';
        //if(dtH_cand > dt) std::cout << "dth bigger than dt!!" << dtH_cand << ' ' << dt << '\n';

        bool valid = (dt >= 0) ? (dtH_cand >= 0 && dtH_cand <= dt) : (dtH_cand <= 0 && dtH_cand >= dt);

        if (!valid) 
        {
            std::cout << "Invalid dtH_cand: " << dtH_cand << " dt: " << dt << "\n";
        }
    }

    return hit;
}

//1 sticky means unchanged
void CollisionInfinitePlane::handle(const Vector& XS, const Vector& VS, const double& dt, 
    const Vector& XH, const double& dtH, Vector& XR, Vector& VR, float cs, float cr) const
{
    VR = (cs * VS) - ((cs + cr) * normal * (normal * VS));
    double remaining_dt = (dt >= 0) ? (dt - dtH) : (dtH - dt);
    XR = XH + VR * remaining_dt;
}