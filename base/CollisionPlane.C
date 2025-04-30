#include "CollisionPlane.h"

using namespace pba;

CollisionInfinitePlane::CollisionInfinitePlane(const Vector& norm, const Vector& p0) :
  normal(norm),
  P0(p0)
  {}



   
bool CollisionInfinitePlane::hit( const Vector& X0, const Vector& Xu, const Vector& V, const double dt, Vector& XH_cand, double& dtH_cand, float rad) const
{
    bool hit = false;
    float fx0 = normal * ((X0) - P0);
    float fxu = normal * ((Xu) - P0);
    float vn  = normal * V ;

    // if(fxu >=0)
    // {
    //     fxu = (Xu - P0) * normal - rad;
    // }
    // else fxu = (Xu - P0) * normal + rad;

    // if(fx0 >=0)
    // {
    //     fx0 = (X0 - P0) * normal - rad;
    // }
    // else fx0 = (X0 - P0) * normal + rad;
    if (V*V < 1.0e-12)
        return false; // no motion, no collision
    if (fx0 <= 0.0f && vn < 0.0f && V*V > 1.0e-12) 
    {
        dtH_cand = 0;
        XH_cand  = X0;
        hit = true;
        return hit;
    }
    else if((fxu == 0 || fx0 * fxu < 0 ))/*|| (fxu <=0 && fx0 <=0)*/ hit = true;
    if(hit)
    {
        //std::cout << V.X() << ' ' << V.Y() << ' ' << V.Z() << '\n';
        // if(V.isZero() == true)
        // {
        //     std::cout << "Zero? " << V.isZero() << '\n';
        //     //float epsilon = 0.06041;
        //     XH_cand = X0 + (normal * rad); //- (normal*epsilon);
        //     dtH_cand = dt;
        //     return hit;
        // }       
        // if(V.isZero())
        // {
        //    std::cout << V.X() << ' ' << V.Y() << ' ' << V.Z() << '\n';
        //    std::cout << "Zero v: "<< normal*V << '\n';
        // } 
        float nv = normal *V;
        if(nv==0)nv = normal * Vector(1.0e-6, 1.0e-6, 1.0e-6);
        XH_cand = X0 + V * ( (normal * (P0 - X0)) / (nv) );
        dtH_cand = (normal * (P0-X0)) / (nv); // should be <= dt
      
        //std::cout << "Dth: " << dtH_cand << '\n';
        if(dtH_cand > dt) std::cout << "dth bigger than dt!!" << dtH_cand << ' ' << dt << '\n';
    }

    return hit;
}

//1 sticky means unchanged
void CollisionInfinitePlane::handle(const Vector& XS, const Vector& VS, const double& dt, 
    const Vector& XH, const double& dtH, Vector& XR, Vector& VR, float cs, float cr) const
{
    VR = (cs * VS) - ((cs + cr) * normal * (normal * VS));
    XR = XH + VR * (dt - dtH);
}