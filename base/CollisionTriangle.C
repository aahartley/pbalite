#include "CollisionTriangle.h"

using namespace pba;

CollisionTriangleRaw::CollisionTriangleRaw(const Vector& V0, const Vector& V1, const Vector& V2) :
  v0(V0),
  v1(V1),
  v2(V2)
{
    e1 = v1 - v0;
    e2 = v2 - v0;
    e3 = v2 - v1; //e3 = e2-e1;
    compute_normal();

    double xllc = v0.X();
    double yllc = v0.Y();
    double zllc = v0.Z();
    if( v1.X() < xllc ){ xllc = v1.X(); }
    if( v2.X() < xllc ){ xllc = v2.X(); }
    if( v1.Y() < yllc ){ yllc = v1.Y(); }
    if( v2.Y() < yllc ){ yllc = v2.Y(); }
    if( v1.Z() < zllc ){ zllc = v1.Z(); }
    if( v2.Z() < zllc ){ zllc = v2.Z(); }

    double xurc = v0.X();
    double yurc = v0.Y();
    double zurc = v0.Z();
    if( v1.X() > xurc ){ xurc = v1.X(); }
    if( v2.X() > xurc ){ xurc = v2.X(); }
    if( v1.Y() > yurc ){ yurc = v1.Y(); }
    if( v2.Y() > yurc ){ yurc = v2.Y(); }
    if( v1.Z() > zurc ){ zurc = v1.Z(); }
    if( v2.Z() > zurc ){ zurc = v2.Z(); }

    aabb = AABB(Vector(xllc,yllc,zllc),Vector(xurc,yurc,zurc));
}

void CollisionTriangleRaw::compute_normal()
{
    Vector n = e1 ^ e2;
    un_normal = n;
    n.normalize();
    normal = n;
}


   
bool CollisionTriangleRaw::hit( const Vector& X0, const Vector& Xu, const Vector& V, const double dt, Vector& XH_cand, double& dtH_cand, float rad) const
{
    Vector effectiveV = V;

    bool hit = false;
    float fx0 = normal * ((X0) - v0);
    float fxu = normal * ((Xu) - v0);
    float vn  = normal * effectiveV ;

    if (V*V < 1.0e-12) return hit; // no motion, no collision

    // Check collision for forward/backward motion
    bool approaching = (dt >= 0) ? (vn < 0) : (vn > 0);
    //init penetration , prevent tunneling
    if (fx0 <= 0.0f && approaching) 
    {
        dtH_cand = 0;
        XH_cand  = X0;
        hit = true;
        return hit;
    }
    else if((fxu == 0 || fx0 * fxu < 0 ))/*|| (fxu <=0 && fx0 <=0)*/ hit = true;
    if(hit)// plane of tri hit, eval barycentric coords
    {
        float nv = normal * V;  // Use real velocity (no flipping)
        if (std::abs(nv) < 1e-6) nv = (V * normal < 0) ? -1e-6 : 1e-6;
        
        XH_cand = X0 + effectiveV * ( (normal * (v0 - X0)) / nv );
        dtH_cand = (normal * (v0-X0)) / nv; // (positive dt results: should be <= dt )
      
        //eval barycentric coords
        double u = (un_normal * ((XH_cand - v0) ^ e2)) / (un_normal * un_normal);
        double v = (un_normal * (e1 ^ (XH_cand - v0))) / (un_normal * un_normal);
        if( !((u >=0 && u <=1) && (v >=0 && v <=1) && (u+v >=0 && u+v <=1)) ) hit = false; // not in tri
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
void CollisionTriangleRaw::handle(const Vector& XS, const Vector& VS, const double& dt, 
    const Vector& XH, const double& dtH, Vector& XR, Vector& VR, float cs, float cr) const
{
    VR = (cs * VS) - ((cs + cr) * normal * (normal * VS));
    double remaining_dt =  (dt - dtH);
    XR = XH + VR * remaining_dt;
}
CollisionTriangle pba::makeCollisionTriangle(  const Vector& p0, const Vector& p1, const Vector& p2 )
{
   return CollisionTriangle( new CollisionTriangleRaw(p0,p1,p2) );
}