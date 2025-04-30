
#include "CollisionHandler.h"

using namespace pba;

CollisionSurfaceRaw::CollisionSurfaceRaw() :
   visible (true),
   wireframe(true),
   points(false),
   coeff_of_restitution (1.0),
   coeff_of_sticky (1.0)
{}


void CollisionSurfaceRaw::addPlane(const CollisionInfinitePlane& plane)
{
    plane_elements.push_back(plane);
}


bool CollisionSurfaceRaw::hit(const Vector& X0, const Vector& XU, const Vector& V, const double dt, CollisionData& data, float radius) const
{
    //init data
    bool hit = false;
    data.hit_time = 2.0 * dt;
    data.hit_tri = false;
    data.hit_plane = false;
    for( size_t i=0;i<plane_elements.size();i++ )
    {
        double dtH_candidate = data.hit_time;
        Vector XH_candidate;
        //check for coll for all planes 
        if(plane_elements[i].hit( X0,XU, V, dt, XH_candidate, dtH_candidate, radius))
        {
            hit = true;
            // Find the earliest collision, i.e. the one that happened longest into the past.
            if( std::fabs(dtH_candidate) < std::fabs(data.hit_time) )
            {
                data.hit_time = dtH_candidate;
                data.XH = XH_candidate;
                data.hit_plane = true;
                data.plane = plane_elements[i];
                data.hit_index = i;
            }
        }
    }

    return hit;
}

CollisionSurface pba::makeCollisionSurface()
{
   return std::make_shared<CollisionSurfaceRaw>();
}
