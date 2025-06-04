//*******************************************************************
//
//   CollisionSurface.C
//
//   Collision Surface stores triangles, and returns data for first collisin
//
//
//
//*******************************************************************

#include "CollisionHandler.h"

namespace pba
{

CollisionSurface::CollisionSurface() 
  : visible_ (true),
    wireframe_(true),
    points_(false),
    coeff_of_restitution_(1.0),
    coeff_of_sticky_(1.0) {}


void CollisionSurface::AddPlane(const CollisionInfinitePlane& plane)
{
    plane_elements_.push_back(plane);
}

void CollisionSurface::AddTriangle(const CollisionTriangle& tri)
{
    tri_elements_.push_back(tri);
}

bool CollisionSurface::Hit(
    const Vector& x_0, const Vector& x_u, const Vector& v_0, const double dt, CollisionData& data, const float rad)
{
    //init data
    bool hit = false;
    data.hit_time = 2.0 * dt;
    data.hit_tri = false;
    data.hit_plane = false;
    for (size_t i = 0; i < tri_elements_.size(); ++i)
    {
        double dt_candidate = data.hit_time;
        Vector x_h_candidate;
        //check for coll for all planes 
        if (tri_elements_[i].Hit(x_0, x_u, v_0, dt, x_h_candidate, dt_candidate, rad))
        {
            hit = true;
            // Find the earliest collision, i.e. the one that happened longest into the past.
            if (std::fabs(dt_candidate) < std::fabs(data.hit_time))
            {
                data.hit_time = dt_candidate;
                data.x_h = x_h_candidate;
                data.hit_tri = true;
                data.tri = &tri_elements_[i];
                data.hit_index = i;
                // std::cout << "hit_tri: " << data.hit_tri << '\n';
                // std::cout << "hit_index: " << data.hit_index << '\n';
            }
        }
    }

    // hit = false;
    // data.hit_time = 2.0 * dt;
    // data.hit_plane = false;
    // for (size_t i = 0; i < plane_elements_.size(); ++i)
    // {
    //     double dt_candidate = data.hit_time;
    //     Vector x_h_candidate;
    //     //check for coll for all planes 
    //     if (plane_elements_[i].Hit(x_0, x_u, v_0, dt, x_h_candidate, dt_candidate, rad))
    //     {
    //         hit = true;
    //         // Find the earliest collision, i.e. the one that happened longest into the past.
    //         if (std::fabs(dt_candidate) < std::fabs(data.hit_time))
    //         {
    //             data.hit_time = dt_candidate;
    //             data.x_h = x_h_candidate;
    //             data.hit_plane = true;
    //             data.plane = &plane_elements_[i];
    //             data.hit_index = i;
    //             // std::cout << "hit_tri: " << data.hit_tri << '\n';
    //             // std::cout << "hit_index: " << data.hit_index << '\n';
    //         }
    //     }
    // }

    return hit;
}

//plane logic, useful?
// bool CollisionSurfaceRaw::hit(const Vector& X0, const Vector& XU, const Vector& V, const double dt, CollisionData& data, float radius) const
// {
//     //init data
//     bool hit = false;
//     data.hit_time = 2.0 * dt;
//     data.hit_tri = false;
//     data.hit_plane = false;
//     for( size_t i = 0; i < plane_elements.size(); i++ )
//     {
//         double dtH_candidate = data.hit_time;
//         Vector XH_candidate;
//         //check for coll for all planes 
//         if(plane_elements[i].hit( X0,XU, V, dt, XH_candidate, dtH_candidate, radius))
//         {
//             hit = true;
//             // Find the earliest collision, i.e. the one that happened longest into the past.
//             if( std::fabs(dtH_candidate) < std::fabs(data.hit_time) )
//             {
//                 data.hit_time = dtH_candidate;
//                 data.XH = XH_candidate;
//                 data.hit_plane = true;
//                 data.plane = plane_elements[i];
//                 data.hit_index = i;
//             }
//         }
//     }

//     return hit;
// }

// TODO: add custom inits?
CollisionSurfacePtr CreateCollisionSurface()
{
   return std::make_unique<CollisionSurface>();
}

}// end of pba namespace