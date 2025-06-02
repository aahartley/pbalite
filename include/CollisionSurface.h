//*******************************************************************
//
//   CollisionSurface.h
//
//   Collision Surface stores triangles, and returns data for first collisin
//
//
//
//*******************************************************************

#ifndef __PBA_COLLISIONSURFACE_H__
#define __PBA_COLLISIONSURFACE_H__

#include "CollisionTriangle.h"
#include "CollisionPlane.h"

#include <vector>
#include <memory>
  
namespace pba
{
  
struct CollisionData
{
  double hit_time;
  CollisionTriangle* tri;
  CollisionInfinitePlane plane;
  bool hit_tri;
  bool hit_plane;
  size_t hit_index;
  Vector x_h; 
};

/*!
  Collision Surface returns collision data for first collision
  Controls graphics of surface
  Surface owns list of triangles
  TODO: if tris decay toggle bool, delete and rebuild bvh later
 */  
class CollisionSurface
{
  public:
    CollisionSurface();
    ~CollisionSurface(){}
  
    void AddTriangle(const CollisionTriangle& t);
    void AddPlane(const CollisionInfinitePlane& p);
    void Clear() { plane_elements_.clear(); tri_elements_.clear(); } 

    size_t triangle_size() const { return tri_elements_.size(); }
    size_t plane_size() const { return plane_elements_.size(); }
  
    bool Hit(const Vector& X0, const Vector& XU, const Vector& V, double tmax, CollisionData& t, float radius);
  
    CollisionTriangle& get_triangle(size_t i) { return tri_elements_[i]; }
    const CollisionTriangle& get_triangle(size_t i) const { return tri_elements_[i]; }

    CollisionInfinitePlane& get_plane(size_t i) { return plane_elements_[i]; }
    const CollisionInfinitePlane& get_plane(size_t i) const { return plane_elements_[i]; }

    void toggle_visible() { visible_ = !visible_; }
    bool is_visible() const { return visible_; }
    void toggle_wireframe() { wireframe_ = !wireframe_; }
    bool use_wireframe() const { return wireframe_; }
    void toggle_points() { points_ = !points_; }
    bool use_points() const { return points_; }
  
    void set_coeff_restitution(double v) { coeff_of_restitution_ = v; }
    const double& coeff_restitution() const { return coeff_of_restitution_; }
  
    void set_coeff_sticky(double v) { coeff_of_sticky_ = v; }
    const double& coeff_sticky() const { return coeff_of_sticky_; }
    
  private:
    bool visible_;
    bool wireframe_;
    bool points_;
  
    double coeff_of_restitution_;
    double coeff_of_sticky_;
  
    std::vector<CollisionTriangle> tri_elements_;
    std::vector<CollisionInfinitePlane> plane_elements_;
 
};
  
using CollisionSurfacePtr = std::unique_ptr<CollisionSurface>;

CollisionSurfacePtr CreateCollisionSurface();
  
  
} //end of pba namespace
  
#endif