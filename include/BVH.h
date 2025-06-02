//*******************************************************************
//
//   BVH.h
//
//  Bounding Volume Hierarchy for continous collison detection
//
//
//
//*******************************************************************

#ifndef ____PBA_BVH_H____
#define ____PBA_BVH_H____

#include "Vector.h"
#include "CollisionSurface.h"
#include "AABB.h"

#include <memory>
#include <vector>

namespace pba
{

/*!
  BVH is an Bounding Volume Hierarchy, tree data structure.
  Reset BVH if destryoing collision surface
 */
class BVH
{
  public:
    //! Define AABB, current level (init to 0), max level for divison, min number of objects per leaf
    BVH(const Vector& llc, const Vector& urc, int level, int max_levels, int min_objects);
    ~BVH();

    //! Divide BVH until max level depth, only splits nodes with triangles
    void Divide();

    //! Add triangles from Collision Surface
    void AddObject(CollisionSurface& s);
    void AddObject(CollisionTriangle& s);
    void AddObject(CollisionTriangle* s);


    //! Test if position line segment intersects bvh, return Collision Data with triangle of first hit
    bool Hit(const Ray& ray, const Vector& X0, const Vector& Xu, const Vector& V, double max_dt, CollisionData& c_data) const;
    //bool hit( const RigidBodyState& s, const size_t i, const double tmax, CollisionData& t ) const;

    //! Return number of objects in leaf
    size_t nbObjects() const { return object_list_.size(); }

  private:
    AABB aabb_;
    BVH* node1_;
    BVH* node2_;
    int level_;
    int max_levels_;
    int min_objects_;
    std::vector<CollisionTriangle*> object_list_;
};

using BVHPtr = std::unique_ptr<BVH>;

inline BVHPtr CreateBVH(const Vector& llc, const Vector& urc, int lvl, int maxlvl, int minobj) 
{
    return std::make_unique<BVH>(llc, urc, lvl, maxlvl, minobj);
}

}//end of pba namespace

#endif