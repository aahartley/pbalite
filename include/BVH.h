
#ifndef ____PBA_BVH_H____
#define ____PBA_BVH_H____

#include "Vector.h"
#include "CollisionSurface.h"
#include "AABB.h"
#include <memory>
#include <vector>
#include <iostream>

namespace pba
{

class BVHRaw
{
  public:
    BVHRaw( const Vector& llc, const Vector& urc, const int lvl, const int maxlvl, const int minobj );
    ~BVHRaw();
    void Divide();
    void addObject(const CollisionSurface& s );
    void addObject(const CollisionTriangle& s );
    bool hit(const Ray& ray, const Vector& X0, const Vector& Xu, const Vector& V, const double max_dt, CollisionData& c_data ) const;
    //bool hit( const RigidBodyState& s, const size_t i, const double tmax, CollisionData& t ) const;
    const size_t nbObjects() const { return object_list.size(); }
  private:
    AABB aabb;
    BVHRaw* node1;
    BVHRaw* node2;
    int level;
    int max_levels;
    int min_objects;
    std::vector<CollisionTriangle> object_list;
};
typedef std::shared_ptr<BVHRaw> BVH;
BVH makeBVH(const Vector& llc, const Vector& urc, const int lvl, const int maxlvl, const int minobj);

}
#endif