//*******************************************************************
//
//   AABB.h
//
//  Axis Aligned Bounding Box for intersection testing
//
//
//
//*******************************************************************

#ifndef ____PBA_AABB_H____
#define ____PBA_AABB_H____

#include "Vector.h"
#include "Ray.h"

namespace pba
{

/*!
  AABB is an axis-aligned bounding box, suitable to designating rectangular regions of space.
 */
class AABB
{
  public:
    //! Note that Vector class inits to (0,0,0)
    AABB();
    AABB(const Vector& LLC, const Vector& URC);
    ~AABB() = default;

    //! Test whether the input AABB intersects this AABB
    bool Intersects(const AABB& aabb) const;

    // //! Return an AABB for the region of intersection of the input AABB and this AABB.
    // AABB Intersection( const AABB& aabb ) const;

    // //! Return the union AABB for the input AABB and this AABB
    // AABB Union( const AABB& aabb ) const;

    //! Divide this AABB into two equal AABBs along the component direction
    void Split(int component, AABB& aabb1, AABB& aabb2) const;

    // //! Return the volume of this AABB
    // const double volume() const;

    //! Test if the ray intersects the AABB in the ranges of [t0, t1]
    bool RayIntersects(const Ray &r, float t0, float t1) const ;

    //! Get lower left corner
    const Vector& llc() const { return llc_;}
    //! Get upper right corner
    const Vector& urc() const { return urc_;}

    //! Helper function for ray intersect to determine urc or llc
    const Vector& getBounds(int sign) const;

    //! Test if point is in the AABB
    bool IsInside(const Vector& P) const;

  private:
    Vector llc_;
    Vector urc_;

};

}//end of pba namespace

#endif // ____PBA_AABB_H____