
#ifndef ____PBA_AABB_H____
#define ____PBA_AABB_H____

#include "Vector.h"
#include "Ray.h"

namespace pba
{

class AABB
{
  public:
    AABB();
    AABB(const Vector& LLC, const Vector& URC);
    ~AABB();

    //! Test whether the input AABB intersects this AABB
    const bool intersects( const AABB& aabb ) const;

    // //! Return an AABB for the region of intersection of the input AABB and this AABB.
    // AABB Intersection( const AABB& aabb ) const;
    // //! Return the union AABB for the input AABB and this AABB
    // AABB Union( const AABB& aabb ) const;

    //! Divide this AABB into two equal AABBs along the component direction
    void split( const int component, AABB& aabb1, AABB& aabb2 ) const;

    // //! Return the volume of this AABB
    // const double volume() const;

    //! Return the distance at which the input ray intersects this AABB
    bool ray_intersect(const Ray &r, float t0, float t1) const ;

    const Vector& getURC() const { return urc;}
    const Vector& getLLC() const { return llc;}

    const Vector& getBounds(int sign) const;

    const bool isInside( const Vector& P ) const;



  private:
    Vector urc;
    Vector llc;



};





}



#endif