//*******************************************************************
//
//   AABB.C
//
//  Axis Aligned Bounding Box for intersection testing
//
//
//
//*******************************************************************

#include "AABB.h"

namespace pba
{

AABB::AABB(){}

AABB::AABB(const Vector& llc, const Vector& urc)
  : llc_(llc),
    urc_(urc) {}

AABB::~AABB(){}

bool AABB::IsInside(const Vector& P) const
{
    for(int i = 0; i < 3; ++i)
    {
       if (P[i] < llc_[i]) { return false; }
       if (P[i] > urc_[i]) { return false; }
    }
    return true;
}

void AABB::Split(int component, AABB& aabb0, AABB& aabb1) const
{
    Vector center = (llc_ + urc_) * 0.5;
    //lower
    Vector urc0 = urc_;
    urc0[component] = center[component];
    aabb0.llc_ = llc_;
    aabb0.urc_ = urc0;
    //upper
    Vector llc1 = llc_;
    llc1[component] = center[component];
    aabb1.llc_ = llc1;
    aabb1.urc_ = urc_;
}

bool AABB::Intersects(const AABB& aabb) const
{
    if (llc_[0] > aabb.urc_[0]) { return false; }
    if (llc_[1] > aabb.urc_[1]) { return false; }
    if (llc_[2] > aabb.urc_[2]) { return false; }
    if (aabb.llc_[0] > urc_[0]) { return false; }
    if (aabb.llc_[1] > urc_[1]) { return false; }
    if (aabb.llc_[2] > urc_[2]) { return false; }
    return true;
}

const Vector& AABB::getBounds(int sign) const
{
    if (sign) return urc_;
    else return llc_;
}

//https://people.csail.mit.edu/amy/papers/box-jgt.pdf
bool AABB::RayIntersects(const Ray& r, float t0, float t1) const 
{
    float tmin = (getBounds(r.sign[0]).X() - r.origin.X()) * r.inv_direction.X();
    float tmax = (getBounds(1 - r.sign[0]).X() - r.origin.X()) * r.inv_direction.X();

    float tymin = (getBounds(r.sign[1]).Y() - r.origin.Y()) * r.inv_direction.Y();
    float tymax = (getBounds(1 - r.sign[1]).Y() - r.origin.Y()) * r.inv_direction.Y();

    if ((tmin > tymax) || (tymin > tmax))
        return false;

    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;

    float tzmin = (getBounds(r.sign[2]).Z() - r.origin.Z()) * r.inv_direction.Z();
    float tzmax = (getBounds(1 - r.sign[2]).Z() - r.origin.Z()) * r.inv_direction.Z();

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;

    return (tmin < t1) && (tmax > t0);
}

} //end of pba namespace