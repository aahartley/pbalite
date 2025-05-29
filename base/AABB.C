#include "AABB.h"

using namespace pba;

AABB::AABB(){}

AABB::AABB(const Vector& LLC, const Vector& URC): urc(URC), llc(LLC){}

AABB::~AABB(){}

const bool AABB::isInside( const Vector& P ) const
{
    for( int i =0; i < 3; i++ )
    {
       if( P[i] < llc[i] ){ return false; }
       if( P[i] > urc[i] ){ return false; }
    }
    return true;
}

void AABB::split( const int component, AABB& aabb0, AABB& aabb1 ) const
{
    Vector center = (llc + urc) * 0.5;
    //lower
    Vector urc0 = urc;
    urc0[component] = center[component];
    aabb0.llc = llc;
    aabb0.urc = urc0;
    //upper
    Vector llc1 = llc;
    llc1[component] = center[component];
    aabb1.llc = llc1;
    aabb1.urc = urc;
}

const bool AABB::intersects( const AABB& aabb ) const
{
    if( llc[0] > aabb.urc[0] ){ return false; }
    if( llc[1] > aabb.urc[1] ){ return false; }
    if( llc[2] > aabb.urc[2] ){ return false; }
    if( aabb.llc[0] > urc[0] ){ return false; }
    if( aabb.llc[1] > urc[1] ){ return false; }
    if( aabb.llc[2] > urc[2] ){ return false; }
    return true;
}

const Vector& AABB::getBounds(int sign) const
{
    if(sign) return urc;
    else return llc;
}

//https://people.csail.mit.edu/amy/papers/box-jgt.pdf
bool AABB::ray_intersect(const Ray& r, float t0, float t1) const 
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

