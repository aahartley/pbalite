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
