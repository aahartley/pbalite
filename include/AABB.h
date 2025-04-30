
#ifndef ____PBA_AABB_H____
#define ____PBA_AABB_H____

#include "Vector.h"

namespace pba
{

class AABB
{
  public:
    AABB();
    AABB(Vector& URC, Vector& LLC);
    ~AABB();

    const Vector& getURC() const { return urc;}
    const Vector& getLLC() const { return llc;}

    const bool isInside( const Vector& P ) const;



  private:
    Vector urc;
    Vector llc;



};





}



#endif