#ifndef __PBA_FORCE_H__
#define __PBA_FORCE_H__

#include "DynamicalState.h"
#include <memory>

namespace pba
{
class ForceBase
{
  public:
    ForceBase(){}
    //! compute forces on the dynamical state and update the accel attribute(s)
    virtual void compute(DynamicalState& s, const double dt) = 0;
    virtual void compute(SPHState& s, const double dt) = 0;

    virtual ~ForceBase(){};
};
typedef std::shared_ptr<pba::ForceBase> Force;

}



#endif