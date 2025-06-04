//*******************************************************************
//
//   Force.h
//
//  Base class for computing forces on different states
//
//
//
//*******************************************************************

#ifndef __PBA_FORCE_H__
#define __PBA_FORCE_H__

#include "DynamicalState.h"
#include "SPHState.h"
#include "SoftBodyState.h"
#include "RigidBodyState.h"

#include <memory>

namespace pba
{

/*!
  Base class for computing forces on different states
 */
template<typename StateType>
class ForceBase
{
  public:
    ForceBase() {}
    //! compute forces on the dynamical state type and update the accel attribute(s)
    virtual void Compute(StateType& s, double dt) = 0;
    virtual ~ForceBase() {};
};
template<typename StateType>
using ForcePtr = std::unique_ptr<ForceBase<StateType>>;
using ForceDynamicsPtr = std::unique_ptr<ForceBase<DynamicalStateData>>;
using ForceSPHPtr = std::unique_ptr<ForceBase<SPHStateData>>;
using ForceSBDPtr = std::unique_ptr<ForceBase<SoftBodyStateData>>;
using ForceRBDPtr = std::unique_ptr<ForceBase<RigidBodyStateData>>;

} //end of pba namesapce


#endif //__PBA_FORCE_H__