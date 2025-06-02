//*******************************************************************
//
//   ExplicitDynamics.h
//
//   Partial Solvers for pos and vel.
//
//
//
//*******************************************************************

#ifndef __PBA_INTEGRATORS_H__
#define __PBA_INTEGRATORS_H__

#include "DynamicalState.h"
#include "SPHState.h"
#include "SoftBodyState.h"
#include "RigidBodyState.h"
#include "ForceLibrary.h"
#include "CollisionHandler.h"
#include "GISolver.h"
#include "Force.h"

namespace pba
{

class AdvancePosition : public GISolverBase
{
  public:
    AdvancePosition(DynamicalStateData& pq, CollisionHandlerBase* coll = nullptr)
      : pq_(pq), 
        coll_handler_(coll), 
        use_coll_(coll_handler_ != nullptr) {}
    ~AdvancePosition() {}
    void Init();
    void Solve(double dt);
  private:
    DynamicalStateData& pq_;
    CollisionHandlerBase* coll_handler_;
    bool use_coll_;
};

template<typename StateType>
class AdvanceVelocity : public GISolverBase
{
  public:
    AdvanceVelocity(StateType& pq, ForceBase<StateType>& f, float acceleration_clamp, float velocity_clamp);
    ~AdvanceVelocity() {}
    void Init() {}
    void Solve(double dt);

    float get_velocity_clamp() const { return velocity_clamp_; }
    void set_velocity_clamp(float v) { velocity_clamp_ = v; }

    float get_acceleration_clamp() const { return acceleration_clamp_; }
    void set_acceleration_clamp(float a) { acceleration_clamp_ = a; }
  private:
    StateType& pq_; 
    ForceBase<StateType>& force_; 
    float acceleration_clamp_;
    float velocity_clamp_;

};


template<typename StateType>
GISolverPtr CreateAdvanceVelocity(StateType& pq, ForceBase<StateType>& f, float a, float v);

inline GISolverPtr CreateAdvancePosition(DynamicalStateData& pq, CollisionHandlerBase* c_h = nullptr)
{
  return std::make_unique<AdvancePosition>(pq, c_h);

}



}

#endif
