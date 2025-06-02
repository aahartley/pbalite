//*******************************************************************
//
//   ExplicitDynamics.C
//
//   Partial Solvers for pos and vel.
//
//
//
//*******************************************************************

#include "ExplicitDynamics.h"
#include "GISolver.h"

#include <iostream>

namespace pba
{
void AdvancePosition::Init(){}
void AdvancePosition::Solve(double dt)
{
    #pragma omp parallel for
    for (size_t i = 0; i < pq_.nb(); ++i)
    {
      pq_.set_pos(i, pq_.pos(i) + pq_.vel(i)*dt);
    }
    if (use_coll_ && coll_handler_)
    {
      coll_handler_->HandleCollisions(dt, pq_);
    } 
    pq_.set_stage(POSITION_UPDATED);

}

template<typename StateType>
AdvanceVelocity<StateType>::AdvanceVelocity(StateType& pq, ForceBase<StateType>& f, float a, float v)
  : pq_(pq),
    force_(f),
    acceleration_clamp_(a),
    velocity_clamp_(v) {}

template<typename StateType>
void AdvanceVelocity<StateType>::Solve(double dt)
{

  force_.Compute(pq_, dt); // computes the force and stores (force/mass) in the state vector acceleration member
  #pragma omp parallel for
  for (size_t i = 0; i < pq_.nb(); ++i)
  {    
    Vector acc = pq_.accel(i);
    float a_mag = acc.magnitude();
    if (a_mag > acceleration_clamp_)
    {
      acc *= acceleration_clamp_ / a_mag;
    }
    Vector vel = pq_.vel(i) + acc*dt;
    float v_mag = vel.magnitude();
    if (v_mag > velocity_clamp_)
    {
      vel *= velocity_clamp_ / v_mag;
    }
    pq_.set_vel(i, vel);
  }
  pq_.set_stage(VELOCITY_UPDATED);

}

template<typename StateType>
GISolverPtr pba::CreateAdvanceVelocity(StateType& pq, ForceBase<StateType>& f, float a, float v)
{
  return std::make_unique<AdvanceVelocity<StateType>>(pq, f, a, v);

}

template GISolverPtr CreateAdvanceVelocity<DynamicalStateData>(
  DynamicalStateData&, ForceBase<DynamicalStateData>&, float, float);
template GISolverPtr CreateAdvanceVelocity<SPHStateData>(
  SPHStateData&, ForceBase<SPHStateData>&, float, float);
template GISolverPtr CreateAdvanceVelocity<SoftBodyStateData>(
  SoftBodyStateData&, ForceBase<SoftBodyStateData>&, float, float);
template GISolverPtr CreateAdvanceVelocity<RigidBodyStateData>(
  RigidBodyStateData&, ForceBase<RigidBodyStateData>&, float, float);

// Explicit instantiation
template class AdvanceVelocity<DynamicalStateData>;
template class AdvanceVelocity<SPHStateData>;
template class AdvanceVelocity<SoftBodyStateData>;
template class AdvanceVelocity<RigidBodyStateData>;


}
