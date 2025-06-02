//*******************************************************************
//
//   ForceLibrary.C
//
//   Library of forces for different states.
//
//
//
//*******************************************************************

#include "ForceLibrary.h"

namespace pba
{

template<typename StateType>
void GravityForce<StateType>::Compute(StateType& s, const double dt)
{
  #pragma omp parallel for
  for (size_t p = 0; p < s.nb(); ++p)
  {
      s.set_accel(p, s.accel(p) + gravity_);
  }
}

void TaitPressureForce::Compute(SPHStateData& s, const double dt)
{
  #pragma omp parallel for
  for (size_t p = 0; p < s.nb(); ++p)
  {
    Vector pressure;
    const Vector& pos = s.pos(p);
    std::vector<size_t> neighbors;
    s.neighbors_list(neighbors, pos, s.get_neighbor_parallel());
    float dens_a = s.get_float_attr("density", p);
    dens_a = std::fmax(dens_a, s.get_density0());
    for(size_t a = 0; a < neighbors.size(); ++a)
    {
      size_t pid = neighbors[a]; 
      float dens_b = s.get_float_attr("density", pid);
      dens_b = std::fmax(dens_b, s.get_density0());

      double pa = strength_ * (std::pow(dens_a/rho_0_, gamma_) - 1.0);
      double pb = strength_ * (std::pow(dens_b/rho_0_, gamma_) - 1.0);
      pressure += s.GradWeight(pid, pos) * s.mass(pid)
                * ((pa/(dens_a*dens_a)) + (pb/(dens_b*dens_b))) ;    
    }
    s.set_accel(p, s.accel(p) - pressure);
  }
}

template<typename StateType>
void HarmonicOscillatorForce<StateType>::Compute(StateType& s, const double dt)
{
    #pragma omp parallel for
    for(size_t p = 0; p < s.nb(); ++p)
    {
        s.set_accel(p, s.accel(p) - kd_ * s.pos(p) / s.mass(p));
    }
}

template<typename StateType>
void AccumulatingForce<StateType>::AddForce(ForceBase<StateType>* f) { forces.push_back(f); }

template<typename StateType>
void AccumulatingForce<StateType>::Compute(StateType& s, const double dt)
{
    #pragma omp parallel for
    for (size_t p = 0; p < s.nb(); ++p)
    {
        s.set_accel(p, Vector(0,0,0)); // initialize accelerations to zero before we start accumulating
    }
    for (size_t p = 0; p < forces.size(); ++p)
    {
        forces[p]->Compute(s, dt);
    }
}

void AccumulatingStrutForce::Compute(SoftBodyStateData& s, const double dt) 
{
  #pragma omp parallel for
  for (size_t i = 0; i < s.nb_pairs(); ++i)
  {
    const SoftEdge& se = s.get_connected_pair(i);
    size_t inode = se.get_first_node();
    size_t jnode = se.get_second_node();
    Vector d_ij = s.pos(jnode) - s.pos(inode);
    Vector v_ij = s.vel(jnode) - s.vel(inode);
    float r_mass = (s.mass(inode) * s.mass(jnode)) / (s.mass(inode) + s.mass(jnode));
    if(crit_damp_)
    {
      float critical_damping = 2 * std::sqrt(r_mass * spring_);
      friction_ = critical_damping;
    }
    Vector force;
    //both particles same spot
    //d_ij unnomrlaized is full error
    if (se.get_edge_length() < 0.000001 )
    {
      force = d_ij * spring_;
      force += v_ij * friction_;
    }
    else
    {
      double separation = d_ij.magnitude() - se.get_edge_length();
      if (d_ij.magnitude() > 0) d_ij.normalize();
      force = d_ij * (separation*spring_);
      force += d_ij * (d_ij*v_ij) * friction_;
    }
    #pragma omp critical
    {
      s.set_accel(jnode, s.accel(jnode) - force/s.mass(jnode));
      s.set_accel(inode, s.accel(inode) + force/s.mass(inode));
    }
  }
}

//TODO: refactor to inits in ExplicitDynamics (dont need extra functionalityy per state type)
// keeping as a reference

// GravityForce instantiations
template<>
ForceDynamicsPtr CreateGravityForce(const Vector& g) {
    return std::make_unique<GravityForce<DynamicalStateData>>(g);
}

template<>
ForceSPHPtr CreateGravityForce(const Vector& g) {
    return std::make_unique<GravityForce<SPHStateData>>(g);
}

template<>
ForceSBDPtr CreateGravityForce(const Vector& g) {
    return std::make_unique<GravityForce<SoftBodyStateData>>(g);
}

template<>
ForceRBDPtr CreateGravityForce(const Vector& g) {
    return std::make_unique<GravityForce<RigidBodyStateData>>(g);
}

// TaitPressureForce (non-template, already specialized for SPH)
ForceSPHPtr CreateTaitPressureForce(float strength, float rest_dens, float gamma) {
    return std::make_unique<TaitPressureForce>(strength, rest_dens, gamma);
}

// HarmonicOscillatorForce instantiations
template<>
ForceDynamicsPtr CreateHarmonicOscillatorForce(double k) {
    return std::make_unique<HarmonicOscillatorForce<DynamicalStateData>>(k);
}

template<>
ForceSPHPtr CreateHarmonicOscillatorForce(double k) {
    return std::make_unique<HarmonicOscillatorForce<SPHStateData>>(k);
}

ForceSBDPtr CreateHarmonicOscillatorForce(double k) {
    return std::make_unique<HarmonicOscillatorForce<SoftBodyStateData>>(k);
}

template<>
ForceRBDPtr CreateHarmonicOscillatorForce(double k) {
    return std::make_unique<HarmonicOscillatorForce<RigidBodyStateData>>(k);
}

// AccumulatingForce instantiations
template<>
ForceDynamicsPtr CreateAccumulatingForce() {
    return std::make_unique<AccumulatingForce<DynamicalStateData>>();
}

template<>
ForceSPHPtr CreateAccumulatingForce() {
    return std::make_unique<AccumulatingForce<SPHStateData>>();
}

template<>
ForceSBDPtr CreateAccumulatingForce() {
    return std::make_unique<AccumulatingForce<SoftBodyStateData>>();
}

template<>
ForceRBDPtr CreateAccumulatingForce() {
    return std::make_unique<AccumulatingForce<RigidBodyStateData>>();
}

// AccumulatingStrutForce (specialized for SoftBody)
ForceSBDPtr CreateAccumulatingStrutForce(double spring, double friction, bool crit_damp) {
    return std::make_unique<AccumulatingStrutForce>(spring, friction, crit_damp);
}

// Explicit instantiations
template class GravityForce<DynamicalStateData>;
template class GravityForce<SPHStateData>;
template class GravityForce<SoftBodyStateData>;
template class GravityForce<RigidBodyStateData>;

template class HarmonicOscillatorForce<DynamicalStateData>;
template class HarmonicOscillatorForce<SPHStateData>;
template class HarmonicOscillatorForce<SoftBodyStateData>;
template class HarmonicOscillatorForce<RigidBodyStateData>;

template class AccumulatingForce<DynamicalStateData>;
template class AccumulatingForce<SPHStateData>;
template class AccumulatingForce<SoftBodyStateData>;
template class AccumulatingForce<RigidBodyStateData>;

}//end of pba namespace
