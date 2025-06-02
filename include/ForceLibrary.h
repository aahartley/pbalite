//*******************************************************************
//
//   ForceLibrary.h
//
//   Library of forces for different states.
//
//
//
//*******************************************************************

#ifndef __PBA_FORCE_LIBRARY_H__
#define __PBA_FORCE_LIBRARY_H__

#include "Force.h"
#include "DynamicalState.h"
#include "SPHState.h"
#include "SoftBodyState.h"
#include "RigidBodyState.h"

#include <memory>

namespace pba
{

/*!
  Gravity Force for all states
 */
template<typename StateType>
class GravityForce : public ForceBase<StateType>
{
  public:
    GravityForce(const Vector& gravity)
      : gravity_(gravity) {}
    ~GravityForce() {}
    void Compute(StateType& s, double dt);

    void set_gravity(const Vector& gravity) { gravity_ = gravity; }
    const Vector& get_gravity() const { return gravity_; }
  private:
    Vector gravity_; 
};

// add more statetypes if needed
/*!
  Tait pressure Force for SPH
 */
class TaitPressureForce : public ForceBase<SPHStateData>
{
  public:
    TaitPressureForce(double strength, double rest_dens, double gamma)
      : strength_(strength),
        rho_0_(rest_dens),
        gamma_(gamma) {}
    ~TaitPressureForce(){}
    void Compute(SPHStateData& s, double dt);

    void set_strenth(double s) { strength_ = s; }
    double get_strength() const { return strength_; }

    void set_rho0(double s) { rho_0_ = s; }
    double get_rho0() const { return rho_0_; }

    void set_gamma(double s) { gamma_ = s; }
    double get_gamma() const { return gamma_; }
  private:
    double strength_;
    double rho_0_;
    double gamma_;
};

/*!
  Harmonic Oscillator Force for all states
 */
template<typename StateType>
class HarmonicOscillatorForce : public ForceBase<StateType>
{
  public:
    HarmonicOscillatorForce(double kd)
      : kd_(kd) {}
    ~HarmonicOscillatorForce(){}
    void Compute(StateType& s, double dt);

    void set_kd(double v){ kd_ = v; }
    double get_kd() const { return kd_; }
  private:
    double kd_;
};

/*!
  Stores vector of forces for all states
 */
template<typename StateType>
class AccumulatingForce : public ForceBase<StateType>
{
  public:
    AccumulatingForce(){}
    ~AccumulatingForce(){}
    void Compute(StateType& s, double dt);

    //! Build up the collection of forces to accumulate
    void AddForce(ForceBase<StateType>* f);
  private:
    std::vector<ForceBase<StateType>*> forces;
};

/*!
  Strut Force for SoftBody (pairs of verts)
 */
class AccumulatingStrutForce : public pba::ForceBase<SoftBodyStateData>
{
  public:
    AccumulatingStrutForce(double spring, double friction, bool crit_damp) 
    : spring_(spring),
      friction_(friction),
      crit_damp_(crit_damp) {}
    ~AccumulatingStrutForce() {}
    void Compute(SoftBodyStateData& pq, double dt);
  private:
    double spring_;
    double friction_;
    bool crit_damp_;
};

template<typename StateType>
ForcePtr<StateType> CreateGravityForce(const Vector& g);

ForceSPHPtr CreateTaitPressureForce(float strength, float rest_dens, float gamma);

template<typename StateType>
ForcePtr<StateType> CreateHarmonicOscillatorForce(double k);

template<typename StateType>
ForcePtr<StateType> CreateAccumulatingForce();

ForceSBDPtr CreateAccumulatingStrutForce(double spring, double friction, bool crit_damp);

}

#endif