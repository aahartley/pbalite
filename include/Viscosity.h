//*******************************************************************
//
//   Viscosity.h
//
//  Viscosities for SPH
//
//
//
//*******************************************************************

#ifndef __PBA_VISCOSITY_H__
#define __PBA_VISCOSITY_H__

#include "SPHState.h"
#include "Force.h"

#include <memory>

namespace pba
{
//SPH SURVERY 2019 & 2022  
class ExplicitViscosity : public ForceBase<SPHStateData>
{
  public:
    ExplicitViscosity(float dv)
      : dynamic_viscosity_(dv) {}
    ~ExplicitViscosity(){}
    void Compute(SPHStateData& s, double dt);

    void set_visc(float kv)  {dynamic_viscosity_ = kv; }
    float get_visc() const { return dynamic_viscosity_; }
  private:
    float dynamic_viscosity_; 
};

inline ForceSPHPtr CreateExplicitViscosity(float dv)
{
  return std::make_unique<ExplicitViscosity>(dv);
}

} //end of pba namespace

#endif //__PBA_VISCOSITY_H__