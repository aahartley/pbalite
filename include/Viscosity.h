#ifndef __PBA_VISCOSITY_H__
#define __PBA_VISCOSITY_H__

#include "DynamicalState.h"
#include "SPHState.h"
#include "SoftBodyState.h"
#include <memory>
#include "Force.h"
namespace pba
{
class ExplicitViscosity : public ForceBase
{
  public:
    ExplicitViscosity(const float kv):
      dynamic_viscosity (kv)
    {}
    ~ExplicitViscosity(){}
    void compute(DynamicalState& s, const double dt);
    void compute(SPHState& s, const double dt);
    void compute(SoftBodyState& s, const double dt);
    void set_visc(const float kv){dynamic_viscosity=kv;}
    const float get_visc() const{return dynamic_viscosity;}
  private:
    float dynamic_viscosity; 
};
pba::Force CreateExplicitViscosity(const float kv);

}

#endif