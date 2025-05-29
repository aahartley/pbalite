#ifndef ____PBA_RIGID_BODY_STATE_H____
#define ____PBA_RIGID_BODY_STATE_H____

#include "DynamicalState.h"
#include "Matrix.h"
#include "LinearAlgebra.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

namespace pba
{


class RigidBodyStateData : public DynamicalStateData
{
  public:
    RigidBodyStateData( const std::string& nam = "RBDDataNoName");
    ~RigidBodyStateData();
    void compute_RBD_data(); // Initialize COM, total mass, home state
    //! position of a particle in the current, rotated state
    Vector vert_pos( const size_t p ) const;
    Vector center_of_mass; // center of mass position
    Matrix angular_rotation; // rotation matrix R
  private:
    float total_mass;
};
typedef std::shared_ptr<RigidBodyStateData> RigidBodyState;
RigidBodyState CreateRigidBodyState( const std::string& nam);

}
#endif