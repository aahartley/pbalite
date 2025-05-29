#include "RigidBodyState.h"

using namespace pba;

RigidBodyStateData::RigidBodyStateData( const std::string& nam) :
  DynamicalStateData(nam+"RigidBodyStateData")
{

}

RigidBodyStateData::~RigidBodyStateData(){}

void RigidBodyStateData::compute_RBD_data() // Initialize COM, total mass, home state
{
    center_of_mass = Vector(0,0,0);
    total_mass = 0.0;
    for( size_t p=0;p<nb();p++)
    {
        total_mass += mass(p);
        center_of_mass += pos(p)*mass(p);
    }
    center_of_mass /= total_mass;
    for( size_t p=0;p<nb();p++)
    {
        Vector d = pos(p)-center_of_mass;
        set_pos(p,d);
    }
}
//! position of a particle in the current, rotated state
Vector RigidBodyStateData::vert_pos( const size_t p ) const
{
    Vector P = angular_rotation*pos(p) + center_of_mass;
    return P;
}

RigidBodyState pba::CreateRigidBodyState( const std::string& nam)
{
    return std::make_shared<RigidBodyStateData>(nam);
}