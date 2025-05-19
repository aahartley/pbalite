
#ifndef __PBA_PBAUTILS_H__
#define __PBA_PBAUTILS_H__

#include "CollisionSurface.h"
#include "PbaThing.h"
#include "SoftBodyState.h"
#include "ObjLoader.h"
#include <iostream>


using namespace std;

namespace pba{


// void AddCollisionSurface( pba::CollisionSurface& s, pba::PbaThing& p ); 

void Display( pba::CollisionSurface& s );

void GeoToSoftBody(const std::string& file, SoftBodyState& s );




};


#endif


