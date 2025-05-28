
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
void CreateInfiniteCube(pba::CollisionSurface& s, double x, double y, double z );
void CreateCube(pba::CollisionSurface& s, double x, double y, double z);

void DisplayInfinitePlanes( pba::CollisionSurface& s );
void Display( pba::CollisionSurface& s );


void GeoToSoftBody(std::vector<Vector>& verts, std::vector<Triangle>& tris, const std::string& file, SoftBodyState& s );




};


#endif


