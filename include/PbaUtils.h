
#ifndef __PBA_PBAUTILS_H__
#define __PBA_PBAUTILS_H__

#include "CollisionSurface.h"
#include "PbaThing.h"
#include "SoftBodyState.h"
#include "ObjLoader.h"


namespace pba
{

// void AddCollisionSurface( pba::CollisionSurface& s, pba::PbaThing& p ); 
void CreateInfiniteCube(CollisionSurface& s, double x, double y, double z);
void CreateCube(CollisionSurface& s, double x, double y, double z);

void DisplayInfinitePlanes(CollisionSurface* s);
void Display(CollisionSurface* s);


void GeoToSoftBody(std::vector<Vector>& verts, std::vector<Triangle>& tris, const std::string& file, SoftBodyStateData& s);

}

#endif


