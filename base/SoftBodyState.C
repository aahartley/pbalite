#include "SoftBodyState.h"

using namespace pba;

SoftEdge pba::CreateSoftEdge( const size_t i, const size_t j, const double edgelength )
{
    return std::make_shared<SoftEdgeData>(i, j, edgelength);
}


SoftBodyStateData::SoftBodyStateData( const std::string& nam) :
  DynamicalStateData(nam+"SoftBodyStateData")
{

}

SoftBodyStateData::~SoftBodyStateData(){}

void SoftBodyStateData::add_pair(size_t i, size_t j, size_t ind)
{
    double L_ab = (pos(j) - pos(i)).magnitude();
    //if(L_ab > 0.2) return;
    SoftEdge edge = CreateSoftEdge(i, j, L_ab);
    connected_pairs[ind] = edge;
}

SoftBodyState pba::CreateSoftBodyState( const std::string& nam)
{
    return std::make_shared<SoftBodyStateData>(nam);
}

//For N vertices, there will be N(N − 1)/2 edges.