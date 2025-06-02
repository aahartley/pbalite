#include "SoftBodyState.h"

namespace pba
{


SoftBodyStateData::SoftBodyStateData(const std::string& name) 
  : DynamicalStateData(name+"SoftBodyStateData")
{
    create_attr("pin", false);

}

SoftBodyStateData::~SoftBodyStateData(){}

void SoftBodyStateData::AddPair(size_t i, size_t j, size_t ind)
{
    double l_ab = (pos(j) - pos(i)).magnitude();
    //if(L_ab > 0.2) return;
    SoftEdge edge(i, j, l_ab);
    connected_pairs_[ind] = edge;
}


//For N vertices, there will be N(N − 1)/2 edges.

}