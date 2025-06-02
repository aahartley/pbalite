//*******************************************************************
//
//   SoftBodyState.h
//
//  Bounding Volume Hierarchy for continous collison detection
//
//
//
//*******************************************************************

#ifndef ____PBA_SOFT_BODY_STATE_H____
#define ____PBA_SOFT_BODY_STATE_H____

#include "DynamicalState.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <memory>

namespace pba
{

class SoftEdge
{
  public:
    SoftEdge() {}
    SoftEdge(size_t i, size_t j, double edgelength) 
      : inode_(i),
        jnode_(j),
        edge_length_(edgelength) {}
    ~SoftEdge() {}
    size_t get_first_node() const { return inode_; }
    size_t get_second_node() const { return jnode_; }
    double get_edge_length() const { return edge_length_; }
  private:
    size_t inode_, jnode_; //indices of the two particles (i.e., a and b)
    double edge_length_; // Lab
};





class SoftBodyStateData : public DynamicalStateData
{
  public:
    SoftBodyStateData(const std::string& name = "SoftBodyDataNoName");
    ~SoftBodyStateData();
    // copy constructor and assignment
    // SoftBodyStateData( const SoftBodyStateData& d );
    // SoftBodyStateData& operator= ( const SoftBodyStateData& d );

    void set_num_pairs(size_t n) { connected_pairs_.resize(n + nb_pairs()); }
    const SoftEdge& get_connected_pair(size_t p) const { return connected_pairs_[p]; }
    size_t nb_pairs() const { return connected_pairs_.size(); }

    bool Empty() const { return connected_pairs_.empty(); }
    void ClearPairs() { connected_pairs_.clear(); }
    void AddPair(size_t i, size_t j, size_t ind);

  private:
    std::vector<SoftEdge> connected_pairs_;
};
using SoftBodyStatePtr = std::unique_ptr<SoftBodyStateData>;
inline SoftBodyStatePtr CreateSoftBodyState(const std::string& name="SoftBodyDataNoName")
{
  return std::make_unique<SoftBodyStateData>(name);
}





}

#endif