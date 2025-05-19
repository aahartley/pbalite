#ifndef ____PBA_SOFT_BODY_STATE_H____
#define ____PBA_SOFT_BODY_STATE_H____

#include "DynamicalState.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

namespace pba
{

class SoftEdgeData
{
  public:
    SoftEdgeData(const size_t i, const size_t j, const double edgelength) :
      inode(i),
      jnode(j),
      edge_length(edgelength)
    {}
    ~SoftEdgeData(){}
    const size_t& get_first_node() const { return inode; }
    const size_t& get_second_node() const { return jnode; }
    const double get_edge_length() const { return edge_length; }
  private:
    size_t inode, jnode; //indices of the two particles (i.e., a and b)
    double edge_length; // Lab
};
typedef std::shared_ptr<SoftEdgeData> SoftEdge;
SoftEdge CreateSoftEdge( const size_t i, const size_t j, const double edgelength );




class SoftBodyStateData : public DynamicalStateData
{
  public:
    SoftBodyStateData( const std::string& nam = "SoftBodyDataNoName");
    ~SoftBodyStateData();
    // copy constructor and assignment
    // SoftBodyStateData( const SoftBodyStateData& d );
    // SoftBodyStateData& operator= ( const SoftBodyStateData& d );

    void set_num_pairs(size_t n) {connected_pairs.resize(n+nb_pairs());}
    const SoftEdge& get_connected_pair( size_t p ) const { return connected_pairs[p]; }
    size_t nb_pairs() const { return connected_pairs.size(); }

    bool empty() const { return connected_pairs.empty(); }
    void clear_pairs() { connected_pairs.clear(); }
    void add_pair( size_t i, size_t j, size_t ind );

  private:
    std::vector<SoftEdge> connected_pairs;
};
typedef std::shared_ptr<SoftBodyStateData> SoftBodyState;
SoftBodyState CreateSoftBodyState( const std::string& nam="SoftBodyDataNoName");





}

#endif