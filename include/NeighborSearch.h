
#ifndef ____PBA_NEIGHBORSEARCH_H____
#define ____PBA_NEIGHBORSEARCH_H____

#include <vector>

#include "AABB.h"
#include "DynamicalState.h"

namespace pba
{

class NeighborSearch
{
  public:
    NeighborSearch();
    NeighborSearch(const AABB& b, float r);
    ~NeighborSearch();

    void populate(const DynamicalStateData& state);
    void clear_grid();

    size_t index( const Vector& P ) const;
    size_t index( size_t i, size_t j, size_t k ) const;
    void anti_index( const size_t ind, size_t& i, size_t& j, size_t& k ) const;


    void neighbors_list(std::vector<size_t>& neighbors, const Vector& pos, const bool use_parallel);


    void in_bounds(const Vector& pos, size_t i);
    void compute_size();
    void set_cellsize(const float r);
    void set_bounds(const AABB& b);



  private:
    AABB bounds;
    float radius;
    std::vector< std::vector<size_t> > grid;
    float cell_size;
    size_t size;
    Vector L;
    size_t nx;
    size_t ny;
    size_t nz;



};





}



#endif