
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
    NeighborSearch(AABB& b, float r);
    ~NeighborSearch();

    void populate(const DynamicalStateData& state);
    void clear();

    size_t index( const Vector& P ) const;
    size_t index( size_t i, size_t j, size_t k ) const;
    void anti_index( const size_t ind, size_t& i, size_t& j, size_t& k ) const;


    void neighborsList(std::vector<size_t>& neighbors, const Vector& pos);


    void inBounds(const Vector& pos, size_t i);



  private:
    AABB bounds;
    float radius;
    std::vector< std::vector<size_t> > grid;
    float cell_size;
    int size;
    Vector L;
    size_t nx;
    size_t ny;
    size_t nz;



};





}



#endif