//*******************************************************************
//
//   NeighborSearch.h
//
//   Library of forces for different states.
//
//
//
//*******************************************************************

#ifndef ____PBA_NEIGHBORSEARCH_H____
#define ____PBA_NEIGHBORSEARCH_H____

#include "AABB.h"
#include "DynamicalState.h"

#include <vector>

namespace pba
{
/*!
  Neighorhood search in a given radius for SPH, return adjacent cells based on kernel radius
 */
class NeighborSearch
{
  public:
    NeighborSearch();
    //! bounds of grid, SPH kernel radius
    NeighborSearch(const AABB& b, float radius);
    ~NeighborSearch();
    //! populate grid with particles
    void Populate(const DynamicalStateData& state);
    void ClearGrid();

    size_t Index(const Vector& P) const;
    size_t Index(size_t i, size_t j, size_t k) const;
    void AntiIndex(size_t ind, size_t& i, size_t& j, size_t& k) const;
    //! populate neighbors vector 
    void NeighborsList(std::vector<size_t>& neighbors, const Vector& pos, bool use_parallel);
    //! check if particle is inside grid
    void InBounds(const Vector& pos, size_t i);
    //! size the grid based on dims
    void ComputeSize();
    void set_cellsize(float r);
    void set_bounds(const AABB& b);

  private:
    AABB bounds_;
    float radius_;
    std::vector< std::vector<size_t> > grid_;
    float cell_size_;
    size_t size_;
    Vector L_;
    size_t nx_;
    size_t ny_;
    size_t nz_;

};


} //end of pba namespace

#endif //____PBA_NEIGHBORSEARCH_H____