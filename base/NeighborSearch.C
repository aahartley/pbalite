#include "NeighborSearch.h"

using namespace pba;

NeighborSearch::NeighborSearch(){}

NeighborSearch::NeighborSearch(AABB& b, float r): bounds(b), radius(r)
{
    cell_size = 2 * radius;
    L = bounds.getURC() - bounds.getLLC();
    nx = L.X()/cell_size + 1; 
    ny = L.Y()/cell_size + 1; 
    nz = L.Z()/cell_size + 1; 
    size = nx * ny * nz;
    grid = std::vector< std::vector<size_t> >(size);

}

NeighborSearch::~NeighborSearch(){}

void NeighborSearch::populate(const DynamicalStateData& state)
{   
    for(size_t i = 0; i < state.nb(); i++)
    {
        inBounds(state.pos(i), i);
    }

}

void NeighborSearch::neighborsList(std::vector<size_t>& neighbors, const Vector& pos)
{
    size_t grid_ind = index(pos);
    size_t i,j,k;
    anti_index(grid_ind, i, j ,k);
    for(size_t kk = k-1; kk <= k+1; kk++ )
    {
        for(size_t jj = j-1; jj <= j+1; jj++ )
        {
            for(size_t ii = i-1; ii <= i+1; ii++ )
            {
                size_t ind = index(ii,jj,kk);
                if( (int)ind < size )
                    neighbors.insert(neighbors.end(), grid[ind].begin(), grid[ind].end());
            }

        }

    }
}


void NeighborSearch::inBounds(const Vector& pos, size_t p)
{
    if(bounds.isInside(pos))
    {
        size_t ind = index(pos);
        grid[ind].push_back(p);
    }
}

void NeighborSearch::clear()
{
    for (std::vector<size_t>& cell : grid) 
    {
        cell.clear();
    }
}

void NeighborSearch::anti_index( const size_t ind, size_t& i, size_t& j, size_t& k ) const
{
   k = ind/(nx*ny);
   j = (ind-k*nx*ny)/nx;
   i = (ind-k*nx*ny-j*nx);
}

size_t NeighborSearch::index( size_t i, size_t j, size_t k ) const
{
   if(i>=nx){ return nx*ny*nz; }
   if(j>=ny){ return nx*ny*nz; }
   if(k>=nz){ return nx*ny*nz; }
   return i + nx*(j+ny*k);
}

size_t NeighborSearch::index( const Vector& pos ) const
{
    Vector rel_pos = pos - bounds.getLLC();
    Vector grid_pos = rel_pos/cell_size;
    size_t i = (size_t)grid_pos.X();
    size_t j = (size_t)grid_pos.Y();
    size_t k = (size_t)grid_pos.Z(); 
    return index(i,j,k);
}