#include "NeighborSearch.h"

using namespace pba;

NeighborSearch::NeighborSearch(){}

NeighborSearch::NeighborSearch(const AABB& b, float r): bounds(b), radius(r)
{
    cell_size = 2 * radius;
    L = bounds.getURC() - bounds.getLLC();
    nx = L.X()/cell_size + 1; 
    ny = L.Y()/cell_size + 1; 
    nz = L.Z()/cell_size + 1; 
    size = nx * ny * nz;
    grid = std::vector< std::vector<size_t> >(size);

}

void NeighborSearch::compute_size()
{
    nx = L.X()/cell_size + 1; 
    ny = L.Y()/cell_size + 1; 
    nz = L.Z()/cell_size + 1; 
    size = nx * ny * nz;
    grid.clear();
    grid.resize(size);
}

void NeighborSearch::set_cellsize(const float c)
{
    cell_size = c;
    compute_size();
}

void NeighborSearch::set_bounds(const AABB& b)
{
    bounds = b;
    L = bounds.getURC() - bounds.getLLC();
    compute_size();
}

NeighborSearch::~NeighborSearch(){}

void NeighborSearch::populate(const DynamicalStateData& state)
{   
    clear_grid();
    #pragma omp parallel for
    for(size_t i = 0; i < state.nb(); i++)
    {
        in_bounds(state.pos(i), i);
    }

}

void NeighborSearch::neighbors_list(std::vector<size_t>& neighbors, const Vector& pos, const bool use_parallel=false)
{
    //todo: pre accolator mem on avg cell size?
    size_t grid_ind = index(pos);
    if(grid_ind < size)
    {
    size_t i,j,k;
    anti_index(grid_ind, i, j ,k);
    if(use_parallel)
    {
        #pragma omp parallel
        {
        std::vector<size_t> thread_neighbors;

        #pragma omp for collapse(3) nowait
        for(size_t kk = k-1; kk <= k+1; kk++ )
        {
            for(size_t jj = j-1; jj <= j+1; jj++ )
            {
                for(size_t ii = i-1; ii <= i+1; ii++ )
                {
                    if (ii < 0 || jj < 0 || kk < 0) continue;
                    size_t ind = index(ii,jj,kk);
                    if( (ind < size ))
                        thread_neighbors.insert(thread_neighbors.end(), grid[ind].begin(), grid[ind].end());
                }

            }

        }
        #pragma omp critical
            neighbors.insert(neighbors.end(), thread_neighbors.begin(), thread_neighbors.end());
        }
    }
    else
    {
        for (int kk = (int)k - 1; kk <= (int)k + 1; kk++)
        {
            for (int jj = (int)j - 1; jj <= (int)j + 1; jj++) 
            {
                for (int ii = (int)i - 1; ii <= (int)i + 1; ii++) 
                {
                    if (ii < 0 || jj < 0 || kk < 0) continue;
                    size_t ind = index(ii, jj, kk);
                    if (ind < size) 
                        neighbors.insert(neighbors.end(), grid[ind].begin(), grid[ind].end());
                    
                }
            }
        }
    }
    }
}


void NeighborSearch::in_bounds(const Vector& pos, size_t p)
{
    if(bounds.isInside(pos))
    {
        size_t ind = index(pos);
        if(ind < size)
        {
            #pragma omp critical
            {
            grid[ind].push_back(p);
            }
        }
    }
}

void NeighborSearch::clear_grid()
{
    #pragma omp parallel for
    for (size_t i = 0; i < grid.size(); ++i)
    {
        grid[i].clear();
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
    for(size_t q=0;q<3;q++)
    {
       if(rel_pos[q] < 0.0 ){ return nx*ny*nz;}
    }
    Vector grid_pos = rel_pos/cell_size;
    for(size_t q=0;q<3;q++)
    {
       if(grid_pos[q] < 0.0 ){ return nx*ny*nz;}
    }
    size_t i = (size_t)grid_pos.X();
    size_t j = (size_t)grid_pos.Y();
    size_t k = (size_t)grid_pos.Z(); 
    return index(i,j,k);
}