//*******************************************************************
//
//   NeighborSearch.C
//
//   Library of forces for different states.
//
//
//
//*******************************************************************

#include "NeighborSearch.h"

using namespace pba;

NeighborSearch::NeighborSearch(){}

NeighborSearch::NeighborSearch(const AABB& b, const float radius)
  : bounds_(b),
    radius_(radius)
{
    cell_size_ = 2 * radius_;
    L_ = bounds_.urc() - bounds_.llc();
    nx_ = L_.X()/cell_size_ + 1; 
    ny_ = L_.Y()/cell_size_ + 1; 
    nz_ = L_.Z()/cell_size_ + 1; 
    size_ = nx_ * ny_ * nz_;
    grid_ = std::vector< std::vector<size_t> >(size_);

}

void NeighborSearch::ComputeSize()
{
    nx_ = L_.X()/cell_size_ + 1; 
    ny_ = L_.Y()/cell_size_ + 1; 
    nz_ = L_.Z()/cell_size_ + 1; 
    size_ = nx_ * ny_ * nz_;
    grid_.clear();
    grid_.resize(size_);
}

void NeighborSearch::set_cellsize(const float c)
{
    cell_size_ = c;
    ComputeSize();
}

void NeighborSearch::set_bounds(const AABB& b)
{
    bounds_ = b;
    L_ = bounds_.urc() - bounds_.llc();
    ComputeSize();
}

NeighborSearch::~NeighborSearch(){}

void NeighborSearch::Populate(const DynamicalStateData& state)
{   
    ClearGrid();
    #pragma omp parallel for
    for (size_t i = 0; i < state.nb(); ++i)
    {
        InBounds(state.pos(i), i);
    }

}

void NeighborSearch::NeighborsList(std::vector<size_t>& neighbors, const Vector& pos, const bool use_parallel=false)
{
    //todo: pre accolator mem on avg cell size_?
    size_t grid_ind = Index(pos);
    if (grid_ind < size_)
    {
    size_t i,j,k;
    
    AntiIndex(grid_ind, i, j ,k);
    if (use_parallel)
    {
        #pragma omp parallel
        {
            std::vector<size_t> thread_neighbors;

            #pragma omp for collapse(3) nowait
            for (size_t kk = k-1; kk <= k+1; ++kk)
            {
                for (size_t jj = j-1; jj <= j+1; ++jj)
                {
                    for (size_t ii = i-1; ii <= i+1; ++ii)
                    {
                        if (ii < 0 || jj < 0 || kk < 0) continue;
                        size_t ind = Index(ii,jj,kk);
                        if (ind < size_ )
                            thread_neighbors.insert(thread_neighbors.end(), grid_[ind].begin(), grid_[ind].end());
                    }
                }
            }
            #pragma omp critical
            {
                neighbors.insert(neighbors.end(), thread_neighbors.begin(), thread_neighbors.end());
            }
        }
    }
    else
    {
        for (size_t kk = k - 1; kk <= k + 1; ++kk)
        {
            for (size_t jj = j - 1; jj <= j + 1; ++jj) 
            {
                for (size_t ii = i - 1; ii <= i + 1; ++ii) 
                {
                    if (ii < 0 || jj < 0 || kk < 0) continue;
                    size_t ind = Index(ii, jj, kk);
                    if (ind < size_) 
                        neighbors.insert(neighbors.end(), grid_[ind].begin(), grid_[ind].end());  
                }
            }
        }
    }
    }
}


void NeighborSearch::InBounds(const Vector& pos, const size_t p)
{
    if (bounds_.IsInside(pos))
    {
        size_t ind = Index(pos);
        if (ind < size_)
        {
            #pragma omp critical
            {
                grid_[ind].push_back(p);
            }
        }
    }
}

void NeighborSearch::ClearGrid()
{
    #pragma omp parallel for
    for (size_t i = 0; i < grid_.size(); ++i)
    {
        grid_[i].clear();
    }
}

void NeighborSearch::AntiIndex(const size_t ind, size_t& i, size_t& j, size_t& k) const
{
   k = ind/(nx_*ny_);
   j = (ind-k*nx_*ny_)/nx_;
   i = (ind-k*nx_*ny_-j*nx_);
}

size_t NeighborSearch::Index(const size_t i, const size_t j, const size_t k) const
{
   if (i >= nx_) { return nx_*ny_*nz_; }
   if (j >= ny_) { return nx_*ny_*nz_; }
   if (k >= nz_) { return nx_*ny_*nz_; }
   return i + nx_*(j+ny_*k);
}

size_t NeighborSearch::Index(const Vector& pos) const
{
    Vector rel_pos = pos - bounds_.llc();
    for (size_t q = 0; q < 3; ++q)
    {
       if (rel_pos[q] < 0.0 ) { return nx_*ny_*nz_;}
    }
    Vector grid_pos = rel_pos/cell_size_;
    for (size_t q =0; q < 3; ++q)
    {
       if (grid_pos[q] < 0.0 ) { return nx_*ny_*nz_;}
    }
    size_t i = static_cast<size_t>(grid_pos.X());
    size_t j = static_cast<size_t>(grid_pos.Y());
    size_t k = static_cast<size_t>(grid_pos.Z());
    return Index(i,j,k);
}