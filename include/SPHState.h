//*******************************************************************
//
//   SPHState.h
//
//   State attributes for SPH.
//
//
//
//*******************************************************************

//(weakly incompressible & DFSPH).
//  DFSPH functions are based off the 2015 paper: https://animation.rwth-aachen.de/media/papers/2015-SCA-DFSPH.pdf
//  2017 paper, where the equations are written differently (still equivalent): https://animation.rwth-aachen.de/media/papers/2017-TVCG-ViscousDFSPH.pdf

#ifndef ____PBA_SPHSTATE_H____
#define ____PBA_SPHSTATE_H____

#include "DynamicalState.h"
#include "NeighborSearch.h"

namespace pba
{
/*!
  SPH state attribs, density functions, neighbor search
 */
class SPHStateData : public DynamicalStateData, public NeighborSearch
{
  public:
    SPHStateData(const AABB& bounds, double h, const std::string& name = "SPHDataNoName" );
    ~SPHStateData();

    float Weight(size_t p, const Vector& P) const;
    Vector GradWeight(size_t p, const Vector& P) const;

    void ComputeDensity();
    void ComputePredictedDensity(size_t p, double dt);
    void ComputeDensityDerivative(size_t p);
    void ComputeFactor();
    void Populate();

    // float average_density(); //const;
    // float average_density_derivative(); //const;
    // float average_predicted_density(); //const;
    float MaxVelocity() const;

    float get_radius() const { return radius_; }
    void set_radius(float v);

    float get_particle_radius() const { return particle_radius_; }
    float get_density0() const { return density0_; }
    void set_density0(const float d0);

    int get_max_iter() const { return max_iter_; }
    void set_max_iter(int maxi);

    bool get_use_user_dt() const { return use_user_dt_; }
    void set_use_user_dt(bool uudt);

    float get_meps() const { return m_eps_;}
    float get_m_max_error() const { return m_max_error_;}
    void set_m_max_error(int mmx);

    float get_max_error() const { return max_error_;}
    void set_max_error(int mx);

    bool get_dd_clamp() const { return dd_clamp_;}
    void set_dd_clamp(bool cl);

    bool get_neighbor_parallel() const { return neighbor_parallel_;}
    void set_neighbor_parallel(bool cl);

  private:
    float radius_;
    float particle_radius_;
    float density0_;
    float m_eps_;
    float max_error_; //divergence
    float m_max_error_; //density
    int max_iter_;
    bool dd_clamp_;
    bool use_user_dt_;
    bool neighbor_parallel_;

};

using SPHStatePtr =  std::unique_ptr<SPHStateData>;

inline SPHStatePtr CreateSPH(const AABB& bounds, double h, const std::string& name = "SPHDataNoName")
{
  return std::make_unique<SPHStateData>(bounds, h, name);
}


} //end of pba namespace

#endif //____PBA_SPHSTATE_H____