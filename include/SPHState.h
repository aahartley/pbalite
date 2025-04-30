
//(weakly incompressible & DFSPH).
//  DFSPH functions are based off the 2015 paper: https://animation.rwth-aachen.de/media/papers/2015-SCA-DFSPH.pdf
//  2017 paper, where the equations are written differently (still equivalent): https://animation.rwth-aachen.de/media/papers/2017-TVCG-ViscousDFSPH.pdf


#ifndef ____PBA_SPHSTATE_H____
#define ____PBA_SPHSTATE_H____

#include "DynamicalState.h"
#include "NeighborSearch.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

namespace pba
{


class SPHStateData : public DynamicalStateData, public NeighborSearch
{
  public:

    SPHStateData( const AABB& bounds, const double h, const std::string& nam = "SPHDataNoName" );
   ~SPHStateData();


    const float get_radius() const { return radius; }
    void set_radius( const float& v );

    const float get_particle_radius() const { return particle_radius; }
    const float get_density0() const { return density0; }
    void set_density0(const float d0);
    const int get_maxIter() const { return maxIter; }
    void set_maxIter(int maxi);
    const bool get_useUserDT() const { return useUserDT; }
    void set_useUserDT(bool uudt);
    const float get_meps() const { return m_eps;}
    const float get_mMaxError() const { return m_maxError;}
    void set_mMaxError(int mmx);
    const float get_maxError() const { return maxError;}
    void set_maxError(int mx);
    const bool get_ddClamp() const { return dd_clamp;}
    void set_ddClamp(bool cl);

    const float weight( size_t p, const Vector& P ) const;
    const Vector grad_weight( size_t p, const Vector& P ) const;

    void compute_density();
    void compute_predicted_density(const size_t p, const double dt);
    void compute_density_derivative(size_t p);
    void compute_factor();
    void populate();

    float average_density(); //const;
    float average_density_derivative(); //const;
    float average_predicted_density(); //const;
    float max_velocity() const;

  
  private:

    float radius;
    float particle_radius;
    float density0;
    float m_eps;
    float maxError; //divergence
    float m_maxError; //density
    int maxIter;
    bool dd_clamp;
    bool useUserDT;



};



typedef std::shared_ptr<SPHStateData> SPHState;

SPHState CreateSPH( const AABB& bounds, const double h, const std::string& nam = "SPHDataNoName" );


}
#endif