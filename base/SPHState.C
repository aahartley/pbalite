//*******************************************************************
//
//   SPHState.C
//
//   State attributes for SPH.
//
//
//
//*******************************************************************

#include "SPHState.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

namespace pba
{

SPHStateData::SPHStateData(const AABB& bounds, const double h, const std::string& name) 
  : DynamicalStateData(name+"SPHStateData"),
    NeighborSearch(bounds, h),
    radius_(h),
    particle_radius_(radius_/4.0f),
    density0_(1000.f),
    m_eps_(1.0e-5f),
    max_error_(0.1f),
    m_max_error_(0.01f),
    max_iter_(100),
    dd_clamp_(true),
    use_user_dt_(false),
    neighbor_parallel_(false)
{
    float particle_diamter = particle_radius_ * 2;
    float mV = 0.8f * particle_diamter * particle_diamter * particle_diamter;
 
    float default_value = 0.0f;
    Vector default_vector(0.f, 0.f, 0.f);
 
    create_attr("density", default_value);
    create_attr("predicted_density", default_value);
    create_attr("density_derivative", default_value);
    create_attr("pressure", default_value );
    create_attr("divergence", default_value );
    create_attr("factor", default_value);
    create_attr("volume", mV);
    create_attr("k_i", default_value);
    create_attr("kv_i", default_value);
    create_attr("pressure_acc", default_vector);
}

SPHStateData::~SPHStateData(){}

//cubic spline kernel
float SPHStateData::Weight(const size_t p, const Vector& position) const
{
    float res = 0.0;
    Vector dist = (position - pos(p));
    float x = dist.magnitude();
    float h3 = radius_ * radius_ * radius_;
    float m_k = static_cast<float>(8.0) / (M_PI*h3);
    const float q = x / radius_;
    if (q <= 1.0)
    {
        if (q <= 0.5)
        {
            const float q2 = q*q;
            const float q3 = q2*q;
            res = m_k * (static_cast<float>(6.0)*q3 - static_cast<float>(6.0)*q2 + static_cast<float>(1.0));
        }
        else
        {
            res = m_k * (static_cast<float>(2.0)*pow(static_cast<float>(1.0) - q, static_cast<float>(3.0)));
        }
    }
    //std::cout << "Weight "<< p << ": " << res << '\n';
    return res;
}

const Vector SPHStateData::GradWeight(const size_t p, const Vector& position) const
{
    Vector res; 
    Vector dist = (position- pos(p));
    const float x = dist.magnitude();
    const float q = x / radius_;
    float h3 = radius_ * radius_ * radius_;
    float m_l = static_cast<float>(48.0) / (M_PI*h3);
    if ((x > 1.0e-9) && (q <= 1.0))
    {
        Vector gradq = dist / x;
        gradq /= radius_;
        if (q <= 0.5)
        {
            res = m_l*q*((float) 3.0*q - static_cast<float>(2.0))*gradq;
        }
        else
        {
            const float factor = static_cast<float>(1.0) - q;
            res = m_l*(-factor*factor)*gradq;
        }
    }
    else
        res.set(0,0,0);

    return res;
}

void SPHStateData::ComputeDensity()
{
    //file.open("./densf/dens"+std::to_string(densiter)+".txt");
    //file << densiter << '\n';
    #pragma omp parallel for
    for (size_t p = 0; p <nb(); ++p)
    {
        float density = 0.0;
        const Vector& position = pos(p);
        std::vector<size_t> neighbors;
        neighbors_list(neighbors, position, neighbor_parallel_);
        //std::cout << p << ": " <<cells.size() << '\n';
        for (size_t a = 0; a < neighbors.size(); a++)
        {
            size_t pid = neighbors[a]; 

            density += get_float_attr("volume", pid) * Weight(pid, position);
            //density += mass(pid) * Weight(pid,P);
        }

        density *= density0_;
        //#pragma omp critical
        // file << density << "\n";
        set_attr("density", p, density);
        // TODO: DEBUG FLAG
        //#pragma omp critical
        //   {
        //   if(density == 0) std::cout <<"cells size: "<<neighbors.size() << " density0: " << density0<<" particle " << p << ": density is 0\n";
        //   if(std::isnan(density) || std::isinf(density)) std::cout << "density bad\n";
        //   }
    } 
    //densiter++;
    // file.close();
}

//DFSPH 2015
//Equation above (12)
//ρ∗i = ρi + ∆t Dρi/Dt = ρi + ∆t ∑j mj (v∗i − v∗j )∇Wi
void SPHStateData::ComputePredictedDensity(const size_t p, const double dt)
{

    float density = get_float_attr("density", p);
    float pdensity = 0.0f;
    const Vector& position = pos(p);
    const Vector& velocity = vel(p);
    std::vector<size_t> neighbors;
    neighbors_list(neighbors, position, neighbor_parallel_);

    for (size_t a = 0; a < neighbors.size(); ++a)
    {
        size_t pid = neighbors[a]; 

        const Vector& v_j = vel(pid);
        pdensity += ((velocity - v_j) * GradWeight(pid, position));
    }

    pdensity *= get_float_attr("volume", p) * density0_;
    //pdensity *= mass(p);

    //pdensity = (density / density0) + dt *pdensity;
    pdensity = (density / 1) + dt *pdensity;

    // if(pdensity <= 0)
    // {
    //    #pragma omp critical
    //    //std::cout << "pdens: " << pdensity << " pdens*dt: " << pdensity*dt << '\n';
    // }
    //if(std::isnan(pdensity) || std::isinf(pdensity)) std::cout << "pdensity bad\n";
    set_attr("predicted_density", p, pdensity);
   
}

//DFSPH 2015
//Equation (9)
//Dρi/Dt = ∑j mj (vi − vj )∇Wi
void SPHStateData::ComputeDensityDerivative(const size_t p)
{
    float density_change = 0.0;
    const Vector& position = pos(p);
    const Vector& velocity = vel(p);
    std::vector<size_t> neighbors;
    neighbors_list(neighbors, position, neighbor_parallel_);

    for (size_t a = 0; a < neighbors.size(); ++a)
    {
        size_t pid = neighbors[a]; 

        const Vector& v_j = vel(pid);
        density_change += ((velocity - v_j) * GradWeight(pid, position));
    }
    density_change *= mass(p) ; //all fluid have constant volume

    // density_change = std::max(density_change, 0.0f);//mannually do this inside solver
    //#pragma omp critical
    //std::cout << "Dens deriv: " << density_change << '\n';
    if (std::isnan(density_change) || std::isinf(density_change)) std::cout << "densitychange bad\n";
    set_attr("density_derivative", p,  density_change);
}

//DFSPH 2015
//Equation (9)
//αi = ρi / |∑j mj ∇Wi|^2 + ∑j |mj ∇Wi|^2
//since alpha is used in ∇p, ρi is cancalled out, i.e. 1 / |∑j mj ∇Wi|^2 + ∑j |mj ∇Wi|^2
void SPHStateData::ComputeFactor()
{
    #pragma omp parallel for   
    for (size_t p = 0; p < nb(); ++p)
    {
        float factor = 0;
        float sum_grad_p = 0.0f; // sum of grad pi and pj

        Vector grad_p_i(0,0,0); //pressure gradient of ith particle
        const Vector& position = pos(p);
        std::vector<size_t> neighbors;
        neighbors_list(neighbors, position, neighbor_parallel_);

        for (size_t a = 0; a < neighbors.size(); ++a)
        {         
            size_t pid = neighbors[a]; 
       
            const Vector& grad_p_j = get_float_attr("volume",pid) * GradWeight(pid, position) * density0_;

            //dot product of a vector by itself is its magnitude squared
            sum_grad_p += grad_p_j * grad_p_j; 
            grad_p_i += grad_p_j;
        }
        sum_grad_p += grad_p_i * grad_p_i;

        if (sum_grad_p > m_eps_)
        {
            //factor  =  get_float_attr("density",p) / sum_grad_p;
            factor  =  1.0 / sum_grad_p;
        }
        else
            factor = 0.f;
        if (std::isnan(factor) || std::isinf(factor)) std::cout << "factor bad\n";
        set_attr("factor", p, factor);
    }
}

void SPHStateData::Populate()
{
   NeighborSearch::Populate(*this);
}

void SPHStateData::set_radius(const float h)
{
   radius_ = h;
   particle_radius_ = radius_/4.0f;
   float particle_diamter = particle_radius_ * 2;
   float m_v = 0.8f * particle_diamter * particle_diamter * particle_diamter;
   #pragma omp parallel for
   for(size_t p = 0; p < nb(); ++p)
   {
      set_attr("volume", p, m_v);
   }
   set_cellsize(2.0*radius_); 
}

void SPHStateData::set_max_iter(int maxi)
{
   if (maxi < 2) maxi = 2;
   max_iter_ = maxi;
}

void SPHStateData::set_use_user_dt(const bool uudt)
{
   use_user_dt_ = uudt;
}

void SPHStateData::set_m_max_error(int mmx)
{
   if (mmx < 0.01f ) mmx = 0.01f;
   m_max_error_ = mmx;
}

void SPHStateData::set_max_error(int mx)
{
   if (mx < 0.1f) mx = 0.1f;
   max_error_ = mx;
}

void SPHStateData::set_density0(float d0)
{
   if (d0 == 0) d0 =1;
   density0_ = d0;
}
void SPHStateData::set_dd_clamp(bool cl)
{
   dd_clamp_ = cl;
}
void SPHStateData::set_neighbor_parallel(bool np)
{
   neighbor_parallel_ = np;
}

float SPHStateData::MaxVelocity() const
{
    if (nb_items_ == 0) return 0;
    float vel_m = vel(0).magnitude();
    #pragma omp parallel for reduction(max:vel_m)
    for (size_t i = 1; i < nb(); ++i) {
        float curr_vel = vel(i).magnitude();
        // if (curr_vel > vel_m)
        //    vel_m = curr_vel;
        vel_m = std::fmax(vel_m, curr_vel);
    }
    return vel_m;
}

// //TODO: add pragmas
// float SPHStateData::average_density() //const
// {
//    //file.open("./densavgf/densavg"+std::to_string(densavgiter)+".txt");  
//    //file << densavgiter << '\n';
//    float density = 0.0;
//    for( size_t p=0;p<nb();p++)
//    {
//       density += get_float_attr("density",p);
//    }
//    density /= nb();
//    //file << density << '\n';
//    //file.close();
//    densavgiter++;
//    return density;
// }

// float SPHStateData::average_predicted_density() //const
// {
//    //file.open("./pdensavgf/pdensavg"+std::to_string(pdensavgiter)+".txt");  
//    //file << pdensavgiter << '\n';
//    float density = 0.0;
//    int count = 0;
//    for( size_t p=0;p<nb();p++)
//    {
//       float d = get_float_attr("predicted_density",p) * density0;

//       density += d;
//    }
//    density /= nb();
//   // file << density << '\n';
//    //file.close();
//    pdensavgiter++;
//    return density;
// }

// float SPHStateData::average_density_derivative() //const
// {
//    //file.open("./ddensavgf/ddensavg"+std::to_string(ddensavgiter)+".txt");  
//    //file << ddensavgiter << '\n';
//    float density = 0.0;
//    for( size_t p=0;p<nb();p++)
//    {
//       density += get_float_attr("density_derivative",p);
//    }
//    density /= nb();
//    //file << density << '\n';
//    //file.close();
//    ddensavgiter++;
//    return density;
// }

}//end of pba namespace