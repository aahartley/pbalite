#include "SPHState.h"

using namespace pba;

SPHStateData::SPHStateData( const AABB& bounds, const double h, const std::string& nam ) :
DynamicalStateData(nam+"SPHStateData"),
NeighborSearch( bounds, h ),
radius(h),
particle_radius(radius/4.0f),
density0(1000.f),
m_eps(1.0e-5f),
maxError(0.1f),
m_maxError(0.01f),
maxIter(100),
dd_clamp(true),
useUserDT(false),
neighbor_parallel(false)
{
    float particle_diamter = particle_radius * 2;
    float mV = 0.8f * particle_diamter * particle_diamter * particle_diamter;
 
    float default_value = 0.0f;
    Vector default_vector(0.f,0.f,0.f);
 
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

SPHState pba::CreateSPH( const AABB& bounds, const double h, const std::string& nam )
{
   return std::make_shared<SPHStateData>(bounds, h, nam);
}

//cubic spline kernel
const float SPHStateData::weight(size_t p, const Vector& P) const
{
    float res = 0.0;
    Vector X = (P - pos(p));
    float x = X.magnitude();
    float h3 = radius * radius * radius;
    float m_k = static_cast<float>(8.0) / (M_PI*h3);
    const float q = x / radius;
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
//std::cout << "weight "<< p << ": " << res << '\n';
    return res;
}
const Vector SPHStateData::grad_weight(size_t p, const Vector& P) const
{
    Vector res; 
    Vector X = (P- pos(p));
    const float x = X.magnitude();
    const float q = x/ radius;
    float h3 = radius * radius * radius;
    float m_l = static_cast<float>(48.0) / (M_PI*h3);
    if ((x > 1.0e-9) && (q <= 1.0))
    {
        Vector gradq = X / x;
        gradq /= radius;
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

void SPHStateData::compute_density()
{

  // file.open("./densf/dens"+std::to_string(densiter)+".txt");
  // file << densiter << '\n';
    #pragma omp parallel for
    for( size_t p=0;p<nb();p++)
    {
      float density = 0.0;
      const Vector P = pos(p);
      std::vector<size_t> neighbors;
      neighbors_list(neighbors, P, neighbor_parallel);
   
      //std::cout << p << ": " <<cells.size() << '\n';
      for(size_t a = 0; a < neighbors.size(); a++)
      {
         size_t pid = neighbors[a]; 

         density += get_float_attr("volume", pid) * weight(pid,P);
         //density += mass(pid) * weight(pid,P);

         
      }

      density *= density0;
      //#pragma omp critical
     // file << density << "\n";
      set_attr("density", p, density);
      //todo:: DEBUG FLAG
    //   #pragma omp critical
    //   {
    //   if(density == 0) std::cout <<"cells size: "<<neighbors.size() << " density0: " << density0<<" particle " << p << ": density is 0\n";
    //   if(std::isnan(density) || std::isinf(density)) std::cout << "density bad\n";
    //   }

//#pragma omp critical
//      {
//	 std::cout << "Particle " << p << std::endl;
//	 std::cout << "\tposition " << P.X() << " " << P.Y() << " " << P.Z() << std::endl;
//         std::cout << "\tpindex " << pindex << std::endl;
//         std::cout << "\tijk " << i << " " << j << " " << k << std::endl;
//         std::cout << "\tNumber of cells " << cells.size() << std::endl;
//         std::cout << "\tNumber of neighbors " << neighbors.size() << std::endl;
//         std::cout << "\tMass " << mass(p) << std::endl;
//         std::cout << "\tDensity " << density << std::endl;
//      }
   } 
   //densiter++;
  // file.close();


}
//DFSPH 2015
//Equation above (12)
//ρ∗i = ρi + ∆t Dρi/Dt = ρi + ∆t ∑j mj (v∗i − v∗j )∇Wi
void SPHStateData::compute_predicted_density(size_t p, const double dt)
{

    float density = get_float_attr("density", p);
    float pdensity = 0.0f;
    const Vector P = pos(p);
    const Vector V = vel(p);
    std::vector<size_t> neighbors;
    neighbors_list(neighbors, P, neighbor_parallel);

    for(size_t a = 0; a < neighbors.size(); a++)
    {
        size_t pid = neighbors[a]; 

         const Vector V_j = vel(pid);
         pdensity += ((V-V_j) * grad_weight(pid,P));

    }

    pdensity *= get_float_attr("volume", p) * density0;
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
void SPHStateData::compute_density_derivative(size_t p)
{


    float density_change = 0.0;
    const Vector P = pos(p);
    const Vector V = vel(p);
    std::vector<size_t> neighbors;
    neighbors_list(neighbors, P, neighbor_parallel);

    for(size_t a = 0; a < neighbors.size(); a++)
    {
        size_t pid = neighbors[a]; 

        const Vector V_j = vel(pid);
        density_change += ((V-V_j) * grad_weight(pid,P));

    }

    density_change *= mass(p) ; //all fluid have constant volume


    // density_change = std::max(density_change, 0.0f);//mannually do this inside solver
    //#pragma omp critical
    //std::cout << "Dens deriv: " << density_change << '\n';
    if(std::isnan(density_change) || std::isinf(density_change)) std::cout << "densitychange bad\n";
    set_attr("density_derivative", p,  density_change);
   
 
}
//DFSPH 2015
//Equation (9)
//αi = ρi / |∑j mj ∇Wi|^2 + ∑j |mj ∇Wi|^2
//since alpha is used in ∇p, ρi is cancalled out, i.e. 1 / |∑j mj ∇Wi|^2 + ∑j |mj ∇Wi|^2
void SPHStateData::compute_factor()
{
    #pragma omp parallel for   
    for( size_t p =0; p < nb(); p++ )
    {
        float factor = 0;
        float sum_grad_p = 0.0f; // sum of grad pi and pj

        Vector grad_p_i(0,0,0); //pressure gradient of ith particle
        const Vector P = pos(p);
        std::vector<size_t> neighbors;
        neighbors_list(neighbors, P, neighbor_parallel);

        for(size_t a = 0; a < neighbors.size(); a++)
        {         
            size_t pid = neighbors[a]; 
       
            const Vector grad_p_j = get_float_attr("volume",pid) * grad_weight(pid, P) * density0;

            //dot product of a vector by itself is its magnitude squared
            sum_grad_p += grad_p_j * grad_p_j; 
            grad_p_i += grad_p_j;
        }
    

        sum_grad_p += grad_p_i * grad_p_i;


        if(sum_grad_p > m_eps)
        {
            //factor  =  get_float_attr("density",p) / sum_grad_p;
            factor  =  1.0 / sum_grad_p;
        }
        else
            factor = 0.f;
        if(std::isnan(factor) || std::isinf(factor)) std::cout << "factor bad\n";
        set_attr("factor", p, factor);
    }
}

void SPHStateData::populate()
{
   NeighborSearch::populate(*this);
}

void SPHStateData::set_radius( const float& v )
{
   radius = v;
   particle_radius = radius/4.0f;
   float particle_diamter = particle_radius * 2;
   float mV = 0.8f * particle_diamter * particle_diamter * particle_diamter;
   #pragma omp parallel for
   for(size_t p = 0; p < nb(); p++)
   {
      set_attr("volume", p, mV);

   }
   set_cellsize(2.0*radius); 
}

void SPHStateData::set_maxIter(int maxi)
{
   if(maxi < 2) maxi = 2;
   maxIter = maxi;
}

void SPHStateData::set_useUserDT(bool uudt)
{
   useUserDT = uudt;
}

void SPHStateData::set_mMaxError(int mmx)
{
   if(mmx < 0.01f ) mmx = 0.01f;
   m_maxError = mmx;
}

void SPHStateData::set_maxError(int mx)
{
   if(mx < 0.1f) mx = 0.1f;
   maxError = mx;
}

void SPHStateData::set_density0(float d0)
{
   if(d0 == 0) d0 =1;
   density0 = d0;
}
void SPHStateData::set_ddClamp(bool cl)
{
   dd_clamp = cl;
}
void SPHStateData::set_neighborParallel(bool cl)
{
   neighbor_parallel = cl;
}

float SPHStateData::max_velocity() const
{
    if(nb_items == 0) return 0;
    float vel_m = vel(0).magnitude();
    #pragma omp parallel for reduction(max:vel_m)
    for (size_t i = 1; i < nb(); i++) {
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