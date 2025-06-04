//*******************************************************************
//
//   CollisionHandler.C
//
//   Continous collison detection, resolve all hits during one time step
//
//
//
//*******************************************************************

#include "CollisionHandler.h"

#include <cmath>
#include <iostream>

namespace pba
{

void CollisionHandlerBase::set_collision_surface(CollisionSurface* c)
{
    surf_ = c;
}

void CollisionHandlerBase::set_bvh(const BVH* b)
{
    bvh_ = b;
}

ElasticCollisionHandler::ElasticCollisionHandler(){}
ElasticCollisionHandler::~ElasticCollisionHandler(){}
// Check if positions imply a collision has already taken place within the allotted time.
// If so, backs the position up along the velocity direction to the point of impact, then
// does an elastic bounce, repat until end of time step
void ElasticCollisionHandler::HandleCollisions(const double dt, DynamicalStateData& pq) 
{
    if (surf_ == nullptr) return;
    #pragma omp parallel for
    for (size_t i = 0; i < pq.nb(); ++i)
    {
        Vector v_0 = pq.vel(i); 
        //get start pos before integration (pos 1st, pos 1st (-dt), vel 1st, vel 1st (-dt))
        Vector x_0 = (pq.pos(i) - v_0*dt);
        //coll results
        Vector x_r = pq.pos(i); //end pos (current pos)
        Vector v_r = v_0;
        float radius = pq.rad(i);
        CollisionData data;
        double running_dt = dt;
        bool keep_checking_for_hits = true;
        int count = 0;
        while (keep_checking_for_hits)
        {
            keep_checking_for_hits = false;
            //check for coll
            if (use_bvh_ && bvh_)
            {
                Vector direction = (x_r - x_0).unitvector();
                Ray ray(x_0, x_r, direction);
                if (bvh_->Hit(ray, x_0, x_r, v_0, running_dt, data))
                {
                    //std::cout << "hit worked\n";
                    keep_checking_for_hits = true;
                    if (data.hit_tri)
                    {
                        data.tri->Handle(v_0, running_dt, data.x_h, data.hit_time, x_r, 
                                         v_r, surf_->coeff_sticky(), surf_->coeff_restitution());
                    }
                    x_0 = data.x_h;
                    v_0 = v_r;
                    running_dt = running_dt - data.hit_time;
        
                    if ((dt >= 0 && running_dt <= 0) || (dt < 0 && running_dt >= 0))
                    {
                        keep_checking_for_hits = false;
                    }              
                }
            }
            else
            {
                if (surf_->Hit(x_0, x_r, v_0, running_dt, data , radius))
                {
                    //std::cout << "hit worked\n";
                    keep_checking_for_hits = true;
                    // Handle collision on plane, with the smallest dtH at hit point XH
                    //std::cout << data.hit_time << '\n';
                    if (data.hit_plane)
                    {
                        surf_->get_plane(data.hit_index).Handle(v_0, running_dt, data.x_h, data.hit_time,
                                                                    x_r, v_r, surf_->coeff_sticky(), surf_->coeff_restitution());
                        // std::cout << "hit_tri:(handle) " << data.hit_tri << '\n';
                        // std::cout << "hit_index:(handle) " << data.hit_index << '\n';
                    }             
                    if (data.hit_tri)
                    {
                        surf_->get_triangle(data.hit_index).Handle(v_0, running_dt, data.x_h, data.hit_time,
                                                                    x_r, v_r, surf_->coeff_sticky(), surf_->coeff_restitution());
                        // std::cout << "hit_tri:(handle) " << data.hit_tri << '\n';
                        // std::cout << "hit_index:(handle) " << data.hit_index << '\n';
                    }
                    // if(data.hit_time==0)
                    // {   
                    //     //std::cout << "hittime is 0\n";
                    //     keep_checking_for_hits=false;
                    // }
                    x_0 = data.x_h;
                    v_0 = v_r;
                    running_dt = running_dt - data.hit_time;
                    count++;
                    if (count == 10) 
                    {
                        //keep_checking_for_hits = false;
                        std::cout << "big hits\n";
                    }
                    //std::cout << "hit " << count << " : "<< running_dt << '\n';
                    if ((dt >= 0 && running_dt <= 0) || (dt < 0 && running_dt >= 0))
                    {
                        keep_checking_for_hits = false;
                    }              
                }
            }
        }
        pq.set_pos(i, x_r);
        pq.set_vel(i, v_r);
    }
}

}//end of pba namespace