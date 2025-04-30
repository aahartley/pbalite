#include "CollisionHandler.h"
#include "LinearAlgebra.h"
#include <cmath>
#include <iostream>


using namespace pba;


void CollisionHandler::set_collision_surface(CollisionSurface& c)
{
    surf = c;
}

ElasticCollisionHandler::ElasticCollisionHandler(){}
ElasticCollisionHandler::~ElasticCollisionHandler(){}
// Check if positions imply a collision has already taken place within the allotted time.
// If so, backs the position up along the velocity direction to the point of impact, then
// does an elastic bounce, repat until end of time step
void ElasticCollisionHandler::handle_collisions(const double dt, DynamicalState& PQ) 
{
    #pragma omp parallel for
    for( size_t i=0;i<PQ->nb();i++ )
    {
        Vector V0 = PQ->vel(i); //curr vel
        Vector X0 = PQ->pos(i) - V0*dt; //get start pos before integration
        //coll results
        Vector XR = PQ->pos(i); //end pos
        Vector VR = V0;
        float radius = PQ->rad(i);
        CollisionData data;
        double running_dt = dt;
        bool keep_checking_for_hits = true;
        //int count = 0;
        while(keep_checking_for_hits)
        {
            keep_checking_for_hits = false;
            //check for coll
            if( surf->hit( X0, XR, V0, running_dt, data , radius))
            {
                keep_checking_for_hits = true;
                // Handle collision on plane, with the smallest dtH at hit point XH
                //std::cout << data.hit_time << '\n';
                if(data.hit_plane)
                    surf->get_plane(data.hit_index).handle( X0, V0, running_dt, data.XH, data.hit_time, XR, VR, surf->coeff_sticky(), surf->coeff_restitution());
                // if(data.hit_time==0)
                // {   
                //     //std::cout << "hittime is 0\n";
                //     keep_checking_for_hits=false;
                // }
                X0 = data.XH;
                V0 = VR;
                //std::cout << "vel: "<<V0.Y() << '\n';
                running_dt = running_dt - data.hit_time;
                //std::cout << "hit " << count << " : "<< running_dt << '\n';
                //count++;
                if( running_dt <= 0.0 ){ keep_checking_for_hits = false; }
                //break;
            }
        }
        PQ->set_pos( i, XR );
        PQ->set_vel( i, VR );
    }
}