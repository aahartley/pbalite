//*******************************************************************
//
//   BVH.C
//
//  Bounding Volume Hierarchy for continous collison detection
//
//
//
//*******************************************************************

#include "BVH.h"

#include <iostream>

namespace pba
{

BVH::BVH(const Vector& llc, const Vector& urc, int level_, int max_levels_, int min_objects_) 
  : aabb_        (AABB(llc, urc)),
    node1_       (nullptr),
    node2_       (nullptr),
    level_       (level_),
    max_levels_  (max_levels_),
    min_objects_ (min_objects_) {}


BVH::~BVH()
{
   if (node1_) { delete node1_; }
   if (node2_) { delete node2_; }
}

void BVH::AddObject(CollisionSurface& o) 
{ 
    for (size_t i = 0; i < o.triangle_size(); ++i)
    {
        CollisionTriangle& tri = o.get_triangle(i);
        object_list_.push_back(&tri);
    } 
}

// TODO: add parallel ability with init sized object list
void BVH::AddObject(CollisionTriangle& o) 
{ 
    object_list_.push_back(&o);
}

void BVH::AddObject(CollisionTriangle* o) 
{ 
    object_list_.push_back(o);
}

bool BVH::Hit(const Ray& ray, const Vector& x_0, const Vector& x_u, const Vector& v, const double dt, CollisionData& c_data) const
{
    c_data.hit_tri = false;
    if (aabb_.RayIntersects(ray, 0, ray.length))
    {
        if (node1_ == nullptr && node2_ == nullptr)
        {
            //init leaf collision data
            c_data.hit_time = 2.0 * dt;
            c_data.hit_tri = false;
            //get first collision
            for (size_t i = 0; i < object_list_.size(); ++i)
            {
                Vector x_h_candidate;
                double dt_h_candidate = c_data.hit_time;
                if (object_list_[i]->Hit(x_0, x_u, v, dt, x_h_candidate, dt_h_candidate, 0)) //no radius for now
                {
                    if (std::fabs(dt_h_candidate) < std::fabs(c_data.hit_time))
                    {
                        c_data.hit_time = dt_h_candidate;
                        c_data.x_h = x_h_candidate;
	                    c_data.tri = object_list_[i];
                        if (!c_data.tri) std::cout<<"errorrrbvh\n";
                        c_data.hit_tri = true;
                    }
                }
            }  
            return c_data.hit_tri;
        }
        else
        {
            //traverse both nodes
            if (node1_ && node2_)
            {
                //each leaf gets own collision data
                CollisionData hit1;
                CollisionData hit2;
                bool result1 = node1_->Hit(ray, x_0, x_u, v, dt, hit1);
                bool result2 = node2_->Hit(ray, x_0, x_u, v, dt, hit2);
                if (!result1 && !result2) { return false; }
                if (result1 && !result2) { c_data = hit1;  return c_data.hit_tri; }
                if (!result1 && result2) { c_data =  hit2; return c_data.hit_tri; }
                //determine which leaf has first collision
                if (std::fabs(hit1.hit_time) > std::fabs(hit2.hit_time))
                {
                    c_data = hit2;
                    return c_data.hit_tri;
                }
                else
                {
                    c_data = hit1; 
                    return c_data.hit_tri;
                }
            }
            //triangles only in node2_
            else if(node1_ == nullptr)
            {
                return node2_->Hit(ray, x_0, x_u, v, dt, c_data );
            }
            //triangles only in node1_
            else if(node2_ == nullptr)
            {
                return node1_->Hit(ray, x_0, x_u, v, dt, c_data );
            }
        }
   }
   return c_data.hit_tri;
}


void BVH::Divide()
{
    if (level_ >= max_levels_ || object_list_.size() <= static_cast<size_t>(min_objects_)) return;
    AABB aabb1, aabb2;
    aabb_.Split(level_%3, aabb1, aabb2);
    if (node1_) delete node1_;
    node1_ = new BVH(aabb1.llc(), aabb1.urc(), level_+1, max_levels_, min_objects_);
    for (size_t i = 0; i < object_list_.size(); ++i)
    {
        if(aabb1.Intersects(object_list_[i]->aabb()))
        {
            node1_->AddObject(object_list_[i]);
        }
    }
    if (node1_->nbObjects() == 0) { delete node1_; node1_ = nullptr; }

    if (node2_) delete node2_;
    node2_ = new BVH(aabb2.llc(), aabb2.urc(), level_+1, max_levels_, min_objects_);
    for (size_t i = 0; i < object_list_.size(); ++i)
    {
        if (aabb2.Intersects(object_list_[i]->aabb()))
        {
            node2_->AddObject(object_list_[i]);
        }
    }
    if (node2_->nbObjects() == 0) { delete node2_; node2_ = nullptr; }

    size_t totalobjects = 0;
    if (node1_)
    { 
        totalobjects += node1_->nbObjects(); 
    }
    if (node2_)
    { 
        totalobjects += node2_->nbObjects(); 
    }
    if (totalobjects < nbObjects())
    {
        std::cout << "BVH Divide mistake at level_ " << level_ << ":  nb to be divided=" << nbObjects() << "   actual number " << totalobjects << std::endl;
    }
    //trianngles stored in leaf(s)
    object_list_.clear();
    if (node1_) node1_->Divide();
    if (node2_) node2_->Divide();
}


}//end of namespace pba