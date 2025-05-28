#include "BVH.h"

using namespace pba;

BVHRaw::BVHRaw( const Vector& llc, const Vector& urc, const int lvl, const int maxlvl, const int minobj ) :
   aabb        (AABB(llc,urc)),
   node1       (0),
   node2       (0),
   level       (lvl),
   max_levels  (maxlvl),
   min_objects (minobj)
{}


BVHRaw::~BVHRaw()
{
   if( node1 ){ delete node1; }
   if( node2 ){ delete node2; }
}

void BVHRaw::addObject(const CollisionSurface& o ) 
{ 
    for(size_t i = 0; i < o->triangle_size(); i++)
    {
        const CollisionTriangle& tri = o->get_triangle(i);
        if(!tri)std::cout << "Error\n";
        object_list.push_back(tri);
    } 
}

void BVHRaw::addObject(const CollisionTriangle& o ) 
{ 
    if(!o)std::cout << "Error\n";
    object_list.push_back(o);
}

bool BVHRaw::hit(const Ray& ray, const Vector& X0, const Vector& Xu, const Vector& V, const double dt, CollisionData& c_data ) const
{
    c_data.hit_tri = false;
    if(aabb.ray_intersect(ray, 0, ray.length))
    {
        if(node1 == 0 && node2 == 0)
        {
            c_data.hit_time = 2.0 * dt;
            c_data.hit_tri = false;
            //get first collision
            for( size_t i = 0; i < object_list.size(); i++ )
            {
                Vector XH_candidate;
                double dtH_candidate = c_data.hit_time;
                if( object_list[i]->hit( X0, Xu, V, dt, XH_candidate, dtH_candidate, 0 )   )
                {
                    if( std::fabs(dtH_candidate) < std::fabs(c_data.hit_time) )
                    {
                        c_data.hit_time = dtH_candidate;
                        c_data.XH = XH_candidate;
	                    c_data.tri = object_list[i];
                        if(!c_data.tri)std::cout<<"errorrrbvh\n";
                        c_data.hit_tri = true;
                    }
                }
            }  
            return c_data.hit_tri;
        }
        else
        {
            //traverse both nodes
            if( node1 != 0 && node2 != 0 )
            {
                CollisionData hit1;
                CollisionData hit2;
                bool result1 = node1->hit(ray, X0, Xu, V, dt, hit1 );
                bool result2 = node2->hit(ray, X0, Xu, V, dt, hit2 );
                if( !result1 && !result2 ){ return false; }
                if( result1 && !result2 ){ c_data = hit1;  return c_data.hit_tri; }
                if( !result1  && result2 ){ c_data =  hit2; return c_data.hit_tri; }
                //determine which leaf has first collision
                if( std::fabs(hit1.hit_time) > std::fabs(hit2.hit_time) )
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
            //trinagles only in node2
            else if( node1 == 0 )
            {
                return node2->hit(ray, X0, Xu, V, dt, c_data );
            }
            //trinagles only in node1
            else if( node2 == 0 )
            {
                return node1->hit(ray, X0, Xu, V, dt, c_data );
            }
        }
   }
   return c_data.hit_tri;
}


void BVHRaw::Divide()
{
    if(level >= max_levels || object_list.size() <= (size_t)min_objects) return;
    AABB aabb1, aabb2;
    aabb.split(level%3, aabb1, aabb2);
    if(node1 != 0) delete node1;
    node1 = new BVHRaw(aabb1.getLLC(), aabb1.getURC(), level+1, max_levels, min_objects);
    for( size_t i =0 ; i < object_list.size(); i++ )
    {
        if( aabb1.intersects( object_list[i]->getAABB() ) )
        {
            node1->addObject( object_list[i] );
        }
    }
    if( node1->nbObjects() == 0 ){ delete node1; node1 = 0; }

    if(node2 != 0) delete node2;
    node2 = new BVHRaw(aabb2.getLLC(), aabb2.getURC(), level+1, max_levels, min_objects);
    for( size_t i =0 ; i < object_list.size(); i++ )
    {
        if( aabb2.intersects( object_list[i]->getAABB() ) )
        {
            node2->addObject( object_list[i] );
        }
    }
    if( node2->nbObjects() == 0 ){ delete node2; node2 = 0; }
    size_t totalobjects = 0;
    if( node1 )
    { 
        totalobjects += node1->nbObjects(); 
    }
    if( node2 )
    { 
        totalobjects += node2->nbObjects(); 
    }

    if( totalobjects < nbObjects() )
    {
        std::cout << "BVH Divide mistake at level " << level << ":  nb to be divided=" << nbObjects() << "   actual number " << totalobjects << std::endl;
    }
    //trianngles stored in node(s)
    object_list.clear();
    if(node1) node1->Divide();
    if(node2) node2->Divide();

    
}

BVH pba::makeBVH(const Vector& llc, const Vector& urc, const int lvl, const int maxlvl, const int minobj)
{
    return std::make_shared<BVHRaw>(llc, urc, lvl, maxlvl, minobj);
}
