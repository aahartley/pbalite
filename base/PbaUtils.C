#include "PbaUtils.h"
#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.

using namespace pba;

void pba::CreateInfiniteCube(pba::CollisionSurface& s, double x, double y, double z )
{
    CollisionInfinitePlane bottom(Vector(0,1,0),Vector(0,-y,0));
    CollisionInfinitePlane top(Vector(0,-1,0),Vector(0,y,0));
    CollisionInfinitePlane right(Vector(-1,0,0),Vector(x,0,0));
    CollisionInfinitePlane left(Vector(1,0,0),Vector(-x,0,0));
    CollisionInfinitePlane front(Vector(0,0,-1),Vector(0,0,z)); //closest to screen
    CollisionInfinitePlane back(Vector(0,0,1),Vector(0,0,-z)); //(z points to screen)
    s->addPlane(bottom);
    s->addPlane(top);
    s->addPlane(right);
    s->addPlane(left);
    s->addPlane(back);
    s->addPlane(front);
}
//clockwise inward normals
void pba::CreateCube(pba::CollisionSurface& s, double x, double y, double z)
{
    Vector p000(-x, -y, -z);
    Vector p001(-x, -y,  z);
    Vector p010(-x,  y, -z);
    Vector p011(-x,  y,  z);
    Vector p100( x, -y, -z);
    Vector p101( x, -y,  z);
    Vector p110( x,  y, -z);
    Vector p111( x,  y,  z);

    // Bottom face (+Y normal)
    pba::CollisionTriangle b1 = makeCollisionTriangle(p100, p000, p101);
    std::cout << "b1: " << b1->getNormal().X() << ", " << b1->getNormal().Y() << ", " << b1->getNormal().Z() << "\n";
    s->addTriangle(b1);

    pba::CollisionTriangle b2 = makeCollisionTriangle(p101, p000, p001);
    std::cout << "b2: " << b2->getNormal().X() << ", " << b2->getNormal().Y() << ", " << b2->getNormal().Z() << "\n";
    s->addTriangle(b2);

    // Top face (-Y normal)
    pba::CollisionTriangle t1 = makeCollisionTriangle(p111, p011, p010);
    std::cout << "t1: " << t1->getNormal().X() << ", " << t1->getNormal().Y() << ", " << t1->getNormal().Z() << "\n";
    s->addTriangle(t1);

    pba::CollisionTriangle t2 = makeCollisionTriangle(p110, p111, p010);
    std::cout << "t2: " << t2->getNormal().X() << ", " << t2->getNormal().Y() << ", " << t2->getNormal().Z() << "\n";
    s->addTriangle(t2);

    // Left face (+X normal)
    pba::CollisionTriangle l1 = makeCollisionTriangle(p000, p010, p011);
    std::cout << "l1: " << l1->getNormal().X() << ", " << l1->getNormal().Y() << ", " << l1->getNormal().Z() << "\n";
    s->addTriangle(l1);

    pba::CollisionTriangle l2 = makeCollisionTriangle(p000, p011, p001);
    std::cout << "l2: " << l2->getNormal().X() << ", " << l2->getNormal().Y() << ", " << l2->getNormal().Z() << "\n";
    s->addTriangle(l2);

    // Right face (-X normal)
    pba::CollisionTriangle r1 = makeCollisionTriangle(p100, p101, p111);
    std::cout << "r1: " << r1->getNormal().X() << ", " << r1->getNormal().Y() << ", " << r1->getNormal().Z() << "\n";
    s->addTriangle(r1);

    pba::CollisionTriangle r2 = makeCollisionTriangle(p100, p111, p110);
    std::cout << "r2: " << r2->getNormal().X() << ", " << r2->getNormal().Y() << ", " << r2->getNormal().Z() << "\n";
    s->addTriangle(r2);

    //  Front -Z)
    pba::CollisionTriangle f1 = makeCollisionTriangle(p010, p110, p100);
    std::cout << "f1: " << f1->getNormal().X() << ", " << f1->getNormal().Y() << ", " << f1->getNormal().Z() << "\n";
    s->addTriangle(f1);

    pba::CollisionTriangle f2 = makeCollisionTriangle(p010, p100, p000);
    std::cout << "f2: " << f2->getNormal().X() << ", " << f2->getNormal().Y() << ", " << f2->getNormal().Z() << "\n";
    s->addTriangle(f2);

    // Back +Z) 
    pba::CollisionTriangle ba1 = makeCollisionTriangle(p011, p101, p111);
    std::cout << "ba1: " << ba1->getNormal().X() << ", " << ba1->getNormal().Y() << ", " << ba1->getNormal().Z() << "\n";
    s->addTriangle(ba1);

    pba::CollisionTriangle ba2 = makeCollisionTriangle(p011, p001, p101);
    std::cout << "ba2: " << ba2->getNormal().X() << ", " << ba2->getNormal().Y() << ", " << ba2->getNormal().Z() << "\n";
    s->addTriangle(ba2);


}



void pba::DisplayInfinitePlanes( pba::CollisionSurface& s )
{
    if(!s->is_visible()) { return; }
    float scale = 0.004;
    float dims = 1024*scale;
    Vector v1(1,0,0);
    for(size_t p =0; p < s->plane_size(); p++)
    {
        CollisionInfinitePlane plane = s->get_plane(p);
        Vector normal = plane.getNormal();
        Vector point = plane.getP0();
        if(normal.X() == 1 || normal.X() == -1) v1.set(0,1,0);
        Vector perp1 = v1^normal;
        perp1.normalize();
        //std::cout << perp1.X() << ' '<< perp1.Y()<< ' '<< perp1.Z() << '\n';
        Vector perp2 = normal^perp1;
        perp2.normalize();
        //std::cout << perp2.X() << ' '<< perp2.Y()<< ' '<< perp2.Z() << '\n';
        Vector p1 = point + (perp1 * dims) + (perp2 * dims);  
        Vector p2 = point - (perp1 * dims) + (perp2 * dims);  
        Vector p3 = point - (perp1 * dims) - (perp2 * dims);  
        Vector p4 = point + (perp1 * dims) - (perp2 * dims);
        GLenum mode = s->use_wireframe() ? GL_LINE_LOOP : GL_QUADS;
        glBegin(mode);  
        glColor3f(0.7f, 0.7f, 0.7f);
        glVertex3f(p1.X(), p1.Y(), p1.Z());
        glVertex3f(p2.X(), p2.Y(), p2.Z());
        glVertex3f(p3.X(), p3.Y(), p3.Z());
        glVertex3f(p4.X(), p4.Y(), p4.Z());
        glEnd(); 
        v1.set(1,0,0);

    }
}

void pba::Display( pba::CollisionSurface& s )
{
    if(!s->is_visible()) { return; }
    for(size_t p =0; p < s->triangle_size(); p++)
    {
        const CollisionTriangle& tri = s->get_triangle(p);
        Vector normal = tri->getNormal();
        Vector v0 = tri->getV0();
        Vector v1 = tri->getV1();
        Vector v2 = tri->getV2();
        GLenum mode = s->use_wireframe() ? GL_LINE_LOOP : GL_TRIANGLES;
        glBegin(mode);  
        glColor3f(0.7f, 0.7f, 0.7f);
        glVertex3f(v0.X(), v0.Y(), v0.Z());
        glVertex3f(v1.X(), v1.Y(), v1.Z());
        glVertex3f(v2.X(), v2.Y(), v2.Z());
        glEnd(); 
    }
}

void pba::GeoToSoftBody(std::vector<Vector>& verts, std::vector<Triangle>& tris, const std::string& file, SoftBodyState& s )
{
    ObjLoader::loadObj(verts, tris, file);
    s->add(verts.size());
    s->set_num_pairs((verts.size()*(verts.size()-1))/2);
    std::cout << "uniq edges "<<(verts.size()*(verts.size()-1))/2 << '\n';
    #pragma omp parallel for
    for(size_t i = 0; i < verts.size(); i++)
    {
        s->set_id(i, i);
        s->set_pos(i,Vector(verts[i].X(), verts[i].Y(), verts[i].Z()));
        s->set_mass(i,1);
        s->set_vel(i, Vector(0,0,0));
        s->set_ci(i, Color(0,0,1,1));;
    }
    #pragma omp parallel for schedule(dynamic)
    for(size_t i = 0; i < verts.size(); i++)
    {
        size_t start_count = (i * (2 * verts.size() - i - 1)) / 2;

        for(size_t j = i+1; j < verts.size(); j++)
        {
            size_t count = start_count + (j - i - 1);
            s->add_pair(i, j, count);
        }
    }
}
