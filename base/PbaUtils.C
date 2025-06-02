#include "PbaUtils.h"

#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.

#include <iostream>

namespace pba
{

void CreateInfiniteCube(CollisionSurface& s, double x, double y, double z)
{
    CollisionInfinitePlane bottom(Vector(0,1,0),Vector(0,-y,0));
    CollisionInfinitePlane top(Vector(0,-1,0),Vector(0,y,0));
    CollisionInfinitePlane right(Vector(-1,0,0),Vector(x,0,0));
    CollisionInfinitePlane left(Vector(1,0,0),Vector(-x,0,0));
    CollisionInfinitePlane front(Vector(0,0,-1),Vector(0,0,z)); //closest to screen
    CollisionInfinitePlane back(Vector(0,0,1),Vector(0,0,-z)); //(z points to screen)
    s.AddPlane(bottom);
    s.AddPlane(top);
    s.AddPlane(right);
    s.AddPlane(left);
    s.AddPlane(back);
    s.AddPlane(front);
}
//clockwise inward normals
void CreateCube(CollisionSurface& s, double x, double y, double z)
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
    pba::CollisionTriangle b1(p100, p000, p101);
    std::cout << "b1: " << b1.normal().X() << ", " << b1.normal().Y() << ", " << b1.normal().Z() << "\n";
    s.AddTriangle(b1);

    pba::CollisionTriangle b2(p101, p000, p001);
    std::cout << "b2: " << b2.normal().X() << ", " << b2.normal().Y() << ", " << b2.normal().Z() << "\n";
    s.AddTriangle(b2);

    // Top face (-Y normal)
    pba::CollisionTriangle t1(p111, p011, p010);
    std::cout << "t1: " << t1.normal().X() << ", " << t1.normal().Y() << ", " << t1.normal().Z() << "\n";
    s.AddTriangle(t1);

    pba::CollisionTriangle t2(p110, p111, p010);
    std::cout << "t2: " << t2.normal().X() << ", " << t2.normal().Y() << ", " << t2.normal().Z() << "\n";
    s.AddTriangle(t2);

    // Left face (+X normal)
    pba::CollisionTriangle l1(p000, p010, p011);
    std::cout << "l1: " << l1.normal().X() << ", " << l1.normal().Y() << ", " << l1.normal().Z() << "\n";
    s.AddTriangle(l1);

    pba::CollisionTriangle l2(p000, p011, p001);
    std::cout << "l2: " << l2.normal().X() << ", " << l2.normal().Y() << ", " << l2.normal().Z() << "\n";
    s.AddTriangle(l2);

    // Right face (-X normal)
    pba::CollisionTriangle r1(p100, p101, p111);
    std::cout << "r1: " << r1.normal().X() << ", " << r1.normal().Y() << ", " << r1.normal().Z() << "\n";
    s.AddTriangle(r1);

    pba::CollisionTriangle r2(p100, p111, p110);
    std::cout << "r2: " << r2.normal().X() << ", " << r2.normal().Y() << ", " << r2.normal().Z() << "\n";
    s.AddTriangle(r2);

    //  Front -Z)
    pba::CollisionTriangle f1(p010, p110, p100);
    std::cout << "f1: " << f1.normal().X() << ", " << f1.normal().Y() << ", " << f1.normal().Z() << "\n";
    s.AddTriangle(f1);

    pba::CollisionTriangle f2(p010, p100, p000);
    std::cout << "f2: " << f2.normal().X() << ", " << f2.normal().Y() << ", " << f2.normal().Z() << "\n";
    s.AddTriangle(f2);

    // Back +Z) 
    pba::CollisionTriangle ba1(p011, p101, p111);
    std::cout << "ba1: " << ba1.normal().X() << ", " << ba1.normal().Y() << ", " << ba1.normal().Z() << "\n";
    s.AddTriangle(ba1);

    pba::CollisionTriangle ba2(p011, p001, p101);
    std::cout << "ba2: " << ba2.normal().X() << ", " << ba2.normal().Y() << ", " << ba2.normal().Z() << "\n";
    s.AddTriangle(ba2);


}



void DisplayInfinitePlanes(CollisionSurface* s)
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

void Display(CollisionSurface* s)
{
    if(!s->is_visible()) { return; }
    for(size_t p =0; p < s->triangle_size(); p++)
    {
        const CollisionTriangle& tri = s->get_triangle(p);
        //const Vector& normal = tri.normal();
        const Vector& v0 = tri.v0();
        const Vector& v1 = tri.v1();
        const Vector& v2 = tri.v2();
        GLenum mode = s->use_wireframe() ? GL_LINE_LOOP : GL_TRIANGLES;
        glBegin(mode);  
        glColor3f(0.7f, 0.7f, 0.7f);
        glVertex3f(v0.X(), v0.Y(), v0.Z());
        glVertex3f(v1.X(), v1.Y(), v1.Z());
        glVertex3f(v2.X(), v2.Y(), v2.Z());
        glEnd(); 
    }
}

void GeoToSoftBody(std::vector<Vector>& verts, std::vector<Triangle>& tris, const std::string& file, SoftBodyStateData& s )
{
    ObjLoader::loadObj(verts, tris, file);
    s.Add(verts.size());
    s.set_num_pairs((verts.size()*(verts.size()-1))/2);
    std::cout << "uniq edges "<<(verts.size()*(verts.size()-1))/2 << '\n';
    #pragma omp parallel for
    for(size_t i = 0; i < verts.size(); i++)
    {
        s.set_id(i, i);
        s.set_pos(i,Vector(verts[i].X(), verts[i].Y(), verts[i].Z()));
        s.set_mass(i,1);
        s.set_vel(i, Vector(0,0,0));
        s.set_ci(i, Color(0,0,1,1));;
    }
    #pragma omp parallel for schedule(dynamic)
    for(size_t i = 0; i < verts.size(); i++)
    {
        size_t start_count = (i * (2 * verts.size() - i - 1)) / 2;

        for(size_t j = i+1; j < verts.size(); j++)
        {
            size_t count = start_count + (j - i - 1);
            s.AddPair(i, j, count);
        }
    }
}

}