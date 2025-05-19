#include "PbaUtils.h"
#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.

using namespace pba;

void pba::Display( pba::CollisionSurface& s )
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

void pba::GeoToSoftBody(const std::string& file, SoftBodyState& s )
{
    std::vector<Vector> verts = ObjLoader::loadObjVert(file);
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
