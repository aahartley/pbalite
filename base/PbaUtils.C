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