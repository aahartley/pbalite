#include "SoftBodyThing.h"

#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.

#include <iostream>
#include <cstdlib>


namespace pba
{
//memeory flow: if remake state- must reset solvers
SoftBodyThing::SoftBodyThing(const std::string nam) 
             : PbaThingyDingy (nam),
               emit           (false)
{
    box = CreateCollisionSurface();
    bvh = CreateBVH(Vector(-3,-3,-3), Vector(3,3,3), 0, 20, 1);
    //pba::CreateInfiniteCube(box, 2, 1, 2);
    pba::CreateCube(*box, 2, 0.5, 2);
    bvh->AddObject(*box);
    bvh->Divide();
    //if deleting surface, set collisons surf and bvh * to null, reset bvh 
    AddCollisionSurface(*box);
    state = CreateSoftBodyState("SoftBody");
    Reset();
    // int inc = 2;
    // state->add(inc);
    // state->set_num_pairs(1);
    std::cout << "Emit: Total Points " << state->nb() << std::endl;
    emitter = ParticleEmitter();
  
    // Vector p(0,2.9,0);
    // Vector v(0,0,0);
    // Color c(0,0,1,1);
    // state->set_pos(0,p);
    // state->set_vel(0,v);
    // state->set_ci(0,c);
    // state->set_rad(0, 0.5f); //0.075

    // state->set_pos(1,Vector(0,2.7,0));
    // state->set_vel(1,v);
    // state->set_ci(1,c);
    // state->set_rad(1, 0.5f); //0.075
    // state->add_pair(0,1,0);
    verts.resize(1); tris.resize(1);
    GeoToSoftBody(verts[0], tris[0], "models/bunny_superlo_scaled.obj", *state);
    //GeoToSoftBody("models/bunny_lo_scaled.obj", state);

    std::cout << state->nb_pairs() << '\n';
    force = CreateAccumulatingForce<SoftBodyStateData>();

    gravityforce = CreateGravityForce<SoftBodyStateData>(Vector(0,-9.81f,0));
    struts = CreateAccumulatingStrutForce(20, 2, false); //10, 0.5
    //if removing a Force, reset accumulatingforce first
    AccumulatingForce<SoftBodyStateData>* f = dynamic_cast<AccumulatingForce<SoftBodyStateData>*>(force.get()); 
	f->AddForce(gravityforce.get());
    f->AddForce(struts.get());
    //if resetting state/force/collisions, reset solvers first
    //if resetting solvers setup, reset solver first
    a = CreateAdvancePosition(*state, &collisions);
    b = CreateAdvanceVelocity<SoftBodyStateData>(*state, *force, 1.0e9, 1.0e9);
    GISolverPtr basicsolver = CreateLeapFrogSolver(a.get(), b.get());
    basicsolver = CreateGISolverSixthOrder(std::move(basicsolver));
    solver = CreateGISolverSubstep(std::move(basicsolver), 2);
    std::cout << name << " constructed\n";

}

SoftBodyThing::~SoftBodyThing(){}

void SoftBodyThing::Init(const std::vector<std::string>& args) 
{
    SetSimulationTimestep(0.005);
}
    
void SoftBodyThing::Display() 
{
    pba::Display(box.get());
    // glPointSize(5.0);
    // glBegin(GL_POINTS);
    // for( size_t i=0;i<state->nb();i++ )
    // {
    //     const Vector& P = state->pos(i);
    //     const Color& ci = state->ci(i);
    //     glColor3f( ci.red(), ci.green(), ci.blue() );
    //     glVertex3f( P.X(), P.Y(), P.Z() );
    // }
    // glEnd();
    glBegin(GL_TRIANGLES);
    for (size_t j = 0; j < tris.size(); ++j)
    {
        for (size_t i = 0; i < tris[j].size(); ++i)
        {
            const Vector& v1 = state->pos(tris[j][i].v1);
            const Vector& v2 = state->pos(tris[j][i].v2);
            const Vector& v3 = state->pos(tris[j][i].v3);

            const Color& ci = Color(0,0,1,1);
            glColor3f( ci.red(), ci.green(), ci.blue() );
            glVertex3f( v1.X(), v1.Y(), v1.Z() );
            glVertex3f( v2.X(), v2.Y(), v2.Z() );
            glVertex3f( v3.X(), v3.Y(), v3.Z() );

        }
    }
    glEnd();
}

void SoftBodyThing::Keyboard(unsigned char key, int x, int y)
{
    PbaThingyDingy::Keyboard(key,x,y);
    if (key == 'v') { box->toggle_visible(); }
    if (key == 'w') { box->toggle_wireframe(); }
    if (key == 'e') { emit = !emit; }
    if (key == 'z')
    {
        auto* f = dynamic_cast<GravityForce<SoftBodyStateData>*>(gravityforce.get()); 
        Vector wind = f->get_gravity() + Vector(2,0,0);
        f->set_gravity(wind );
        std::cout << "gravity: " << f->get_gravity().X() << f->get_gravity().Y() << std::endl;

    }
    if (key == 'g')
    {
        auto* f = dynamic_cast<GravityForce<SoftBodyStateData>*>(gravityforce.get()); 
        f->set_gravity(f->get_gravity()/1.1);
        std::cout << "gravity: " << f->get_gravity().Y() << std::endl;
    }
    if (key == 'G')
    { 
        auto* f = dynamic_cast<GravityForce<SoftBodyStateData>*>(gravityforce.get()); 
        f->set_gravity( f->get_gravity()*1.1 );
        std::cout << "gravity: " << f->get_gravity().Y() << std::endl;
    }
    if (key == 'c')
    {
        box->set_coeff_restitution( box->coeff_restitution()/1.1 );
        std::cout << "coefficient of restituion: " << box->coeff_restitution() << std::endl;
    }
    if (key == 'C')
    { 
        box->set_coeff_restitution( box->coeff_restitution()*1.1 );
        std::cout << "coefficient of restituion: " << box->coeff_restitution() << std::endl;
    }
    if (key == 's')
    {
        box->set_coeff_sticky( box->coeff_sticky()/1.1 );
        std::cout << "coefficient of sticky: " << box->coeff_sticky() << std::endl;
    }
    if (key == 'S')
    { 
        box->set_coeff_sticky( box->coeff_sticky()*1.1 );
        std::cout << "coefficient of sticky: " << box->coeff_sticky() << std::endl;
    }
    if (key == 'l')
    {
        solver.reset();
        a = CreateAdvancePosition(*state, &collisions);
        b = CreateAdvanceVelocity<SoftBodyStateData>(*state, *force, 1.0e9, 1.0e9);
        solver = CreateLeapFrogSolver(a.get(), b.get());
        std::cout << "Using Leap Frog solver" << std::endl;
    }
//     if( key == 'n' )
//     {
//         GISolver solvera = CreateAdvancePositionColl( state, collisions );
//         GISolver solverb = CreateAdvanceVelocity(state,force);
//         solver = CreateForwardEulerSolver(solvera,solverb); // forward
//    std::cout << "Using Forward Euler solver" << std::endl;
//     }
//     if( key == 'b' )
//     {
//         GISolver solverb = CreateAdvanceVelocity(state,force);
//         GISolver solvera = CreateAdvancePositionColl( state, collisions );
//         solver = CreateForwardEulerSolver(solverb,solvera); //backward
//         std::cout << "Using Backward Euler solver" << std::endl;
//     }

}


void SoftBodyThing::solve()
{
    if (emit)
    {
        //emitter.emitCube(state, 6, Vector(0,0,0));
    }
    solver->Solve(dt);
}

void SoftBodyThing::Reset()
{
    state->Clear();
}

void SoftBodyThing::Usage()
{
   PbaThingyDingy::Usage();
   cout << "=== " << name << " ===\n";
   cout << "v            toggle visibility of collision surface\n";
   cout << "w            toggle wireframe/solid display of collision surface\n";
   cout << "g/G          reduce/increase gravityy \n";
   cout << "e            toggle particle emission on/off\n";
   cout << "c/C          reduce/increase coefficient of restitution\n";
   cout << "s/S          reduce/increase coefficient of sticky\n";
   cout << "l            use leap frog solver\n";
   cout << "n            use forward euler solver\n";
   cout << "b            use backward euler solver\n";
}

void SoftBodyThing::AddCollisionSurface(CollisionSurface& s)
{
    std::cout << "Add CollisionSurface\n";
    if (box == nullptr)
    {
        std::cout << "adding new surface\n";
        box = std::make_unique<CollisionSurface>(s);
    }
    s.set_coeff_restitution(1);
    //s->set_coeff_sticky(0.1);
    //if deleting bvh/surf, set collision bvh* to nullptr, same for surface
    collisions.set_collision_surface(box.get());
    collisions.set_bvh(bvh.get());
}


pba::PbaThing CreateSoftBodyThing() { return PbaThing(new SoftBodyThing()); }

}//end of pba namespace


