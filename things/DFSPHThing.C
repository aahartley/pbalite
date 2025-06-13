

#include "DFSPHThing.h"

#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.

#include <iostream>
#include <cstdlib>
#include <chrono>

namespace pba
{

DFSPHThing::DFSPHThing(const std::string nam) 
  : PbaThingyDingy(nam),
    emit(false),
    solve_time_(0)
{
    box = CreateCollisionSurface();
    float x = 0.3;
    float y = 0.7;
    float z = 0.3;
    pba::CreateCube(*box, x, y, z);
    bvh = CreateBVH(Vector(-5,-5,-5), Vector(5,5,5), 0, 20, 1);
    bvh->AddObject(*box);
    bvh->Divide();
    AddCollisionSurface(*box);
    state = CreateSPH(AABB(Vector(-5,-5,-5), Vector(5,5,5)), 0.1, "DFSPHState");
    Reset();
    // int inc = 1;
    // state->add(inc);

    // std::cout << "Emit: Total Points " << state->nb() << std::endl;
    // emitter = ParticleEmitter();
    // for(size_t i = 0; i < state->nb(); i++)
    // {
    //     Vector p(0,2.9,0);
    //     Vector v(0,0,0);
    //     Color c(0,0,1,1);
    //     state->set_pos(i,p);
    //     state->set_vel(i,v);
    //     state->set_ci(i,c);
    //     state->set_mass(i, state->get_float_attr("volume", i) * state->get_density0());
    //     state->set_id(i,i);

    // }

    force = CreateAccumulatingForce<SPHStateData>();

    gravityforce = CreateGravityForce<SPHStateData>(Vector(0,-9.81f,0));
    viscosity = CreateExplicitViscosity(0.01);

    AccumulatingForce<SPHStateData>* f = dynamic_cast<AccumulatingForce<SPHStateData>*>(force.get()); 
	f->AddForce(gravityforce.get());
    f->AddForce(viscosity.get());
    // GISolver a = CreateAdvancePositionCollSPH(state, collisions);
    // GISolver b = CreateAdvanceVelocitySPH(state, force, 1.0e9, 1.0e9);
    // GISolver b_euler = CreateBackwardEulerSolver(a, b);
    solver = CreateDFSPHSolver(*state, *force, &collisions, 1.0e9, 1.0e9);
    std::cout << name << " constructed\n";

}

DFSPHThing::~DFSPHThing(){}

void DFSPHThing::Init( const std::vector<std::string>& args ) 
{
    SetSimulationTimestep(0.01);
}
    
void DFSPHThing::Display() 
{
//     pba::Display(box.get());
//     // glPointSize(5.0);
//     // glBegin(GL_POINTS);
//     // for( size_t i=0;i<state->nb();i++ )
//     // {
//     //     const Vector& P = state->pos(i);
//     //     const Color& ci = state->ci(i);
//     //     glColor3f( ci.red(), ci.green(), ci.blue() );
//     //     glVertex3f( P.X(), P.Y(), P.Z() );
//     // }
//     // glEnd();
//     for (size_t i = 0; i < state->nb(); ++i)
//     {
//        const Color& ci = state->ci(i);
//        const pba::Vector& v = state->pos(i);
//        glPushMatrix();
//        glColor3f(ci.red(), ci.green(), ci.blue());
//        glTranslatef(v.X(), v.Y(),v.Z());
//        glutSolidSphere(0.025, 30,30);
//        glPopMatrix();
//     }
}

void DFSPHThing::Keyboard( unsigned char key, int x, int y )
{
    PbaThingyDingy::Keyboard(key,x,y);
    if (key == 'v') { box->toggle_visible(); }
    if (key == 'w') { box->toggle_wireframe(); }
    if (key == 'e') { emit = !emit; }
    if (key == 'z')
    {
        auto* f = dynamic_cast<GravityForce<SPHStateData>*>(gravityforce.get()); 
        Vector wind = f->get_gravity() + Vector(2,0,0);
        f->set_gravity(wind );
    }
    if (key == 'g')
    {
        auto* f = dynamic_cast<GravityForce<SPHStateData>*>(gravityforce.get()); 
        f->set_gravity(f->get_gravity()/1.1);
        
    }
    if(key == 'G')
    { 
        auto* f = dynamic_cast<GravityForce<SPHStateData>*>(gravityforce.get()); 
        f->set_gravity( f->get_gravity()*1.1 );
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
    if  (key == 'S')
    { 
        box->set_coeff_sticky( box->coeff_sticky()*1.1 );
        std::cout << "coefficient of sticky: " << box->coeff_sticky() << std::endl;
    }
//     if( key == 'l' )
//     {
//         GISolver solvera = CreateAdvancePositionColl( state, collisions );
//         GISolver solverb = CreateAdvanceVelcity(state,force);
//         solver = CreateLeapFrogSolver(solvera,solverb);
//         std::cout << "Using Leap Frog solver" << std::endl;
//     }
//     if( key == 'n' )
//     {
//         GISolver solvera = CreateAdvancePositionColl( state, collisions );
//         GISolver solverb = CreateAdvanceVelcity(state,force);
//         solver = CreateForwardEulerSolver(solvera,solverb); // forward
//    std::cout << "Using Forward Euler solver" << std::endl;
//     }
//     if( key == 'b' )
//     {
//         GISolver solverb = CreateAdvanceVelcity(state,force);
//         GISolver solvera = CreateAdvancePositionColl( state, collisions );
//         solver = CreateForwardEulerSolver(solverb,solvera); //backward
//         std::cout << "Using Backward Euler solver" << std::endl;
//     }

}


void DFSPHThing::solve()
{
    auto start = std::chrono::high_resolution_clock::now();
    if (emit)
    {
        emitter.EmitCube(*state, 6, Vector(0,0,0));
        emit = false;
    }
    solver->Solve(dt);
    state->EraseOutsideOfBounds(Vector(-5,-5,-5), Vector(5,5,5));
    ++frame;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    solve_time_ += static_cast<double>(duration.count());
    if (frame % 24 == 0)
    {
        std::cout << "Avg solve time: " << solve_time_/24 << " ms" << std::endl;
        solve_time_ = 0;
    }
}

void DFSPHThing::Reset()
{
    state->Clear();
}

void DFSPHThing::Usage()
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

void DFSPHThing::AddCollisionSurface(CollisionSurface& s)
{
    std::cout << "Add CollisionSurface\n";
    if(box == nullptr)
    {
        box = std::make_unique<CollisionSurface>(s);
    }
    s.set_coeff_restitution(0.5);
    //s->set_coeff_sticky(0.1);
    collisions.set_collision_surface(box.get());
    collisions.set_bvh(bvh.get());
}


pba::PbaThing CreateDFSPHThing() { return PbaThing(new DFSPHThing()); }


}
