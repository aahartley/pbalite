

#include "GravityThing.h"
#include <cstdlib>
#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.
#include <iostream>



using namespace std;

using namespace pba;





GravityThing::GravityThing(const std::string nam) :
 PbaThingyDingy (nam),
 emit       (false)
{
    box = CreateCollisionSurface();
    pba::CreateCube(*box, 2, 2, 2);
    AddCollisionSurface(*box);
    state = CreateDynamicalState("ParticleState");
    Reset();
    int inc = 1;
    state->Add(inc);

    std::cout << "Emit: Total Points " << state->nb() << std::endl;
    emitter = ParticleEmitter();
    for(size_t i = 0; i < state->nb(); i++)
    {
        Vector p(0,1,0);
        Vector v(0,0,0);
        Color c(0,0,1,1);
        state->set_pos(i,p);
        state->set_vel(i,v);
        state->set_ci(i,c);
        state->set_rad(i, 0.5f); //0.075


    }

    force = CreateAccumulatingForce<DynamicalStateData>();

    gravityforce = CreateGravityForce<DynamicalStateData>(Vector(0,-9.81f,0));

    AccumulatingForce<DynamicalStateData>* f = dynamic_cast<AccumulatingForce<DynamicalStateData>*>(force.get()); 
	f->AddForce(gravityforce.get());
    //GISolver a = CreateAdvancePosition(state);
    a = CreateAdvancePosition(*state, &collisions);
    b = CreateAdvanceVelocity(*state, *force, 1.0e9, 1.0e9);
    solver = CreateForwardEulerSolver(a.get(), b.get());
    std::cout << name << " constructed\n";

}

GravityThing::~GravityThing(){}

void GravityThing::Init( const std::vector<std::string>& args ) 
{
    SetSimulationTimestep(0.01);
}
    
void GravityThing::Display() 
{
    pba::Display(box.get());
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for( size_t i=0;i<state->nb();i++ )
    {
        const Vector& P = state->pos(i);
        const Color& ci = state->ci(i);
        glColor3f( ci.red(), ci.green(), ci.blue() );
        glVertex3f( P.X(), P.Y(), P.Z() );
    }
    glEnd();
}

void GravityThing::Keyboard( unsigned char key, int x, int y )
{
    PbaThingyDingy::Keyboard(key,x,y);
    if( key == 'v' ){ box->toggle_visible(); }
    if( key == 'w' ){ box->toggle_wireframe(); }
    if( key == 'e' ){ emit = !emit; }
    if( key == 'z' )
    {
        auto* f = dynamic_cast<GravityForce<DynamicalStateData>*>(gravityforce.get()); 
        Vector wind = f->get_gravity() + Vector(2,0,0);
        f->set_gravity(wind );
    }
    if( key == 'g' )
    {
        auto* f = dynamic_cast<GravityForce<DynamicalStateData>*>(gravityforce.get()); 
        f->set_gravity(f->get_gravity()/1.1);
        
    }
    if( key == 'G' )
    { 
        auto* f = dynamic_cast<GravityForce<DynamicalStateData>*>(gravityforce.get()); 
        f->set_gravity( f->get_gravity()*1.1 );
    }
    if( key == 'c' )
    {
        box->set_coeff_restitution( box->coeff_restitution()/1.1 );
        std::cout << "coefficient of restituion: " << box->coeff_restitution() << std::endl;
    }
    if( key == 'C' )
    { 
        box->set_coeff_restitution( box->coeff_restitution()*1.1 );
        std::cout << "coefficient of restituion: " << box->coeff_restitution() << std::endl;
    }
    if( key == 's' )
    {
        box->set_coeff_sticky( box->coeff_sticky()/1.1 );
        std::cout << "coefficient of sticky: " << box->coeff_sticky() << std::endl;
    }
    if( key == 'S' )
    { 
        box->set_coeff_sticky( box->coeff_sticky()*1.1 );
        std::cout << "coefficient of sticky: " << box->coeff_sticky() << std::endl;
    }
    if( key == 'l' )
    {
        solver.reset();
        b = CreateAdvanceVelocity<DynamicalStateData>(*state, *force, 1.0e9, 1.0e9);
        a = CreateAdvancePosition( *state, &collisions );
        solver = CreateForwardEulerSolver(a.get(),b.get()); 
        std::cout << "Using Leap Frog solver" << std::endl;
    }
    if( key == 'n' )
    {
        solver.reset();
        b = CreateAdvanceVelocity<DynamicalStateData>(*state, *force, 1.0e9, 1.0e9);
        a = CreateAdvancePosition( *state, &collisions );
        solver = CreateForwardEulerSolver(a.get(),b.get()); 
   std::cout << "Using Forward Euler solver" << std::endl;
    }
    if( key == 'b' )
    {
        solver.reset();
        b = CreateAdvanceVelocity<DynamicalStateData>(*state, *force, 1.0e9, 1.0e9);
        a = CreateAdvancePosition( *state, &collisions );
        solver = CreateForwardEulerSolver(a.get(),b.get()); 
        std::cout << "Using Backward Euler solver" << std::endl;
    }

}


void GravityThing::solve()
{
    if(emit)
    {
        emitter.emitCube(*state, 6, Vector(0,0,0));
    }
    solver->Solve(dt);
}

void GravityThing::Reset()
{
    state->Clear();
}

void GravityThing::Usage()
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

void GravityThing::AddCollisionSurface(CollisionSurface& s)
{
    std::cout << "Add CollisionSurface\n";
    if(box == nullptr)
    {
        box = std::make_unique<CollisionSurface>(s);
    }
    s.set_coeff_restitution(0.5);
    //s->set_coeff_sticky(0.1);
    collisions.set_collision_surface(box.get());
}


pba::PbaThing pba::CreateGravityThing(){ return PbaThing( new GravityThing() ); }


