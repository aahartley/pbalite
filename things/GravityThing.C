

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
    box = makeCollisionSurface();
    CollisionInfinitePlane bottom(Vector(0,1,0),Vector(0,-3,0));
    CollisionInfinitePlane top(Vector(0,-1,0),Vector(0,3,0));
    CollisionInfinitePlane right(Vector(-1,0,0),Vector(3,0,0));
    CollisionInfinitePlane left(Vector(1,0,0),Vector(-3,0,0));
    CollisionInfinitePlane front(Vector(0,0,-1),Vector(0,0,3)); //closest to screen
    CollisionInfinitePlane back(Vector(0,0,1),Vector(0,0,-3)); //(z points to screen)
    box->addPlane(bottom);
    box->addPlane(top);
    box->addPlane(right);
    box->addPlane(left);
    box->addPlane(back);
    box->addPlane(front);
    AddCollisionSurface(box);
    state = CreateDynamicalState("ParticleState");
    Reset();
    int inc = 1;
    state->add(inc);

    std::cout << "Emit: Total Points " << state->nb() << std::endl;
    emitter = ParticleEmitter();
    for(size_t i = 0; i < state->nb(); i++)
    {
        Vector p(0,2.9,0);
        Vector v(0,0,0);
        Color c(0,0,1,1);
        state->set_pos(i,p);
        state->set_vel(i,v);
        state->set_ci(i,c);
        state->set_rad(i, 0.5f); //0.075


    }

    force = CreateAccumulatingForce();

    gravityforce = CreateGravityForce(Vector(0,-9.81f,0));

    std::shared_ptr<AccumulatingForce> f = dynamic_pointer_cast<AccumulatingForce>(force); 
	f->add_force(gravityforce);
    //GISolver a = CreateAdvancePosition(state);
    GISolver a = CreateAdvancePositionColl(state, collisions);
    GISolver b = CreateAdvanceVelocity(state, force);
    solver = CreateForwardEulerSolver(a, b);
    std::cout << name << " constructed\n";

}

GravityThing::~GravityThing(){}

void GravityThing::Init( const std::vector<std::string>& args ) 
{
    SetSimulationTimestep(0.01);
}
    
void GravityThing::Display() 
{
    pba::Display(box);
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
        std::shared_ptr<GravityForce> f = dynamic_pointer_cast<GravityForce>(gravityforce); 
        Vector wind = f->get_gravity() + Vector(2,0,0);
        f->set_gravity(wind );
    }
    if( key == 'g' )
    {
        std::shared_ptr<GravityForce> f = dynamic_pointer_cast<GravityForce>(gravityforce); 
        f->set_gravity(f->get_gravity()/1.1);
        
    }
    if( key == 'G' )
    { 
        std::shared_ptr<GravityForce> f = dynamic_pointer_cast<GravityForce>(gravityforce); 
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
        GISolver solvera = CreateAdvancePositionColl( state, collisions );
        GISolver solverb = CreateAdvanceVelocity(state,force);
        solver = CreateLeapFrogSolver(solvera,solverb);
        std::cout << "Using Leap Frog solver" << std::endl;
    }
    if( key == 'n' )
    {
        GISolver solvera = CreateAdvancePositionColl( state, collisions );
        GISolver solverb = CreateAdvanceVelocity(state,force);
        solver = CreateForwardEulerSolver(solvera,solverb); // forward
   std::cout << "Using Forward Euler solver" << std::endl;
    }
    if( key == 'b' )
    {
        GISolver solverb = CreateAdvanceVelocity(state,force);
        GISolver solvera = CreateAdvancePositionColl( state, collisions );
        solver = CreateForwardEulerSolver(solverb,solvera); //backward
        std::cout << "Using Backward Euler solver" << std::endl;
    }

}


void GravityThing::solve()
{
    if(emit)
    {
        emitter.emitCube(state, 6, Vector(0,0,0));
    }
    solver->solve(dt);
}

void GravityThing::Reset()
{
    state->clear();
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

void GravityThing::AddCollisionSurface(pba::CollisionSurface& s)
{
    std::cout << "Add CollisionSurface\n";
    box = s;
    s->set_coeff_restitution(0.5);
    //s->set_coeff_sticky(0.1);
    collisions.set_collision_surface(box);
}


pba::PbaThing pba::CreateGravityThing(){ return PbaThing( new GravityThing() ); }


