#include "Vector.h"
#include "Color.h"
#include "PbaThing.h"
#include "PbaUtils.h"
#include "SoftBodyState.h"
#include "GISolver.h"
#include "ExplicitDynamics.h"
#include "ForceLibrary.h"
#include "CollisionHandler.h"
#include "ParticleEmitter.h"

namespace pba
{

class SoftBodyThing: public PbaThingyDingy
{
  public:

    // Feel free to customize the name of this thing.
    SoftBodyThing(const std::string nam = "Soft_Body_Thing");
    ~SoftBodyThing();

    //! Initialization, including GLUT initialization.
    //! Called once at the beginning.  Could be used
    //! to set up things once.
    void Init(const std::vector<std::string>& args);
   
    /////////////////////////////////////////////////////////////// 
    // CASCADING CALLBACK FUNCTIONS 
    // The methods below are called as part of a bigger set
    // of similar calls.  Most of the other calls take place
    // in the viewer portion of this project.
    ///////////////////////////////////////////////////////////////

    //! Implements a display event
    //! This is where you code the opengl calls to display 
    //! your system.
    void Display();

    //! Implements responses to keyboard events 
    //! This is called when you hit a key
    void Keyboard(unsigned char key, int x, int y);

    //! Implements simulator updates during an idle period
    //! This is where the update process is coded
    //! for your dynamics problem.
    void solve();

    //! Implements reseting parameters and/or state
    //! This is called when you hit the 'r' key
    void Reset();

    //! Displays usage information on stdout
    //! If you set up actions with the Keyboard()
    //! callback, you should include a statement 
    //! here as to what the keyboard option is.
    void Usage();
    void AddCollisionSurface(CollisionSurface& s);
    
  private:
    // flag for whether to create more particles
    bool emit;
    SoftBodyStatePtr state;
    GISolverPtr a;
    GISolverPtr b;
    GISolverPtr solver;
    ForceSBDPtr force;
    ForceSBDPtr gravityforce;
    ForceSBDPtr struts;
 
    ElasticCollisionHandler collisions;
    CollisionSurfacePtr box;
    BVHPtr bvh;
    ParticleEmitter emitter;

    std::vector< std::vector<Vector> > verts;
    std::vector< std::vector<Triangle> > tris;

};

pba::PbaThing CreateSoftBodyThing();

}





