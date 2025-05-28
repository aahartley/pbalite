
#ifndef ____PBA_RAY_H____
#define ____PBA_RAY_H____

namespace pba
{
struct Ray 
{
  Ray(const Vector& o, const Vector& e, const Vector& d) 
  {
      origin = o;
      length = (e-o).magnitude();
      direction = d;
      inv_direction = Vector(1/d.X(), 1/d.Y(), 1/d.Z());
      sign[0] = (inv_direction.X() < 0);
      sign[1] = (inv_direction.Y() < 0);
      sign[2] = (inv_direction.Z() < 0);
  }
  Vector origin;
  double length;
  Vector direction;
  Vector inv_direction;
  int sign[3];
};
}
#endif