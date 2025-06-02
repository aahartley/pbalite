//*******************************************************************
//
//   DynamicalState.C
//
//  Functioning base class to store attribute data
//
//
//
//*******************************************************************

#include "DynamicalState.h"
#include <iostream>

namespace pba
{

DynamicalStateData::DynamicalStateData(const std::string& name)
  : name_(name),
    time_(0.0),
    nb_items_(0),
    stage_(NONE)
{
    // create standard set of attributes:
    //   id, pos, vel, accel, ci, mass
    vector_attributes_["pos"] = DSAttribute<Vector>("pos", Vector(0,0,0));
    vector_attributes_["vel"] = DSAttribute<Vector>("vel", Vector(0,0,0));
    vector_attributes_["accel"] = DSAttribute<Vector>("accel", Vector(0,0,0));
    float_attributes_["mass"] = DSAttribute<float>("mass", 1.0);
    float_attributes_["rad"] = DSAttribute<float>("rad", 0.0);
    int_attributes_["id"] = DSAttribute<int>("id", -1);
    color_attributes_["ci"] = DSAttribute<Color>("ci", Color(1,1,1,0));
    re_find_main_attrs();
 
}

void DynamicalStateData::create_attr(const std::string& name, int def)
{
    //map returnns end() if key not found
    if (int_attributes_.find(name) != int_attributes_.end()) { return; }
    int_attributes_[name] = DSAttribute<int>(name, def);
    int_attributes_[name].expand_to(int_attributes_["id"].size()); //make new attr same size as rest of data
    re_find_main_attrs();

}

void DynamicalStateData::create_attr(const std::string& name, float def)
{
   if (float_attributes_.find(name) != float_attributes_.end()) { return; }
   float_attributes_[name] = DSAttribute<float>(name, def);
   float_attributes_[name].expand_to( int_attributes_["id"].size());
   re_find_main_attrs();

}

void DynamicalStateData::create_attr(const std::string& name, const Vector& def)
{
   if (vector_attributes_.find(name) != vector_attributes_.end()){ return; }
   vector_attributes_[name] = DSAttribute<Vector>(name, def);
   vector_attributes_[name].expand_to(int_attributes_["id"].size());
   re_find_main_attrs();

}

void DynamicalStateData::create_attr(const std::string& name, const Color& def)
{
   if (color_attributes_.find(name) != color_attributes_.end()){ return; }
   color_attributes_[name] = DSAttribute<Color>(name, def);
   color_attributes_[name].expand_to(int_attributes_["id"].size());
   re_find_main_attrs();

}

size_t DynamicalStateData::Add()
{
   size_t add_size = nb_items_ + 1;
   for (std::map<std::string,DSAttribute<int>>::iterator a = int_attributes_.begin(); a != int_attributes_.end(); ++a)
   {
      a->second.expand_to(add_size);
   }
   for (std::map<std::string,DSAttribute<float>>::iterator a = float_attributes_.begin(); a != float_attributes_.end(); ++a)
   {
      a->second.expand_to(add_size);
   }
   for (std::map<std::string,DSAttribute<Vector>>::iterator a = vector_attributes_.begin(); a != vector_attributes_.end(); ++a )
   {
      a->second.expand_to(add_size);
   }
   for (std::map<std::string,DSAttribute<Color>>::iterator a = color_attributes_.begin(); a != color_attributes_.end(); ++a)
   {
      a->second.expand_to(add_size);
   }
   nb_items_ = nb_items_ + 1;
   re_find_main_attrs();

   return add_size-1; // return the index of the new particle
}

size_t DynamicalStateData::Add(size_t nb)
{
   size_t add_size = nb_items_ + nb;
   for (std::map<std::string,DSAttribute<int>>::iterator a = int_attributes_.begin(); a != int_attributes_.end(); ++a)
   {
      a->second.expand_to(add_size);
   }
   for (std::map<std::string,DSAttribute<float>>::iterator a = float_attributes_.begin(); a != float_attributes_.end(); ++a)
   {
      a->second.expand_to(add_size);
   }
   for (std::map<std::string,DSAttribute<Vector>>::iterator a = vector_attributes_.begin(); a != vector_attributes_.end(); ++a)
   {
      a->second.expand_to(add_size);
   }
   for (std::map<std::string,DSAttribute<Color>>::iterator a = color_attributes_.begin(); a != color_attributes_.end(); ++a)
   {
      a->second.expand_to(add_size);
   }
   nb_items_ = nb_items_ + nb;
   re_find_main_attrs();

   return add_size-1; // return the index of the last particle
}

void DynamicalStateData::Clear()
{
   for (std::map<std::string, DSAttribute<int>>::iterator a = int_attributes_.begin(); a != int_attributes_.end(); ++a)
   {
      a->second.clear();
   }
   for (std::map<std::string, DSAttribute<float>>::iterator a = float_attributes_.begin(); a != float_attributes_.end(); ++a)
   {
      a->second.clear();
   }
   for (std::map<std::string, DSAttribute<Vector>>::iterator a = vector_attributes_.begin(); a != vector_attributes_.end(); ++a)
   {
      a->second.clear();
   }
   for (std::map<std::string, DSAttribute<Color>>::iterator a = color_attributes_.begin(); a != color_attributes_.end(); ++a)
   {
      a->second.clear();
   }

   time_ = 0.0;
   nb_items_ = 0.0;
}

int DynamicalStateData::get_int_attr(const std::string& name, size_t p) const
{
   std::map<std::string,DSAttribute<int>>::const_iterator a = int_attributes_.find(name);
   return a->second.get(p);
}

float DynamicalStateData::get_float_attr(const std::string& name, size_t p) const
{
   std::map<std::string,DSAttribute<float>>::const_iterator a = float_attributes_.find(name);
   return a->second.get(p);
}

const Vector& DynamicalStateData::get_vector_attr(const std::string& name, size_t p) const
{
   std::map<std::string,DSAttribute<Vector>>::const_iterator a = vector_attributes_.find(name);
   return a->second.get(p);
}

const Color& DynamicalStateData::get_color_attr(const std::string& name, size_t p) const
{
   std::map<std::string,DSAttribute<Color>>::const_iterator a = color_attributes_.find(name);
   return a->second.get(p);
}

    
int DynamicalStateData::id(size_t p) const
{
   return ids_->second.get(p);
}

float DynamicalStateData::mass(size_t p) const
{
   return masses_->second.get(p);
}

float DynamicalStateData::rad(size_t p) const
{
   return radii_->second.get(p);
}

const Vector& DynamicalStateData::pos(size_t p) const
{
   return positions_->second.get(p);
}

const Vector& DynamicalStateData::vel(size_t p) const
{
   return velocities_->second.get(p);
}

const Vector& DynamicalStateData::accel(size_t p) const
{
   return accelerations_->second.get(p);
}

const Color& DynamicalStateData::ci(size_t p) const
{
   return cis_->second.get(p);
}


void DynamicalStateData::set_attr(const std::string& name, size_t p, int value)
{
   int_attributes_[name].set(p, value);
}

void DynamicalStateData::set_attr(const std::string& name, size_t p, float value)
{
   float_attributes_[name].set(p, value);
}

void DynamicalStateData::set_attr(const std::string& name, size_t p, const Vector& value) 
{
   vector_attributes_[name].set(p, value);
}

void DynamicalStateData::set_attr(const std::string& name, size_t p, const Color& value) 
{
   color_attributes_[name].set(p, value);
}


void DynamicalStateData::set_id(size_t p, int value)
{
   ids_->second.set(p, value);
}

void DynamicalStateData::set_pos(size_t p, const Vector& value)
{
   positions_->second.set(p, value);
}

void DynamicalStateData::set_vel(size_t p, const Vector& value)
{
   velocities_->second.set(p, value);
}

void DynamicalStateData::set_accel(size_t p, const Vector& value)
{
   accelerations_->second.set(p, value);
}

void DynamicalStateData::set_ci(size_t p, const Color& value)
{
   cis_->second.set(p, value);
}

void DynamicalStateData::set_mass(size_t p, float value)
{
   masses_->second.set(p, value);
}

void DynamicalStateData::set_rad(size_t p, float value)
{
   radii_->second.set(p, value);
}

void DynamicalStateData::re_find_main_attrs()
{
   positions_ = vector_attributes_.find( "pos" );
   if( positions_ == vector_attributes_.end() )
   {
      std::cout << "ERROR could not find positions_\n";
   }
   velocities_ = vector_attributes_.find( "vel" );
   if( velocities_ == vector_attributes_.end() )
   {
      std::cout << "ERROR could not find velocities_\n";
   }
   accelerations_ = vector_attributes_.find( "accel" );
   if( accelerations_ == vector_attributes_.end() )
   {
      std::cout << "ERROR could not find accelerations_\n";
   }
   masses_ = float_attributes_.find( "mass" );
   if( masses_ == float_attributes_.end() )
   {
      std::cout << "ERROR could not find masses_\n";
   }
   ids_ = int_attributes_.find( "id" );
   if( ids_ == int_attributes_.end() )
   {
      std::cout << "ERROR could not find ids_\n";
   }
   cis_ = color_attributes_.find( "ci" );
   if( cis_ == color_attributes_.end() )
   {
      std::cout << "ERROR could not find cis_\n";
   }
   radii_ = float_attributes_.find( "rad" );
   if( radii_ == float_attributes_.end() )
   {
      std::cout << "ERROR could not find radii_\n";
   }
}

int DynamicalStateData::EraseOutsideOfBounds( const Vector& llc, const Vector& urc )
{
   AABB bounds(llc, urc);
   //std::cout << "erase outsdie of bounds" << ": llc " << llc.X() << ' ' << llc.Y() << ' ' << llc.Z() <<  " urc " << urc.X() << ' ' << urc.Y() << ' ' << urc.Z() <<'\n';
   size_t p = 0;
   int count = 0;
   while( p < nb_items_ )
   {
      const Vector& P = pos(p);
      if(bounds.IsInside(P))
      {
         p++;
      }
      else
      {
         for( std::map< std::string, DSAttribute<int> >::iterator a = int_attributes_.begin(); a != int_attributes_.end(); a++ )
         {
            a->second.erase(p);
         }
         for( std::map< std::string, DSAttribute<float> >::iterator a = float_attributes_.begin(); a != float_attributes_.end(); a++ )
         {
            a->second.erase(p);
         }
         for( std::map< std::string, DSAttribute<Vector> >::iterator a = vector_attributes_.begin(); a != vector_attributes_.end(); a++ )
         {
            a->second.erase(p);
         }
         for( std::map< std::string, DSAttribute<Color> >::iterator a = color_attributes_.begin(); a != color_attributes_.end(); a++ )
         {
            a->second.erase(p);
         }
         //std::cout << "erased particle: " << p << '\n';
         nb_items_--;
	 count++; 
      }
   }
   if(count > 0)
      std::cout << "Number removed: " << count << std::endl;
   re_find_main_attrs();
   return count;
}

DynamicalStatePtr CreateDynamicalState(const std::string& name)
{
    //return DynamicalState(new DynamicalStateData(name)); //2 memory alloc header
    return std::make_unique<DynamicalStateData>(name);
}

pba::AABB BoundingBox(const DynamicalStateData& d)
{
   Vector llc = d.pos(0);
   Vector urc = d.pos(0);
   for(size_t i = 1; i < d.nb(); i++)
   {
      Vector P = d.pos(i);
      if( P[0] < llc[0] ){ llc[0] = P[0]; }
      if( P[1] < llc[1] ){ llc[1] = P[1]; }
      if( P[2] < llc[2] ){ llc[2] = P[2]; }
      if( P[0] > urc[0] ){ urc[0] = P[0]; }
      if( P[1] > urc[1] ){ urc[1] = P[1]; }
      if( P[2] > urc[2] ){ urc[2] = P[2]; }
   }
   return pba::AABB(llc,urc);
}

}//end of pba namespace