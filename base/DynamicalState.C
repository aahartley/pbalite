#include "DynamicalState.h"

using namespace pba;

DynamicalStateData::DynamicalStateData(const std::string& nam) :
 name(nam),
 time(0.0),
 nb_items(0)
{
    // create standard set of attributes:
    //   id, pos, vel, accel, ci, mass
    vector_attributes["pos"] = DSAttribute<Vector>( "pos", Vector(0,0,0) );
    vector_attributes["vel"] = DSAttribute<Vector>( "vel", Vector(0,0,0) );
    vector_attributes["accel"] = DSAttribute<Vector>( "accel", Vector(0,0,0) );
    float_attributes["mass"] = DSAttribute<float>( "mass", 1.0 );
    float_attributes["rad"] = DSAttribute<float>( "rad", 0.0 );
    int_attributes["id"] = DSAttribute<int>( "id", -1 );
    color_attributes["ci"] = DSAttribute<Color>( "ci", Color(1,1,1,0) );
    re_find_main_attrs();
 
}

void DynamicalStateData::create_attr(const std::string& nam, const int& def)
{
    //map returnns end() if key not found
    if(int_attributes.find(nam) != int_attributes.end()) { return; }
    int_attributes[nam] = DSAttribute<int>(nam, def);
    int_attributes[nam].expand_to(int_attributes["id"].size()); //make new attr same size as rest of data
    re_find_main_attrs();

}

void DynamicalStateData::create_attr( const std::string& nam, const float& def )
{
   if(float_attributes.find(nam) != float_attributes.end()){ return; }
   float_attributes[nam] = DSAttribute<float>(nam, def);
   float_attributes[nam].expand_to( int_attributes["id"].size());
   re_find_main_attrs();

}

void DynamicalStateData::create_attr( const std::string& nam, const Vector& def )
{
   if(vector_attributes.find(nam) != vector_attributes.end()){ return; }
   vector_attributes[nam] = DSAttribute<Vector>(nam, def);
   vector_attributes[nam].expand_to(int_attributes["id"].size());
   re_find_main_attrs();

}

void DynamicalStateData::create_attr( const std::string& nam, const Color& def )
{
   if(color_attributes.find(nam) != color_attributes.end()){ return; }
   color_attributes[nam] = DSAttribute<Color>(nam, def);
   color_attributes[nam].expand_to(int_attributes["id"].size());
   re_find_main_attrs();

}

const size_t DynamicalStateData::add()
{
   size_t add_size = nb_items + 1;
   for(std::map<std::string,DSAttribute<int>>::iterator a = int_attributes.begin(); a != int_attributes.end(); a++)
   {
      a->second.expand_to(add_size);
   }
   for(std::map<std::string,DSAttribute<float>>::iterator a = float_attributes.begin(); a != float_attributes.end(); a++)
   {
      a->second.expand_to(add_size);
   }
   for(std::map<std::string,DSAttribute<Vector>>::iterator a = vector_attributes.begin(); a != vector_attributes.end(); a++)
   {
      a->second.expand_to(add_size);
   }
   for(std::map<std::string,DSAttribute<Color>>::iterator a = color_attributes.begin(); a != color_attributes.end(); a++)
   {
      a->second.expand_to(add_size);
   }
   nb_items = nb_items + 1;
   re_find_main_attrs();

   return add_size-1; // return the index of the new particle
}

const size_t DynamicalStateData::add( const size_t nb )
{
   size_t add_size = nb_items + nb;
   for(std::map<std::string,DSAttribute<int>>::iterator a = int_attributes.begin(); a != int_attributes.end(); a++)
   {
      a->second.expand_to(add_size);
   }
   for(std::map<std::string,DSAttribute<float>>::iterator a = float_attributes.begin(); a != float_attributes.end(); a++)
   {
      a->second.expand_to(add_size);
   }
   for(std::map<std::string,DSAttribute<Vector>>::iterator a = vector_attributes.begin(); a != vector_attributes.end(); a++)
   {
      a->second.expand_to(add_size);
   }
   for(std::map<std::string,DSAttribute<Color>>::iterator a = color_attributes.begin(); a != color_attributes.end(); a++)
   {
      a->second.expand_to(add_size);
   }
   nb_items = nb_items + nb;
   re_find_main_attrs();

   return add_size-1; // return the index of the last particle
}

void DynamicalStateData::clear()
{
   for(std::map<std::string, DSAttribute<int>>::iterator a = int_attributes.begin(); a != int_attributes.end(); a++)
   {
      a->second.clear();
   }
   for(std::map<std::string, DSAttribute<float>>::iterator a = float_attributes.begin(); a != float_attributes.end(); a++)
   {
      a->second.clear();
   }
   for(std::map<std::string, DSAttribute<Vector>>::iterator a = vector_attributes.begin(); a != vector_attributes.end(); a++)
   {
      a->second.clear();
   }
   for(std::map<std::string, DSAttribute<Color>>::iterator a = color_attributes.begin(); a != color_attributes.end(); a++)
   {
      a->second.clear();
   }

   time = 0.0;
   nb_items = 0.0;
}

const int& DynamicalStateData::get_int_attr(const std::string& nam, const size_t p) const
{
   std::map<std::string,DSAttribute<int>>::const_iterator a = int_attributes.find(nam);
   return a->second.get(p);
}

const float& DynamicalStateData::get_float_attr(const std::string& nam, const size_t p) const
{
   std::map<std::string,DSAttribute<float>>::const_iterator a = float_attributes.find(nam);
   return a->second.get(p);
}

const Vector& DynamicalStateData::get_vector_attr(const std::string& nam, const size_t p) const
{
   std::map<std::string,DSAttribute<Vector>>::const_iterator a = vector_attributes.find(nam);
   return a->second.get(p);
}

const Color& DynamicalStateData::get_color_attr(const std::string& nam, const size_t p) const
{
   std::map<std::string,DSAttribute<Color>>::const_iterator a = color_attributes.find(nam);
   return a->second.get(p);
}

    
const int& DynamicalStateData::id(const size_t p) const
{
   return ids->second.get(p);
}

const float& DynamicalStateData::mass(const size_t p) const
{
   return masses->second.get(p);
}

const float& DynamicalStateData::rad(const size_t p) const
{
   return radii->second.get(p);
}

const Vector& DynamicalStateData::pos(const size_t p) const
{
   return positions->second.get(p);
}

const Vector& DynamicalStateData::vel(const size_t p) const
{
   return velocities->second.get(p);
}

const Vector& DynamicalStateData::accel(const size_t p) const
{
   return accelerations->second.get(p);
}

const Color& DynamicalStateData::ci(const size_t p) const
{
   return cis->second.get(p);
}


void DynamicalStateData::set_attr(const std::string& nam, const size_t p, const int& value)
{
   int_attributes[nam].set(p, value);
}

void DynamicalStateData::set_attr(const std::string& nam, const size_t p, const float& value)
{
   float_attributes[nam].set(p, value);
}

void DynamicalStateData::set_attr(const std::string& nam, const size_t p, const Vector& value) 
{
   vector_attributes[nam].set(p, value);
}

void DynamicalStateData::set_attr(const std::string& nam, const size_t p, const Color& value) 
{
   color_attributes[nam].set(p, value);
}


void DynamicalStateData::set_id(const size_t p, const int& value)
{
   ids->second.set(p, value);
}

void DynamicalStateData::set_pos(const size_t p, const Vector& value)
{
   positions->second.set(p, value);
}

void DynamicalStateData::set_vel(const size_t p, const Vector& value)
{
   velocities->second.set(p, value);
}

void DynamicalStateData::set_accel(const size_t p, const Vector& value)
{
   accelerations->second.set(p, value);
}

void DynamicalStateData::set_ci(const size_t p, const Color& value)
{
   cis->second.set(p, value);
}

void DynamicalStateData::set_mass(const size_t p, const float& value)
{
   masses->second.set(p, value);
}

void DynamicalStateData::set_rad(const size_t p, const float& value)
{
   radii->second.set(p, value);
}

void DynamicalStateData::re_find_main_attrs()
{
   positions = vector_attributes.find( "pos" );
   if( positions == vector_attributes.end() )
   {
      std::cout << "ERROR could not find positions\n";
   }
   velocities = vector_attributes.find( "vel" );
   if( velocities == vector_attributes.end() )
   {
      std::cout << "ERROR could not find velocities\n";
   }
   accelerations = vector_attributes.find( "accel" );
   if( accelerations == vector_attributes.end() )
   {
      std::cout << "ERROR could not find accelerations\n";
   }
   masses = float_attributes.find( "mass" );
   if( masses == float_attributes.end() )
   {
      std::cout << "ERROR could not find masses\n";
   }
   ids = int_attributes.find( "id" );
   if( ids == int_attributes.end() )
   {
      std::cout << "ERROR could not find ids\n";
   }
   cis = color_attributes.find( "ci" );
   if( cis == color_attributes.end() )
   {
      std::cout << "ERROR could not find cis\n";
   }
   radii = float_attributes.find( "rad" );
   if( radii == float_attributes.end() )
   {
      std::cout << "ERROR could not find radii\n";
   }
}

int DynamicalStateData::erase_outside_bounds( const Vector& llc, const Vector& urc )
{
   AABB bounds(llc, urc);
   //std::cout << "erase outsdie of bounds" << ": llc " << llc.X() << ' ' << llc.Y() << ' ' << llc.Z() <<  " urc " << urc.X() << ' ' << urc.Y() << ' ' << urc.Z() <<'\n';
   size_t p = 0;
   int count = 0;
   while( p < nb_items )
   {
      const Vector& P = pos(p);
      if(bounds.isInside(P))
      {
         p++;
      }
      else
      {
         for( std::map< std::string, DSAttribute<int> >::iterator a = int_attributes.begin(); a != int_attributes.end(); a++ )
         {
            a->second.erase(p);
         }
         for( std::map< std::string, DSAttribute<float> >::iterator a = float_attributes.begin(); a != float_attributes.end(); a++ )
         {
            a->second.erase(p);
         }
         for( std::map< std::string, DSAttribute<Vector> >::iterator a = vector_attributes.begin(); a != vector_attributes.end(); a++ )
         {
            a->second.erase(p);
         }
         for( std::map< std::string, DSAttribute<Color> >::iterator a = color_attributes.begin(); a != color_attributes.end(); a++ )
         {
            a->second.erase(p);
         }
         //std::cout << "erased particle: " << p << '\n';
         nb_items--;
	 count++; 
      }
   }
   if(count > 0)
      std::cout << "Number removed: " << count << std::endl;
   re_find_main_attrs();
   return count;
}

DynamicalState pba::CreateDynamicalState(const std::string& nam)
{
    //return DynamicalState(new DynamicalStateData(nam)); //2 memory alloc header
    return std::make_shared<DynamicalStateData>(nam);
}

pba::AABB pba::BoundingBox( const DynamicalState& d )
{
   Vector llc = d->pos(0);
   Vector urc = d->pos(0);
   for(size_t i = 1; i < d->nb(); i++)
   {
      Vector P = d->pos(i);
      if( P[0] < llc[0] ){ llc[0] = P[0]; }
      if( P[1] < llc[1] ){ llc[1] = P[1]; }
      if( P[2] < llc[2] ){ llc[2] = P[2]; }
      if( P[0] > urc[0] ){ urc[0] = P[0]; }
      if( P[1] > urc[1] ){ urc[1] = P[1]; }
      if( P[2] > urc[2] ){ urc[2] = P[2]; }
   }
   return pba::AABB(llc,urc);
}
