#ifndef __PBA_DYNAMICALSTATE_H__
#define __PBA_DYNAMICALSTATE_H__

#include "Color.h"
#include "Vector.h"
#include "AABB.h"

#include <vector>
#include <string>
#include <map>
#include <memory>
#include <iostream>

namespace pba
{

template<typename T>
class DSAttribute
{
  public:
    DSAttribute() : name("unknown") {}
    DSAttribute(const std::string& nam, const T& def) : name(nam), defVal(def) {}
    ~DSAttribute(){}

      
    const size_t size() const { return data.size(); }
    const bool empty() const { return data.empty(); }
    void expand_to(size_t n)
    {
       if(data.size() >= n){ return; }
       size_t old_size = data.size();
       data.resize(n);
       for(size_t i=old_size; i<data.size(); i++)
       {
          data[i] = defVal;
       }
    }
    
    void set(size_t i, const T& value)  { data[i] = value; }
    const T& get(size_t i) const { return data[i]; }
    T& get(size_t i) { return data[i]; }
    void clear() { data.clear(); }
    const std::string& attr_name() const { return name; }
    const T& default_value() const { return defVal; }
    void erase(const size_t i)
    {
      typename std::vector<T>::iterator removable = data.begin() + i;
	    data.erase(removable);
    }
    //typename for ::
    typename std::vector<T>::const_iterator cbegin() const { return data.begin(); }
    typename std::vector<T>::const_iterator cend() const { return data.end(); }
    typename std::vector<T>::iterator begin() { return data.begin(); }
    typename std::vector<T>::iterator end() { return data.end(); }

  private:
    std::vector<T> data;
    std::string name;
    T defVal;
};

class DynamicalStateData
{
  public:
    DynamicalStateData(const std::string& nam = "DynamicDataNoName");
    ~DynamicalStateData() {}

    void create_attr(const std::string& nam, const int& def);
    void create_attr(const std::string& nam, const float& def);
    void create_attr(const std::string& nam, const Vector& def);
    void create_attr(const std::string& nam, const Color& def);
 
    // Add a single particle
    const size_t add();
    // Add many particles
    const size_t add(const size_t nb);
    size_t nb() const { return nb_items; }
    void clear();

    const int& get_int_attr(const std::string& nam, const size_t p) const;
    const float& get_float_attr(const std::string& nam, const size_t p) const;
    const Vector& get_vector_attr(const std::string& nam, const size_t p) const;
    const Color& get_color_attr(const std::string& nam, const size_t p) const;

    void set_attr(const std::string& nam, const size_t p, const int& value); 
    void set_attr(const std::string& nam, const size_t p, const float& value); 
    void set_attr(const std::string& nam, const size_t p, const Vector& value); 
    void set_attr(const std::string& nam, const size_t p, const Color& value); 
    
    const Vector& pos(const size_t p) const;
    const Vector& vel(const size_t p) const;
    const Vector& accel(const size_t p) const;
    const float& mass(const size_t p) const;
    const float& rad(const size_t p) const;
    const int& id(const size_t p) const;
    const Color& ci(const size_t p) const;

    void set_pos(const size_t p, const Vector& value);
    void set_vel(const size_t p, const Vector& value);
    void set_accel(const size_t p, const Vector& value);
    void set_mass(const size_t p, const float& value);
    void set_rad(const size_t p, const float& value);
    void set_id(const size_t p, const int& value);
    void set_ci(const size_t p, const Color& value);

    const std::string& Name() const { return name; }
    void re_find_main_attrs();

    int erase_outside_bounds( const Vector& llc, const Vector& urc );



  protected:
    std::string name;
    double time;
    size_t nb_items;

    std::map< std::string, DSAttribute<int> > int_attributes;
    std::map< std::string, DSAttribute<float> > float_attributes;
    std::map< std::string, DSAttribute<Vector> > vector_attributes;
    std::map< std::string, DSAttribute<Color> > color_attributes;

    std::map< std::string, DSAttribute<Vector> >::iterator    positions;
    std::map< std::string, DSAttribute<Vector> >::iterator    velocities;
    std::map< std::string, DSAttribute<Vector> >::iterator    accelerations;
    std::map< std::string, DSAttribute<float> >::iterator     masses;
    std::map< std::string, DSAttribute<float> >::iterator     radii;
    std::map< std::string, DSAttribute<int> >::iterator       ids;
    std::map< std::string, DSAttribute<Color> >::iterator     cis;
};
typedef std::shared_ptr<pba::DynamicalStateData> DynamicalState;
  
DynamicalState CreateDynamicalState(const std::string& nam = "DynamicalDataNoName");
AABB BoundingBox( const DynamicalState& d );

}

#endif