//*******************************************************************
//
//   DynamicalState.h
//
//  Functioning base class to store attribute data
//
//
//
//*******************************************************************
#ifndef __PBA_DYNAMICALSTATE_H__
#define __PBA_DYNAMICALSTATE_H__

#include "Color.h"
#include "Vector.h"
#include "AABB.h"

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <memory>

namespace pba
{
/*!
  DSAttribute stores a vector of type T
 */
template<typename T>
class DSAttribute
{
  public:
    DSAttribute() : name_("unknown") {}
    DSAttribute(const std::string& name, const T& def_val) 
      : name_(name), 
        def_val_(def_val) {}
    ~DSAttribute(){}

    const size_t size() const { return data_.size(); }
    const bool empty() const { return data_.empty(); }
    void expand_to(size_t n)
    {
       if (data_.size() >= n) { return; }
       size_t old_size = data_.size();
       data_.resize(n);
       for (size_t i = old_size; i < data_.size(); ++i)
       {
          data_[i] = def_val_;
       }
    }
    
    void set(size_t i, const T& value)  { data_[i] = value; }
    const T& get(size_t i) const { return data_[i]; }
    T& get(size_t i) { return data_[i]; }
    void clear() { data_.clear(); }
    const std::string& attr_name() const { return name_; }
    const T& default_value() const { return def_val_; }
    void erase(const size_t i)
    {
      typename std::vector<T>::iterator removable = data_.begin() + i;
	    data_.erase(removable);
    }
    //typename for ::
    typename std::vector<T>::const_iterator cbegin() const { return data_.begin(); }
    typename std::vector<T>::const_iterator cend() const { return data_.end(); }
    typename std::vector<T>::iterator begin() { return data_.begin(); }
    typename std::vector<T>::iterator end() { return data_.end(); }

  private:
    std::vector<T> data_;
    std::string name_;
    T def_val_;
};

enum IntegrationStage {
    NONE            ,      
    VELOCITY_UPDATED,
    POSITION_UPDATED
};

/*!
  DynamicalStateData maps of DSAttributes
 */
//! TODO: Add copy and move constructors
class DynamicalStateData
{
  public:
    DynamicalStateData(const std::string& name = "DynamicDataNoName");
    ~DynamicalStateData() {}

    void create_attr(const std::string& name, int def);
    void create_attr(const std::string& name, float def);
    void create_attr(const std::string& name, const Vector& def);
    void create_attr(const std::string& name, const Color& def);
 
    //! Add a single particle
    size_t Add();
    //! Add many particles
    size_t Add(size_t nb);

    size_t nb() const { return nb_items_; }
    void Clear();

    int get_int_attr(const std::string& name, size_t p) const;
    float get_float_attr(const std::string& name, size_t p) const;
    const Vector& get_vector_attr(const std::string& name, size_t p) const;
    const Color& get_color_attr(const std::string& name, size_t p) const;

    void set_attr(const std::string& name, size_t p, int value); 
    void set_attr(const std::string& name, size_t p, float value); 
    void set_attr(const std::string& name, size_t p, const Vector& value); 
    void set_attr(const std::string& name, size_t p, const Color& value); 
    
    const Vector& pos(size_t p) const;
    const Vector& vel(size_t p) const;
    const Vector& accel(size_t p) const;
    float mass(size_t p) const;
    float rad(size_t p) const;
    int id(size_t p) const;
    const Color& ci(size_t p) const;

    void set_pos(size_t p, const Vector& value);
    void set_vel(size_t p, const Vector& value);
    void set_accel(size_t p, const Vector& value);
    void set_mass(size_t p, float value);
    void set_rad(size_t p, float value);
    void set_id(size_t p, int value);
    void set_ci(size_t p, const Color& value);

    const std::string& name() const { return name_; }
    void re_find_main_attrs();
    void set_stage(const IntegrationStage& s) { stage_ = s; }
    const IntegrationStage& get_stage() { return stage_; }
    //! erase particles outside of aabb and return the number
    int EraseOutsideOfBounds(const Vector& llc, const Vector& urc);

  protected:
    std::string name_;
    double time_;
    size_t nb_items_;
    IntegrationStage stage_;

    std::unordered_map< std::string, DSAttribute<int> > int_attributes_;
    std::unordered_map< std::string, DSAttribute<float> > float_attributes_;
    std::unordered_map< std::string, DSAttribute<Vector> > vector_attributes_;
    std::unordered_map< std::string, DSAttribute<Color> > color_attributes_;

    std::unordered_map< std::string, DSAttribute<Vector> >::iterator    positions_;
    std::unordered_map< std::string, DSAttribute<Vector> >::iterator    velocities_;
    std::unordered_map< std::string, DSAttribute<Vector> >::iterator    accelerations_;
    std::unordered_map< std::string, DSAttribute<float> >::iterator     masses_;
    std::unordered_map< std::string, DSAttribute<float> >::iterator     radii_;
    std::unordered_map< std::string, DSAttribute<int> >::iterator       ids_;
    std::unordered_map< std::string, DSAttribute<Color> >::iterator     cis_;
};
using DynamicalStatePtr = std::unique_ptr<pba::DynamicalStateData>;
  
DynamicalStatePtr CreateDynamicalState(const std::string& name = "DynamicalDataNoName");
AABB BoundingBox(const DynamicalStateData& d);

} //end of pba namespace

#endif //__PBA_DYNAMICALSTATE_H__