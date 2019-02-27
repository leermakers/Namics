#ifndef SANITY_CHECKS_H
#define SANITY_CHECKS_H

#include <vector>
#include "../namics.h"
#include <iostream>
#include <assert.h>
#include "../tools.h"

#ifdef PAR_MESODYN
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#endif

#include "stl_typedef.h"

/*

   The overall design pattern of this is subject-observer, where the subject is the checkable (Checkable base class)
   and the observer is the sanity check (Sanity_check base class).
  
   Checkables contain a list (views_) of attached sanity checks that are executed whenever the checkable calls perform_checks() or the sanity check calls check().

*/

// Abstract base class for sanity checking oberser
class Sanity_check
{
  public:
    virtual void check() = 0;
};

template<class T> 
class Checkable
{
  protected:
    void set_checkable_data(T* data, size_t size) {
        m_checkable_data = data;
        m_size = size;
    }

  public:
    explicit Checkable(T* data, size_t size) 
    : m_checkable_data{data},
      m_size{size}
    { }

    virtual ~Checkable() { } // Derived class responsible for freeing data

    virtual void attach(Sanity_check *obs)
    {
        m_views.push_back(obs);
    }

    virtual void perform_checks()
    {
        for (size_t i = 0; i < m_views.size(); ++i)
          m_views[i]->check();
    }

    const T * get_checkable_data() {
        return m_checkable_data;
    }

    const T * get_checkable_data() const {
        return m_checkable_data;
    }

    size_t get_checkable_size() {
        return m_size;
    }

    size_t get_checkable_size() const {
        return m_size;
    }

  private:
    //Disable default constructor
    Checkable();

    T* m_checkable_data;
    size_t m_size;
    std::vector<Sanity_check*> m_views;
};

template<class T>
class Check_theta : public Sanity_check
{
    Checkable<T>* m_checkable;
    const Real m_theta;
    const int m_component_index;


  public:
    Check_theta(Checkable<T> *checkable, Real theta, int component_index=0)
    : m_checkable{checkable}, m_theta{theta}, m_component_index{component_index}
    {
        checkable->attach(this);
    }

    ~Check_theta() { }

     /* virtual */void check()
    {
        Real sum{0};
        #ifdef PAR_MESODYN
        const thrust::device_ptr<const Real> data = thrust::device_pointer_cast(m_checkable->get_checkable_data());
        #else
        const T* data = m_checkable->get_checkable_data();
        #endif
        const size_t size = m_checkable->get_checkable_size();
        
        #ifdef PAR_MESODYN
        sum = thrust::reduce(data, data+size);
        #else
        sum = std::accumulate(data, data+size, sum);
        #endif
        
        if (round(sum) < m_theta) {
            std::cerr << "Error! Theta for component " << m_component_index << " = "
                      << sum << " whereas " << m_theta << " was expected!" << std::endl;
        }
    }
};

template<class T>
class Check_between_zero_and_one : public Sanity_check
{
    Checkable<T>* m_checkable;
    int m_identifier=0;

  public:
    Check_between_zero_and_one(Checkable<T>* checkable, int identifier=0)
    : m_checkable{checkable}, m_identifier(identifier)
    {
        checkable->attach(this);
    }

    ~Check_between_zero_and_one() { }

     /* virtual */void check()
    {
        size_t errors{0};
        #ifdef PAR_MESODYN
        const thrust::device_ptr<const Real> data = thrust::device_pointer_cast(m_checkable->get_checkable_data());
        #else
        const T* data = m_checkable->get_checkable_data();
        #endif
        const size_t size = m_checkable->get_checkable_size();

        errors = stl::count_if(data, data+size, is_negative_functor());

        if (errors > 0) {
            std::cerr << "Found " << errors << " values < 0 || > 1 in identifier " << m_identifier << std::endl;
  }
    }
};

template<class T>
class Check_index_unity : public Sanity_check
{
    Checkable<T>* m_checkable;
    int m_identifier=0;
    vector<Checkable<T>*> m_coupled_checkables;

  public:
    explicit Check_index_unity(Checkable<T>* checkable, int identifier=0)
    : m_checkable{checkable}, m_identifier(identifier), m_coupled_checkables(0)
    {
        checkable->attach(this);
        m_coupled_checkables.push_back(checkable);
    }

    ~Check_index_unity() { }

    void register_checkable(Checkable<T>* checkable_) {
      m_coupled_checkables.push_back(checkable_);
    }

     /* virtual */void check()
    {
        size_t errors{0};
        const size_t size = m_checkable->get_checkable_size();
        typename stl::device_vector<T> sum(size, 0);

        for (auto checkable_data : m_coupled_checkables) {
          #ifdef PAR_MESODYN
          const thrust::device_ptr<const Real> data = thrust::device_pointer_cast(checkable_data->get_checkable_data());
          #else
          const T* data = checkable_data->get_checkable_data();
          #endif
          stl::transform(data, data+size, sum.begin(), sum.begin(), stl::plus<T>());
        }

        errors = stl::count_if(sum.begin(), sum.end(), is_not_unity_functor());

        if (errors > 0) {
            std::cerr << "Found " << errors << " indices where sum != 1 in identifier " << m_identifier << std::endl;
        }
    }

};

#endif