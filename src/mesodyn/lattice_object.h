#ifndef LATTICE_OBJECT_H
#define LATTICE_OBJECT_H

#include <algorithm>
#include <iterator>
#include <cstdio>
#include <vector>
#include <memory>
#include "neighborlist.h"
#include "sanity_checks.h"
#include "../lattice.h"
#include "value_index_pair.h"
#include "lattice_accessor.h"
#include "memento.h"

#ifdef PAR_MESODYN
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#endif

#include "stl_typedef.h"

template <class T>
class Lattice_object : public Lattice_accessor, public Checkable<T>
{

private:
//disable default constructor, because we always need use the fixed size of Lat:
Lattice_object() = delete;
std::vector<shared_ptr<Neighborlist>> m_neighborlist;
shared_ptr< Memento<T> > saved_state;
void set_checkable_data_wrapper() {
      #ifdef PAR_MESODYN
      this->set_checkable_data( thrust::raw_pointer_cast(m_data.data()), system_size );
      #else
      this->set_checkable_data( m_data.data(), system_size );
      #endif
}

public:

//Stores the system minus the boundaries
stl::device_vector<T> m_data;
shared_ptr<Value_index_pair<T>> available_sites;

typedef std::map<Offset_map, shared_ptr<Value_index_pair<T>>> Neighborlist_map;

Neighborlist_map available_neighbors;
const Lattice* m_subject_lattice;

virtual ~Lattice_object() { }

explicit Lattice_object(const Lattice* Lat_, T init=0.0)
: Lattice_accessor{Lat_}, Checkable<T>{ (T*)this, system_size }, m_data(system_size, init), m_subject_lattice{ Lat_ }
    {
      //Because resizing m_data may have changed the memory
      this->set_checkable_data_wrapper();
    }


template <class OtherT>
Lattice_object(const Lattice* Lat_, const stl::device_vector<OtherT>& init)
    : Lattice_accessor{Lat_}, Checkable<T>{ (T*)this, system_size }, m_data(init), m_subject_lattice{ Lat_ }
    {
      //Because resizing m_data may have changed the memory
      set_checkable_data_wrapper();
    }


template <class OtherT>
Lattice_object(const Lattice_object<OtherT>& copy)
    : Lattice_accessor{copy.m_subject_lattice}, Checkable<T>{ (T*)&copy, copy.system_size },
      m_neighborlist{copy.get_neighborlists()}, m_data(system_size), m_subject_lattice{copy.m_subject_lattice}
    { 
      stl::copy(copy.begin(), copy.end(), this->begin()); 
      set_checkable_data_wrapper();
    }

void attach_neighborlist(shared_ptr<Neighborlist> neighborlist, Offset_map offset)

	{
      m_neighborlist.push_back( neighborlist );

      if (available_sites)
        assert_available_sites_equal(neighborlist);
  
      available_sites = make_shared<Value_index_pair<T>>(m_data, neighborlist->get_subject() );
      available_neighbors[offset] = make_shared<Value_index_pair<T>>(m_data, neighborlist->get_neighbors() );
  }

void assert_available_sites_equal(shared_ptr<Neighborlist> neighborlist) {
      auto sites = neighborlist->get_subject();
      if (
        !stl::equal(sites.begin(), sites.end(), available_sites->indices.begin())
      ) throw "Error, tried to attach incompatible neighborlist";
}


void load_array(T* data, size_t size)
    {
      assert( size == m_data.size() );
      #if defined(CUDA) && ! defined(PAR_MESODYN)
        TransferDataToHost( m_data.data(), data, this->size() );
      #else
        stl::copy(data, data + size, m_data.begin());
      #endif
    }

T& operator()(size_t x, size_t y, size_t z) {
    #ifdef PAR_MESODYN
    return m_data[index(x,y,z)];
    #else
    try {
      return m_data.at(index(x,y,z));
    } catch (const std::out_of_range& oor) {
      std::cerr << "Out of Range error in Lattice_object: " << oor.what() << '\n';
      exit(0);
    }
    #endif
  }

const T operator()(size_t x, size_t y, size_t z) const {
    #ifdef PAR_MESODYN
    return m_data[index(x,y,z)];
    #else
    try {
      return m_data.at(index(x,y,z));
    } catch (const std::out_of_range& oor) {
      std::cerr << "Out of Range error in Lattice_object: " << oor.what() << '\n';
      exit(0);
    }
    #endif
  }


T& operator[](size_t x) {
      #ifdef PAR_MESODYN
      return m_data[x];
      #else
      try {
        return m_data.at(x);
      } catch (const std::out_of_range& oor) {
        std::cerr << "Out of Range error in Lattice_object: " << oor.what() << '\n';
        exit(0);
      }
      #endif
    }


const T& operator[](size_t x) const {
      #ifdef PAR_MESODYN
      return m_data[x];
      #else
      try {
        return m_data.at(x);
      } catch (const std::out_of_range& oor) {
        std::cerr << "Out of Range error in Lattice_object: " << oor.what() << '\n';
        exit(0);
      }
      #endif
    }


template <class OtherT>
typename stl::device_vector<T>& operator=(stl::device_vector<OtherT>& rhs)
{
  assert(rhs.size() == m_data.size());
  m_data = rhs;
  set_checkable_data_wrapper();
  return m_data;
}

template <class OtherT>
Lattice_object<T>& operator=(Lattice_object<OtherT>& copy)
{
  assert(copy.size() == m_data.size());
  m_data = copy.m_data;
  m_neighborlist = copy.m_neighborlist;
  available_sites = copy.available_sites;
  available_neighbors = copy.available_neighbors;
  set_checkable_data_wrapper();
  return *this;
}

#ifdef PAR_Mesodyn

template <class OtherT>
typename thrust::device_vector<T>& operator=(std::vector<OtherT>& rhs)
{
    assert(rhs.size() == m_data.size());
    m_data = rhs;
    return m_data;
}


template <class OtherT>
typename thrust::device_vector<T>& operator=(thrust::host_vector<OtherT>& rhs)
{
    assert(rhs.size() == m_data.size());
    m_data = rhs;
    return m_data;
}
#endif


operator T *() {
  #ifdef PAR_MESODYN
    return thrust::raw_pointer_cast(m_data.data());
  #else
    return const_cast<T*>(m_data.data());
  #endif
}


typename stl::device_vector<T>::iterator begin() {
  return m_data.begin();
}


typename stl::device_vector<T>::iterator end() {
  return m_data.end();
}


typename stl::device_vector<T>::const_iterator begin() const {
  return m_data.begin();
}


typename stl::device_vector<T>::const_iterator end() const {
  return m_data.end();
}


size_t size() {
  return (system_size);
}

size_t size() const {
  return (system_size);
}


void clear() {
    m_data.clear();
}

T* data() {
#ifdef PAR_MESODYN
  return (T*)m_data;
#else
  return m_data.data();
#endif
}

std::vector<shared_ptr<Neighborlist>> get_neighborlists() const {
  return m_neighborlist;
}

void save_state() {
    saved_state = make_shared<Memento<T>>(m_data);
}

const stl::device_vector<T>& previous_state() {
  return saved_state->state;
}

const stl::device_vector<T>& previous_state() const {
  return saved_state->state;
}

void reinstate_previous_state() {
  m_data = saved_state->state;
  set_checkable_data_wrapper();
}

stl::host_vector<T> to_host() {
  stl::host_vector<T> host_m_data = m_data;
  return host_m_data;
}

};

#endif