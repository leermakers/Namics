#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "factory.h"
#include "lattice_accessor.h"
#include "value_index_pair.h"
#include "neighborlist.h"
#include "lattice_object.h"
#include <vector>
#include <functional>
#include <map>
#ifdef PAR_MESODYN
  #include <thrust/copy.h>
#endif

#include "stl_typedef.h"

class Boundary1D;
class Boundary2D;
class Boundary3D;

namespace Boundary {

  enum class Type {
      MIRROR,
      PERIODIC,
  };

  typedef std::map< std::string, Boundary::Type> Adapter_type;

  static Boundary::Adapter_type Adapter;

  typedef std::map<Dimension, Boundary::Type> Map;
  typedef Factory_template<Boundary1D, Dimensionality, Lattice_object<size_t>&, Boundary::Map> Factory;

}

class Boundary1D {
  public:
    Boundary1D( Lattice_object<size_t>& mask, Boundary::Map boundary_type );
    virtual ~Boundary1D() { }
    void update_boundaries(stl::device_vector<Real>&);
    void zero_boundaries(stl::device_vector<Real>&);

  protected:

    Lattice_object<size_t>& m_mask;
    Neighborlist m_neighborlist;

    Boundary::Type X_BOUNDARY_TYPE;
    std::map<Boundary::Type, size_t> OFFSET;

    void set_x_neighbors();
};

class Boundary2D : public Boundary1D {
  public:
    Boundary2D(Lattice_object<size_t>& mask, Boundary::Map boundary_type_);
    virtual ~Boundary2D() { }


  private:
    Boundary::Type Y_BOUNDARY_TYPE;
    std::map<Boundary::Type, size_t> OFFSET;

    void set_y_neighbors();
};

class Boundary3D : public Boundary2D {
  public:
    Boundary3D(Lattice_object<size_t>& mask, Boundary::Map boundary_type_);
    ~Boundary3D() { }


  private:
    Boundary::Type Z_BOUNDARY_TYPE;
    std::map<Boundary::Type, size_t> OFFSET;

    void set_z_neighbors();
};

#endif