#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "lattice_accessor.h"
#include "value_index_pair.h"
#include "neighborlist.h"
#include "lattice_object.h"
#include <vector>
#include <functional>
#include <map>
#include <thrust/copy.h>

#include "stl_typedef.h"

enum class Boundary_type {
    MIRROR,
    PERIODIC,
};

using Boundary_map = std::map<Dimension, Boundary_type>;

class Boundary1D {
  public:
    Boundary1D(Lattice_object<size_t>& mask, Boundary_map boundary_type );
    virtual ~Boundary1D() { }
    void update_boundaries(stl::device_vector<Real>&);
    void zero_boundaries(stl::device_vector<Real>&);

  protected:    

    Lattice_object<size_t>& m_mask;
    Neighborlist m_neighborlist;

    Boundary_type X_BOUNDARY_TYPE;
    std::map<Boundary_type, size_t> OFFSET;

    void set_x_neighbors();
};

class Boundary2D : public Boundary1D {
  public:
    Boundary2D(Lattice_object<size_t>& mask, Boundary_map boundary_type_);
    virtual ~Boundary2D() { }

  private:
    Boundary_type Y_BOUNDARY_TYPE;
    std::map<Boundary_type, size_t> OFFSET;

    void set_y_neighbors();
};

class Boundary3D : public Boundary2D {
  public:
    Boundary3D(Lattice_object<size_t>& mask, Boundary_map boundary_type_);
    ~Boundary3D() { }

  private:
    Boundary_type Z_BOUNDARY_TYPE;
    std::map<Boundary_type, size_t> OFFSET;

    void set_z_neighbors();
};

#endif