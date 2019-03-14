#include "boundary_conditions.h"
#include "lattice_accessor.h"
#include <algorithm>
#include <iostream>

Boundary::Adapter_type Adapter {
    {"mirror", Boundary::Type::MIRROR},
    {"periodic", Boundary::Type::PERIODIC}
};

Boundary1D::Boundary1D(Lattice_object<size_t>& mask, Boundary::Map boundary_type_for_dim)
: m_mask{mask}, m_neighborlist(mask),
  X_BOUNDARY_TYPE{
      boundary_type_for_dim[Dimension::X]
  }
{
    OFFSET[Boundary::Type::MIRROR] = 1;
    OFFSET[Boundary::Type::PERIODIC] = mask.MX;

    set_x_neighbors();
}

Boundary2D::Boundary2D(Lattice_object<size_t>& mask, Boundary::Map boundary_type_for_dim)
: Boundary1D(mask, boundary_type_for_dim),
  Y_BOUNDARY_TYPE {
      boundary_type_for_dim[Dimension::Y]
  }
{
    OFFSET[Boundary::Type::MIRROR] = 1;
    OFFSET[Boundary::Type::PERIODIC] = mask.MY;

    set_y_neighbors();
}

Boundary3D::Boundary3D(Lattice_object<size_t>& mask, Boundary::Map boundary_type_for_dim)
: Boundary2D(mask, boundary_type_for_dim),
  Z_BOUNDARY_TYPE{
      boundary_type_for_dim[Dimension::Z]
  }
{
    OFFSET[Boundary::Type::MIRROR] = 1;
    OFFSET[Boundary::Type::PERIODIC] = mask.MZ;
    set_z_neighbors();
}

void Boundary1D::update_boundaries(stl::device_vector<Real>& input) {
    
    Value_index_pair<Real> boundary {
        input,
        Boundary1D::m_neighborlist.get_subject()
    };

    Value_index_pair<Real> system_edge {
        input,
        Boundary1D::m_neighborlist.get_neighbors()
    };

    stl::copy(system_edge.begin(), system_edge.end(), boundary.begin());
}

void Boundary1D::zero_boundaries(stl::device_vector<Real>& input) {

    Value_index_pair<Real> boundary {
        input,
        Boundary1D::m_neighborlist.get_subject()
    };

    stl::fill(boundary.begin(), boundary.end(), 0);
}


void Boundary1D::set_x_neighbors() {

    Neighborlist_config x0 {
        Dimension::X,
        Direction::plus,
        OFFSET[X_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::x0_boundary, m_mask, std::placeholders::_1)
    };

    Neighborlist_config xm {
        Dimension::X,
        Direction::minus,
        OFFSET[X_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::xm_boundary, m_mask, std::placeholders::_1)
    };

    m_neighborlist.register_config(x0);
    m_neighborlist.register_config(xm);
    m_neighborlist.build();
}

void Boundary2D::set_y_neighbors() {

    Neighborlist_config y0 {
        Dimension::Y,
        Direction::plus,
        OFFSET[Y_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::y0_boundary, m_mask, std::placeholders::_1)
    };

    Neighborlist_config ym {
        Dimension::Y,
        Direction::minus,
        OFFSET[Y_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::ym_boundary, m_mask, std::placeholders::_1)
    };

    m_neighborlist.register_config(y0);
    m_neighborlist.register_config(ym);
    m_neighborlist.build();
}

void Boundary3D::set_z_neighbors() {

    Neighborlist_config z0 {
        Dimension::Z,
        Direction::plus,
        OFFSET[Z_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::z0_boundary, m_mask, std::placeholders::_1)
    };

    Neighborlist_config zm {
        Dimension::Z,
        Direction::minus,
        OFFSET[Z_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::zm_boundary, m_mask, std::placeholders::_1)
    };

    m_neighborlist.register_config(z0);
    m_neighborlist.register_config(zm);
    m_neighborlist.build();
}