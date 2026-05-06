#include "boundary_conditions.h"

Register_class<Boundary1D, Boundary1D, Dimensionality, const Lattice_object<size_t>&, Boundary::Map> boundary_one_dimensions(one_D);
Register_class<Boundary1D, Boundary2D, Dimensionality, const Lattice_object<size_t>&, Boundary::Map> boundary_two_dimensions(two_D);
Register_class<Boundary1D, Boundary3D, Dimensionality, const Lattice_object<size_t>&, Boundary::Map> boundary_three_dimensions(three_D);

Boundary1D::Boundary1D(const Lattice_object<size_t>& mask, Boundary::Map boundary_type_for_dim) noexcept
: m_mask{mask}, m_x_neighborlist(mask),
  X_BOUNDARY_TYPE{
      boundary_type_for_dim[Dimension::X]
  }
{
    set_x_neighbors();
}

Boundary2D::Boundary2D(const Lattice_object<size_t>& mask, Boundary::Map boundary_type_for_dim) noexcept
: Boundary1D(mask, boundary_type_for_dim),
  m_y_neighborlist(mask),
  Y_BOUNDARY_TYPE {
      boundary_type_for_dim[Dimension::Y]
  }
{
    set_y_neighbors();
}

Boundary3D::Boundary3D(const Lattice_object<size_t>& mask, Boundary::Map boundary_type_for_dim) noexcept
: Boundary2D(mask, boundary_type_for_dim),
  m_z_neighborlist(mask),
  Z_BOUNDARY_TYPE{
      boundary_type_for_dim[Dimension::Z]
  }
{
    set_z_neighbors();
}

void Boundary1D::apply_boundary(stl::device_vector<Real>& input, Neighborlist& nlist) {

    Value_index_pair<Real> boundary {
        input,
        nlist.get_subject()
    };

    Value_index_pair<Real> system_edge {
        input,
        nlist.get_neighbors()
    };

    stl::copy(EXEC_PAR system_edge.begin(), system_edge.end(), boundary.begin());
}

void Boundary1D::zero_boundary(stl::device_vector<Real>& input, Neighborlist& nlist) {

    Value_index_pair<Real> boundary {
        input,
        nlist.get_subject()
    };

    stl::fill(EXEC_PAR boundary.begin(), boundary.end(), 0);
}

// sequential application: x first, then y, then z.
// this ensures correct cascading for edges and corners:
// after x boundaries are set, y boundaries read from already-updated x boundary cells,
// and z boundaries read from already-updated x and y boundary cells.

void Boundary1D::update_boundaries(stl::device_vector<Real>& input) {
    apply_boundary(input, m_x_neighborlist);
}

void Boundary1D::zero_boundaries(stl::device_vector<Real>& input) {
    zero_boundary(input, m_x_neighborlist);
}

void Boundary2D::update_boundaries(stl::device_vector<Real>& input) {
    Boundary1D::update_boundaries(input);
    apply_boundary(input, m_y_neighborlist);
}

void Boundary2D::zero_boundaries(stl::device_vector<Real>& input) {
    Boundary1D::zero_boundaries(input);
    zero_boundary(input, m_y_neighborlist);
}

void Boundary3D::update_boundaries(stl::device_vector<Real>& input) {
    Boundary2D::update_boundaries(input);
    apply_boundary(input, m_z_neighborlist);
}

void Boundary3D::zero_boundaries(stl::device_vector<Real>& input) {
    Boundary2D::zero_boundaries(input);
    zero_boundary(input, m_z_neighborlist);
}

void Boundary1D::set_x_neighbors() {

    std::map<Boundary::Type, size_t> offset;
    offset[Boundary::Type::MIRROR] = 1;
    offset[Boundary::Type::PERIODIC] = m_mask.MX;

    Neighborlist_config x0 {
        Dimension::X,
        Direction::plus,
        offset[X_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::x0_boundary, m_mask, std::placeholders::_1)
    };

    Neighborlist_config xm {
        Dimension::X,
        Direction::minus,
        offset[X_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::xm_boundary, m_mask, std::placeholders::_1)
    };

    m_x_neighborlist.register_config(x0);
    m_x_neighborlist.register_config(xm);
    m_x_neighborlist.build();
}

void Boundary2D::set_y_neighbors() {

    std::map<Boundary::Type, size_t> offset;
    offset[Boundary::Type::MIRROR] = 1;
    offset[Boundary::Type::PERIODIC] = m_mask.MY;

    Neighborlist_config y0 {
        Dimension::Y,
        Direction::plus,
        offset[Y_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::y0_boundary, m_mask, std::placeholders::_1)
    };

    Neighborlist_config ym {
        Dimension::Y,
        Direction::minus,
        offset[Y_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::ym_boundary, m_mask, std::placeholders::_1)
    };

    m_y_neighborlist.register_config(y0);
    m_y_neighborlist.register_config(ym);
    m_y_neighborlist.build();
}

void Boundary3D::set_z_neighbors() {

    std::map<Boundary::Type, size_t> offset;
    offset[Boundary::Type::MIRROR] = 1;
    offset[Boundary::Type::PERIODIC] = m_mask.MZ;

    Neighborlist_config z0 {
        Dimension::Z,
        Direction::plus,
        offset[Z_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::z0_boundary, m_mask, std::placeholders::_1)
    };

    Neighborlist_config zm {
        Dimension::Z,
        Direction::minus,
        offset[Z_BOUNDARY_TYPE],
        bind(&Lattice_object<size_t>::zm_boundary, m_mask, std::placeholders::_1)
    };

    m_z_neighborlist.register_config(z0);
    m_z_neighborlist.register_config(zm);
    m_z_neighborlist.build();
}
