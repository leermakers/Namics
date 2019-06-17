#include "component.h"

IComponent::IComponent(Lattice* lat_, Lattice_object<Real>& rho_init)
    : Lat{lat_}, rho(rho_init), alpha(lat_)
    { }

Component::Component(Lattice* Lat, shared_ptr<Boundary1D> boundary, Lattice_object<Real> rho)
//TODO: fix alpha size.
    : IComponent(Lat, rho), boundary(boundary) {
  update_boundaries();
}

Component::~Component() {
}

/******* Interface *******/

int Component::update_density(const Lattice_object<Real>& J, int sign) {
  //Explicit update

  if (J.size() != rho.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  stl::transform(J.begin(), J.end(), rho.begin(), rho.begin(), saxpy_functor(sign) );

  return 0;
}

int Component::update_density(const Lattice_object<Real>& J1, Real ratio, int sign) {
  //Implicit update
  if (J1.size() != rho.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }
  // Rho <- A * J1 + Rho
  stl::transform(J1.previous_state().begin(), J1.previous_state().end(), rho.begin(), rho.begin(), saxpy_functor(sign*ratio) );
  stl::transform(J1.begin(), J1.end(), rho.begin(), rho.begin(), saxpy_functor((1.0-ratio)*sign) );

  return 0;
}

void Component::update_boundaries() {
  boundary->update_boundaries( alpha.m_data );
  boundary->update_boundaries( rho.m_data );
}

Real Component::theta() {
  //TODO: update once rho gets a neighborlist
  Real sum{0.0};
  boundary->zero_boundaries( rho.m_data );
  #ifdef PAR_MESODYN
  sum = stl::reduce(rho.begin(), rho.end(), sum);
  #else
  sum = stl::accumulate(rho.begin(), rho.end(), sum);
  #endif
  boundary->update_boundaries( rho.m_data );
  return sum;
}