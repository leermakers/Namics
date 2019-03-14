#include "flux.h"


Register_class<IFlux, Flux1D, Dimensionality, Lattice*, const Real, Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>> flux_one_dimensions(one_D);
Register_class<IFlux, Flux2D, Dimensionality, Lattice*, const Real, Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>> flux_two_dimensions(two_D);
Register_class<IFlux, Flux3D, Dimensionality, Lattice*, const Real, Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>> flux_three_dimensions(three_D);

IFlux::IFlux(Lattice* lat_, shared_ptr<IComponent> A_, shared_ptr<IComponent> B_)
    : J(lat_), component_a{A_}, component_b{B_} {}

Flux1D::Flux1D(Lattice* Lat, const Real D, Lattice_object<size_t>& mask, shared_ptr<IComponent> A, shared_ptr<IComponent> B, shared_ptr<Gaussian_noise> gaussian)
    : IFlux(Lat, A, B), J_plus(Lat), L(Lat), mu(Lat), t_L(Lat), t_mu(Lat), D{D}, gaussian(gaussian)
  {
  
  Neighborlist_config configuration {
    Dimension::X,
    Direction::plus,
    1
  };

  shared_ptr<Neighborlist> x_neighborlist = make_shared<Neighborlist>(mask);
  x_neighborlist->register_config(configuration);
  x_neighborlist->build();
  attach_neighborlists(x_neighborlist, Dimension::X);
}

Flux2D::Flux2D(Lattice* Lat, const Real D, Lattice_object<size_t>& mask, shared_ptr<IComponent> A, shared_ptr<IComponent> B, shared_ptr<Gaussian_noise> gaussian)
    : Flux1D(Lat, D, mask, A, B, gaussian)
  {

  Neighborlist_config configuration {
    Dimension::Y,
    Direction::plus,
    1
  };

  shared_ptr<Neighborlist> y_neighborlist = make_shared<Neighborlist>(mask);
  y_neighborlist->register_config(configuration);
  y_neighborlist->build();
  attach_neighborlists(y_neighborlist, Dimension::Y);
}

Flux3D::Flux3D(Lattice* Lat, const Real D, Lattice_object<size_t>& mask, shared_ptr<IComponent> A, shared_ptr<IComponent> B, shared_ptr<Gaussian_noise> gaussian)
    : Flux2D(Lat, D, mask, A, B, gaussian)
  {
  
  Neighborlist_config configuration {
    Dimension::Z,
    Direction::plus,
    1
  };

  shared_ptr<Neighborlist> z_neighborlist = make_shared<Neighborlist>(mask);
  z_neighborlist->register_config(configuration);
  z_neighborlist->build();

  attach_neighborlists(z_neighborlist, Dimension::Z);
}

Flux1D::~Flux1D() {
}

Flux2D::~Flux2D() {
}

Flux3D::~Flux3D() {
}

void Flux1D::attach_neighborlists(shared_ptr<Neighborlist> neighborlist, Dimension dimension) {
  J.attach_neighborlist(neighborlist, dimension);
  J_plus.attach_neighborlist(neighborlist, dimension);
  L.attach_neighborlist(neighborlist, dimension);
  mu.attach_neighborlist(neighborlist, dimension);
  t_L.attach_neighborlist(neighborlist, dimension);
  t_mu.attach_neighborlist(neighborlist, dimension);
}

void Flux1D::flux() {

  //Zero (with bounds checking) vector J before use
  stl::fill(J.begin(), J.end(), 0);

  if (component_a->rho.size() != J.size()) {
    //already checked: A.alpha.size = B.alpha.size and A.rho.size = B.rho.size
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  std::thread onsager(&Flux1D::onsager_coefficient, this, std::ref(Flux1D::component_a->rho), std::ref(Flux1D::component_b->rho));
  std::thread potential(&Flux1D::potential_difference, this, std::ref(Flux1D::component_a->alpha), std::ref(Flux1D::component_b->alpha));
  onsager.join(); potential.join();

  langevin_flux(Dimension::X);
}

void Flux2D::flux() {
  Flux1D::flux();

  Flux1D::langevin_flux(Dimension::Y);
}

void Flux3D::flux() {
  Flux2D::flux();

  Flux1D::langevin_flux(Dimension::Z);
}

int Flux1D::onsager_coefficient(Lattice_object<Real>& A, Lattice_object<Real>& B) {

  if (A.size() != B.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }
  stl::transform(A.begin(), A.end(), B.begin(), L.begin(), stl::multiplies<Real>());

  return 0;
}

int Flux1D::potential_difference(Lattice_object<Real>& A, Lattice_object<Real>& B) {

  if (A.size() != B.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  stl::transform(A.begin(), A.end(), B.begin(), mu.begin(), stl::minus<Real>());

  gaussian->add_noise(mu);

  return 0;
}

int Flux1D::langevin_flux(Dimension dimension) {
    //J_plus[z] = -D * ((L[z] + L[z + jump]) * (mu[z + jump] - mu[z]))
    //J_minus[z] = -J_plus[z - jump] (substituted into equation below)
    //J = J_plus[z] + J_minus[z]

  stl::fill(J_plus.begin(), J_plus.end(), 0);

  stl::transform(mu.available_neighbors[dimension]->begin(), mu.available_neighbors[dimension]->end(), mu.available_sites->begin(), t_mu.available_sites->begin(), stl::minus<Real>());
  stl::transform(L.available_sites->begin(), L.available_sites->end(), L.available_neighbors[dimension]->begin(), t_L.available_sites->begin(), stl::plus<Real>());
  stl::transform(t_mu.begin(), t_mu.end(), t_L.begin(), J_plus.begin(), const_multiply_functor(-D) );

  stl::transform(J_plus.begin(), J_plus.end(), J.begin(), J.begin(), stl::plus<Real>());
  stl::transform(J.available_neighbors[dimension]->begin(), J.available_neighbors[dimension]->end(), J_plus.available_sites->begin(), J.available_neighbors[dimension]->begin(), stl::minus<Real>());
  
  return 0;
}