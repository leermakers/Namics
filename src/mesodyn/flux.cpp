#include "flux.h"

Register_class<IFlux, Flux1D, Dimensionality, Lattice*, Real, const Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>> flux_one_dimensions(one_D);
Register_class<IFlux, Flux2D, Dimensionality, Lattice*, Real, const Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>> flux_two_dimensions(two_D);
Register_class<IFlux, Flux3D, Dimensionality, Lattice*, Real, const Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>> flux_three_dimensions(three_D);

IFlux::IFlux(Lattice* lat_, shared_ptr<IComponent> A_, shared_ptr<IComponent> B_)
    : J(lat_), m_lat{lat_}, component_a{A_}, component_b{B_} { }

ILangevin_flux::ILangevin_flux(Lattice* Lat, Real D_, shared_ptr<IComponent> A_, shared_ptr<IComponent> B_, shared_ptr<Gaussian_noise> gaussian_)
    : IFlux(Lat, A_, B_), L(Lat), mu(Lat), D{D_}, gaussian(gaussian_) { }

Flux1D::Flux1D(Lattice* Lat, Real D, const Lattice_object<size_t>& mask, shared_ptr<IComponent> A, shared_ptr<IComponent> B, shared_ptr<Gaussian_noise> gaussian)
    : ILangevin_flux(Lat, D, A, B, gaussian), J_plus(Lat),  t_L(Lat), t_mu(Lat)
{
  Neighborlist_config configuration {
    Dimension::X,
    Direction::plus,
    1
  };

  offset = configuration.get_offset();

  shared_ptr<Neighborlist> x_neighborlist = make_shared<Neighborlist>(mask);
  x_neighborlist->register_config(configuration);
  x_neighborlist->build();
  attach_neighborlists(x_neighborlist, offset);
}

Flux2D::Flux2D(Lattice* Lat, Real D, const Lattice_object<size_t>& mask, shared_ptr<IComponent> A, shared_ptr<IComponent> B, shared_ptr<Gaussian_noise> gaussian)
    : Flux1D(Lat, D, mask, A, B, gaussian)
{

  Neighborlist_config configuration {
    Dimension::Y,
    Direction::plus,
    1
  };

  offset = configuration.get_offset();

  shared_ptr<Neighborlist> y_neighborlist = make_shared<Neighborlist>(mask);
  y_neighborlist->register_config(configuration);
  y_neighborlist->build();
  attach_neighborlists(y_neighborlist, offset);
}

Flux3D::Flux3D(Lattice* Lat, Real D, const Lattice_object<size_t>& mask, shared_ptr<IComponent> A, shared_ptr<IComponent> B, shared_ptr<Gaussian_noise> gaussian)
    : Flux2D(Lat, D, mask, A, B, gaussian)
{
  
  Neighborlist_config configuration {
    Dimension::Z,
    Direction::plus,
    1
  };

  offset = configuration.get_offset();

  shared_ptr<Neighborlist> z_neighborlist = make_shared<Neighborlist>(mask);
  z_neighborlist->register_config(configuration);
  z_neighborlist->build();

  attach_neighborlists(z_neighborlist, offset);
}

Flux3D_extended_stencil::Flux3D_extended_stencil(Lattice* Lat, Real D, const Lattice_object<size_t>& mask, shared_ptr<IComponent> A, shared_ptr<IComponent> B, shared_ptr<Gaussian_noise> gaussian)
    : ILangevin_flux(Lat, D, A, B, gaussian), t_L(Lat), t_mu(Lat), t_J(Lat)
{
  for (short x = -1 ; x < 2 ; ++x)
    for ( short y = -1 ; y < 2 ; ++y )
      for ( short z = -1 ; z < 2 ; ++z ) {
        if ( x != 0 or y != 0 or z != 0) {
          shared_ptr<Neighborlist> neighborlist = make_shared<Neighborlist>(mask);

          Offset_map temp {
                {Dimension::X, x},
                {Dimension::Y, y},
                {Dimension::Z, z}
              };

          if (x == 1)
            forward_offsets.emplace_back( temp );
          else if (x == 0 and z == 1)
            forward_offsets.emplace_back( temp );
          else if (z == 0 and x == 0 and y == 1)
            forward_offsets.emplace_back( temp );
          else
            backward_offsets.emplace_back( temp );

          Neighborlist_config_extended config( temp, bind(&Lattice_object<size_t>::skip_bounds, mask, std::placeholders::_1) );

          neighborlist->register_config(config);
          neighborlist->build();
          J.attach_neighborlist(neighborlist, temp);
          L.attach_neighborlist(neighborlist, temp);
          mu.attach_neighborlist(neighborlist, temp);
          t_L.attach_neighborlist(neighborlist, temp);
          t_mu.attach_neighborlist(neighborlist, temp);
          t_J.attach_neighborlist(neighborlist, temp);
        } 
      }
}

Flux1D::~Flux1D() {
}

Flux2D::~Flux2D() {
}

Flux3D::~Flux3D() {
}

Flux3D_extended_stencil::~Flux3D_extended_stencil() {
}

int ILangevin_flux::onsager_coefficient(Lattice_object<Real>& A, Lattice_object<Real>& B) {

  if (A.size() != B.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }
  stl::transform(A.begin(), A.end(), B.begin(), L.begin(), stl::multiplies<Real>());

  return 0;
}

int ILangevin_flux::potential_difference(Lattice_object<Real>& A, Lattice_object<Real>& B) {

  if (A.size() != B.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  stl::transform(A.begin(), A.end(), B.begin(), mu.begin(), stl::minus<Real>());

  gaussian->add_noise(mu);

  m_lat->set_bounds((Real*)mu);

  return 0;
}

void Flux1D::attach_neighborlists(shared_ptr<Neighborlist> neighborlist, const Offset_map& offset_) {
  J.attach_neighborlist(neighborlist, offset_);
  J_plus.attach_neighborlist(neighborlist, offset_);
  L.attach_neighborlist(neighborlist, offset_);
  mu.attach_neighborlist(neighborlist, offset_);
  t_L.attach_neighborlist(neighborlist, offset_);
  t_mu.attach_neighborlist(neighborlist, offset_);
}

void Flux1D::flux() {

  //Zero (with bounds checking) vector J before use
  stl::fill(J.begin(), J.end(), 0.0);

  if (component_a->rho.size() != J.size()) {
    //already checked: A.alpha.size = B.alpha.size and A.rho.size = B.rho.size
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  onsager_coefficient(component_a->rho, component_b->rho);
  potential_difference(component_a->alpha, component_b->alpha);

  langevin_flux(Flux1D::offset);
}

void Flux2D::flux() {
  Flux1D::flux();

  Flux1D::langevin_flux(Flux2D::offset);
}

void Flux3D::flux() {
  Flux2D::flux();

  Flux1D::langevin_flux(Flux3D::offset);
}

void Flux3D_extended_stencil::flux() {
  //Zero (with bounds checking) vector J before use
  stl::fill(J.begin(), J.end(), 0.0);

  if (component_a->rho.size() != J.size()) {
    //already checked: A.alpha.size = B.alpha.size and A.rho.size = B.rho.size
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  onsager_coefficient(component_a->rho, component_b->rho);
  potential_difference(component_a->alpha, component_b->alpha);

  for (auto& offset : forward_offsets)
    langevin_flux_forward(offset);

  for (auto& offset : backward_offsets)
    langevin_flux_backward(offset);
}

int Flux3D_extended_stencil::langevin_flux_forward(const Offset_map& offset_) {
  stl::transform(mu.available_neighbors[offset_]->begin(), mu.available_neighbors[offset_]->end(), mu.available_sites->begin(), t_mu.available_sites->begin(), stl::minus<Real>());
  stl::transform(L.available_neighbors[offset_]->begin(), L.available_neighbors[offset_]->end(), L.available_sites->begin(), t_L.available_sites->begin(), stl::plus<Real>());
  stl::transform(t_mu.available_sites->begin(), t_mu.available_sites->end(), t_L.available_sites->begin(), t_J.available_sites->begin(), binary_norm_functor(-1.0/26.0*D));
  stl::transform(J.available_sites->begin(), J.available_sites->end(), t_J.available_sites->begin(), J.available_sites->begin(), stl::plus<Real>());

  return 0;
}

int Flux3D_extended_stencil::langevin_flux_backward(const Offset_map& offset_) {
  stl::transform(mu.available_neighbors[offset_]->begin(), mu.available_neighbors[offset_]->end(), mu.available_sites->begin(), t_mu.available_sites->begin(), reverse_minus_functor());
  stl::transform(L.available_neighbors[offset_]->begin(), L.available_neighbors[offset_]->end(), L.available_sites->begin(), t_L.available_sites->begin(), stl::plus<Real>());
  stl::transform(t_mu.available_sites->begin(), t_mu.available_sites->end(), t_L.available_sites->begin(), t_J.available_sites->begin(), binary_norm_functor(-1.0/26.0*D) );
  stl::transform(J.available_sites->begin(), J.available_sites->end(), t_J.available_sites->begin(), J.available_sites->begin(), stl::minus<Real>());

  return 0;
}

int Flux1D::langevin_flux(const Offset_map& offset_) {
    //J_plus[z] = -D * ((L[z] + L[z + jump]) * (mu[z + jump] - mu[z]))
    //J_minus[z] = -J_plus[z - jump] (substituted into equation below)
    //J = J_plus[z] + J_minus[z]

  stl::fill(J_plus.begin(), J_plus.end(), 0.0);

  stl::transform(mu.available_neighbors[offset_]->begin(), mu.available_neighbors[offset_]->end(), mu.available_sites->begin(), t_mu.available_sites->begin(), stl::minus<Real>());
  stl::transform(L.available_neighbors[offset_]->begin(), L.available_neighbors[offset_]->end(), L.available_sites->begin(), t_L.available_sites->begin(), stl::plus<Real>());
  stl::transform(t_mu.begin(), t_mu.end(), t_L.begin(), J_plus.begin(), binary_norm_functor(-1.0/6.0*D) );

  stl::transform(J_plus.begin(), J_plus.end(), J.begin(), J.begin(), stl::plus<Real>());
  stl::transform(J.available_neighbors[offset_]->begin(), J.available_neighbors[offset_]->end(), J_plus.available_sites->begin(), J.available_neighbors[offset_]->begin(), stl::minus<Real>());
  
  return 0;
}