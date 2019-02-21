#ifndef FLUX_H
#define FLUX_H

#include "../lattice.h"
#include "lattice_object.h"
#include "neighborlist.h"
#include "component.h"
#include "gaussian_noise.h"
#include "lattice_accessor.h"

class Flux1D;
class Flux2D;
class Flux3D;

typedef Factory<Flux1D, Dimensionality, Lattice*, const Real, Lattice_object<size_t>&, shared_ptr<Component>, shared_ptr<Component>, shared_ptr<Gaussian_noise>> Flux_factory;

class Flux1D {
public:
  Flux1D(Lattice*, const Real, Lattice_object<size_t>&, shared_ptr<Component>, shared_ptr<Component>, shared_ptr<Gaussian_noise>);
  virtual ~Flux1D();

  virtual int langevin_flux();

  enum error {
    ERROR_SIZE_INCOMPATIBLE,
    ERROR_NOT_IMPLEMENTED,
  };

  Lattice_object<Real> J_plus;
  Lattice_object<Real> J;


protected:
  void attach_neighborlists(shared_ptr<Neighborlist>, Dimension);
  int onsager_coefficient(Lattice_object<Real>&, Lattice_object<Real>&);
  int potential_difference(Lattice_object<Real>&, Lattice_object<Real>&);
  int langevin_flux(Dimension);
  //int mask(const stl::host_vector<int>&);

private:
  shared_ptr<Component> A;
  shared_ptr<Component> B;
  Lattice_object<Real> L;
  Lattice_object<Real> mu;
  Lattice_object<Real> t_L;
  Lattice_object<Real> t_mu;
  const Real D;
  shared_ptr<Gaussian_noise> gaussian;
};

class Flux2D : public Flux1D {
public:
  Flux2D(Lattice*, const Real, Lattice_object<size_t>&, shared_ptr<Component>, shared_ptr<Component>, shared_ptr<Gaussian_noise>);
  virtual ~Flux2D();

  virtual int langevin_flux() override;
};

class Flux3D : public Flux2D {
public:
  Flux3D(Lattice*, const Real, Lattice_object<size_t>&, shared_ptr<Component>, shared_ptr<Component>, shared_ptr<Gaussian_noise>);
  ~Flux3D();

  virtual int langevin_flux() override;
};

#endif