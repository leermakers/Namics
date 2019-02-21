#ifndef COMPONENT_H
#define COMPONENT_H

#include "lattice_object.h"
#include "boundary_conditions.h"
#include <memory>

enum error {
  ERROR_PERIODIC_BOUNDARY,
  ERROR_SIZE_INCOMPATIBLE,
  ERROR_NOT_IMPLEMENTED,
  ERROR_FILE_FORMAT,
};

class Component {
public:
  Component(Lattice*, shared_ptr<Boundary1D>, Lattice_object<Real>); //1D
  ~Component();

  Lattice_object<Real> rho;
  Lattice_object<Real> alpha;

  Lattice* Lat;

  int update_density(Lattice_object<Real>&, int = 1);     //Explicit scheme
  int update_density(Lattice_object<Real>&, Lattice_object<Real>&, Real ratio, int = 1); //Implicit scheme
  int update_boundaries();
  Real theta();

private:
  shared_ptr<Boundary1D> boundary;
};

#endif