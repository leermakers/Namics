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

class IComponent {
  public:

    IComponent(Lattice* lat_, Lattice_object<Real> rho_init)
    : Lat{lat_}, rho{rho_init}, alpha{lat_}
    { }

    Lattice* Lat;
    Lattice_object<Real> rho;
    Lattice_object<Real> alpha;

    //Where lattice objects are fluxes, either implicit or explicit (flux at k and k+1)
    virtual int update_density(Lattice_object<Real>&, int = 1) = 0;
    virtual int update_density(Lattice_object<Real>&, Lattice_object<Real>&, Real ratio, int = 1) = 0; //Implicit scheme

    virtual int update_boundaries() = 0;

    //Compute sum of densities, without system boundaries.
    virtual Real theta() = 0;
};

class Component : public IComponent {
public:
  Component(Lattice*, shared_ptr<Boundary1D>, Lattice_object<Real>); //1D
  ~Component();
  
  int update_density(Lattice_object<Real>&, int = 1) override;     //Explicit scheme
  int update_density(Lattice_object<Real>&, Lattice_object<Real>&, Real ratio, int = 1) override; //Implicit scheme
  int update_boundaries() override;
  Real theta() override;

private:
  shared_ptr<Boundary1D> boundary;
};

#endif