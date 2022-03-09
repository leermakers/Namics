#ifndef COMPONENT_H
#define COMPONENT_H

#include "lattice_object.h"
#include "boundary_conditions.h"
#include <memory>
#include <thread>

enum error {
  ERROR_PERIODIC_BOUNDARY,
  ERROR_SIZE_INCOMPATIBLE,
  ERROR_NOT_IMPLEMENTED,
  ERROR_FILE_FORMAT,
};

class IComponent {
  public:

    IComponent(Lattice*, Lattice_object<Real>&);
    virtual ~IComponent() {} // not responsible for deleting Lattice pointer

    Lattice* Lat;
    Lattice_object<Real> rho;
    Lattice_object<Real> alpha;

    //Where lattice objects are fluxes, either implicit or explicit (flux at k and k+1)
    virtual int update_density(const Lattice_object<Real>&, int) = 0;
    virtual int update_density(const Lattice_object<Real>&, Real ratio, int=0) = 0; //Implicit scheme

    virtual void update_boundaries() = 0;

    //Compute sum of densities, without system boundaries.
    virtual Real theta() = 0;
};

class Component : public IComponent {
public:
  Component(Lattice*, shared_ptr<Boundary1D>, Lattice_object<Real>); //1D
  virtual ~Component();
  
  int update_density(const Lattice_object<Real>&, int) override;     //Explicit scheme
  int update_density(const Lattice_object<Real>&, Real ratio, int=0) override; //Implicit scheme
  void update_boundaries() override;
  std::thread par_update_boundaries() {
    return std::thread(&Component::update_boundaries, this);
  }
  Real theta() override;

private:
  shared_ptr<Boundary1D> boundary;
};

#endif