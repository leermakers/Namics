#include "factory.h"
#include "flux.h"
#include "component.h"
#include "lattice_accessor.h"
#include "gaussian_noise.h"
#include "../lattice.h"

Register_class<IFlux, Flux1D, Dimensionality, Lattice*, const Real, Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>> flux_one_dimensions(one_D);
Register_class<IFlux, Flux2D, Dimensionality, Lattice*, const Real, Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>> flux_two_dimensions(two_D);
Register_class<IFlux, Flux3D, Dimensionality, Lattice*, const Real, Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>> flux_three_dimensions(three_D);
Register_class<Boundary1D, Boundary1D, Dimensionality, Lattice_object<size_t>&, Boundary::Map> boundary_one_dimensions(one_D);
Register_class<Boundary1D, Boundary2D, Dimensionality, Lattice_object<size_t>&, Boundary::Map> boundary_two_dimensions(two_D);
Register_class<Boundary1D, Boundary3D, Dimensionality, Lattice_object<size_t>&, Boundary::Map> boundary_three_dimensions(three_D);