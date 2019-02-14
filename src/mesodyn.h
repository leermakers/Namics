#ifndef MESODYNxH
#define MESODYNxH
#include "input.h"
#include "namics.h"
#include "solve_scf.h"
#include "system.h"
#include "output.h"
#include "lattice.h"
#include "tools.h"
#include <random>     // noise generation
#include <ctime>      // output labling
#include <cassert>    // making sure all underlying assumptions (e.g. lattice geometry) are satisfied
#include <functional> // boundary conditions
#include <algorithm>  // transform, copy, find functions
#include <limits>     // output
#include <sstream>
#include <unistd.h>   // output
#include <memory>
#include "mesodyn/lattice_object.h"
#include "mesodyn/boundary_conditions.h"

#ifdef PAR_MESODYN
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/device_ptr.h>
#include <thrust/count.h>
#endif

#include "mesodyn/stl_typedef.h"

class Boundary1D;

enum error {
  ERROR_PERIODIC_BOUNDARY,
  ERROR_SIZE_INCOMPATIBLE,
  ERROR_NOT_IMPLEMENTED,
  ERROR_FILE_FORMAT,
};

class Gaussian_noise {
  //Makes sure that we keep generating new numbers, instead of the same over and over.

public:
  Gaussian_noise(shared_ptr<Boundary1D>, Real, int, Real, Real); // Not seeded (32 bits of randomness)
  Gaussian_noise(shared_ptr<Boundary1D>, Real, int, Real, Real, size_t); // Seeded
  int generate(size_t);
  int add_noise(stl::device_vector<Real>&);
  int add_noise(Lattice_object<Real>&);
  stl::device_vector<Real> noise;

private:
  seed_seq seed;
  mt19937 prng;
  normal_distribution<Real> dist;
  shared_ptr<Boundary1D> boundary;
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
  Lattice_object<Real> J_minus;
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

class Mesodyn : public Lattice_accessor {

private:
  /* Constructor arguments*/
  const string name;
  const vector<Input*> In;
  const vector<Lattice*> Lat;
  const vector<Molecule*> Mol;
  const vector<Segment*> Seg;
  const vector<State*> Sta;
	const vector<Reaction*> Rea;
  const vector<System*> Sys;
  const vector<Solve_scf*> New;
  vector <Output*> Out;
  const string brand;

  std::vector<string> KEYS;
  std::vector<string> PARAMETERS;
  std::vector<string> VALUES;

  bool input_success;

  /* Read from file */
  const Real D; // diffusionconstant
  const Real dt;
  const Real mean; // mean of gaussian noise (should be 0)
  const Real stddev; // stdev of gaussian noise (should be 1*D)
  const Real seed;  // seed of gaussian noise
  const bool seed_specified;
  const int timesteps; // length of the time evolution
  const int timebetweensaves; // how many timesteps before mesodyn writes the current variables to file
  const Real cn_ratio; // how much of the old J gets mixed in the crank-nicolson scheme
  int initialization_mode;
  const size_t component_no; // number of components in the system, read from SysMonMolList


  /* Flow control */
  Real* solve_explicit();
  void explicit_start();
  int noise_flux();
  Real* solve_crank_nicolson();
  void load_alpha(Real*, const size_t);
  void sanity_check();
  Real calculate_order_parameter();

  /* Initialization*/
  enum init {
    INIT_HOMOGENEOUS,
    INIT_FROMPRO,
    INIT_FROMVTK,
  };
  Real system_volume;
  stl::device_vector<Real> rho;
  vector<string> tokenize(string, char);
  string read_filename;
  int initial_conditions();
  vector<Lattice_object<Real>> init_rho_homogeneous();
  int norm_density(vector<Real>& rho, Real theta);
  void set_update_lists();
  vector<vector<int>> update_plus;
  vector<vector<int>> update_minus;
  Real boundaryless_volume;
  std::map<shared_ptr<Component>, Real> theta;

  /* Helper class instances */
  shared_ptr<Gaussian_noise> gaussian;
  shared_ptr<Boundary1D> boundary;
  vector< shared_ptr<Component> > component;
  vector< shared_ptr<Component> > solver_component;
  vector< unique_ptr<Flux1D> > flux;
  vector< unique_ptr<Flux1D> > solver_flux;

  /* Mesodyn specific output */
  ostringstream filename;
  int writes;
  bool write_vtk;
  void set_filename();
  Real order_parameter;

  /* Mathematics */
  int factorial (int);
  int combinations (int, int);


public:
  Mesodyn(int, vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<State*>, vector<Reaction*>, vector<Molecule*>, vector<System*>, vector<Solve_scf*>, string);
  ~Mesodyn();

  bool mesodyn();

  int norm_theta(vector< shared_ptr<Component> >&);

  Real lost;

  /* Inputs / output class interface functions */

  int write_output(int);
  bool CheckInput(const int);


  //Const-correct way of initializing member variables from file.
  template<typename Datatype>
  Datatype initialize(string option, Datatype default_value) {

    for (size_t i = 0 ; i < Mesodyn::PARAMETERS.size(); ++i)
		  if (option==Mesodyn::PARAMETERS[i]) {
          Datatype value;
          std::istringstream buffer{VALUES[i]};
          buffer >> value;
          return value;
      }
     return default_value;
  }
};
#endif
