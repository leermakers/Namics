#ifndef MESODYNxH
#define MESODYNxH
#include "input.h"
#include "solve_scf.h"
#include "system.h"
#include "output.h"
#include "lattice.h"
#include "tools.h"
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
#include "mesodyn/component.h"
#include "mesodyn/flux.h"
#include "mesodyn/collection_procedures.h"


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
  Real* solve_crank_nicolson();
  void load_alpha(Real*, const size_t);
  void sanity_check();

  /* Initialization*/
  enum init {
    INIT_HOMOGENEOUS,
    INIT_FROMPRO,
    INIT_FROMVTK,
  };

  stl::device_vector<Real> rho;
  string read_filename;
  int initial_conditions();
  std::map<size_t, size_t> generate_pairs(size_t);
  Real boundaryless_volume;
  

  /* Helper class instances */
  shared_ptr<Gaussian_noise> gaussian;
  shared_ptr<Boundary1D> boundary;
  vector< shared_ptr<IComponent> > components;
  vector< shared_ptr<IFlux> > fluxes;

  /* Mesodyn specific output */
  ostringstream filename;
  int writes;
  void set_filename();

  unique_ptr<Norm_densities> norm_densities;
  unique_ptr<Order_parameter> order_parameter;


public:
  Mesodyn(int, vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<State*>, vector<Reaction*>, vector<Molecule*>, vector<System*>, vector<Solve_scf*>, string);
  ~Mesodyn();

  bool mesodyn();

  static std::vector<string> KEYS;
  static std::vector<string> PARAMETERS;
  static std::vector<string> VALUES;

  int norm_theta(vector< shared_ptr<IComponent> >&);
  Real calculate_order_parameter();

  /* Inputs / output class interface functions */

  int write_output(int);
  bool CheckInput();
  

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
