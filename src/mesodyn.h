#ifndef MESODYNxH
#define MESODYNxH
#include "input.h"
#include "namics.h"
#include "solve_scf.h"
#include "system.h"
#include "output.h"
#include "lattice.h"
#include <random>     // noise generation
#include <ctime>      // output labling
#include <cassert>    // making sure all underlying assumptions (e.g. lattice geometry) are satisfied
#include <functional> // boundary conditions
#include <algorithm>  // transform, copy, find functions
#include <limits>     // output
#include <unistd.h>   // output
#include <memory>

class Boundary1D;

enum error {
  ERROR_PERIODIC_BOUNDARY,
  ERROR_SIZE_INCOMPATIBLE,
  ERROR_NOT_IMPLEMENTED,
  ERROR_FILE_FORMAT,
};

class Lattice_Access {
public:

  Lattice_Access(Lattice*);
  ~Lattice_Access();

  template <typename T>
  inline T val(vector<T>& v, int x, int y, int z) {
      return v[x * JX + y * JY + z * JZ];
  }

  template <typename T>
  inline T* val_ptr(vector<T>& v, int x, int y, int z) {
      return &v[x * JX + y * JY + z * JZ];
  }

  inline int index(int x, int y, int z);
  inline vector<int> coordinate(int n);
  inline void skip_bounds( function<void(int,int,int)> );
  inline void par_skip_bounds( function<void(int,int,int)> );
  inline void bounds( function<void(int,int,int)> );
  inline void x0_boundary( function<void(int, int, int)>);
  inline void xm_boundary( function<void(int, int, int)>);
  inline void y0_boundary( function<void(int, int, int)>);
  inline void ym_boundary( function<void(int, int, int)>);
  inline void z0_boundary( function<void(int, int, int)>);
  inline void zm_boundary( function<void(int, int, int)>);

  const int dimensions;

private:
  const int JX;
  const int JY;
  const int JZ;
  int setMY(Lattice*);
  int setMZ(Lattice*);

protected:
  const int M;
  const int MX;
  const int MY;
  const int MZ;
};

class Gaussian_noise {
  //Makes sure that we keep generating new numbers, instead of the same over and over.

public:
  Gaussian_noise(Boundary1D*, Real, int, Real, Real); // Not seeded (32 bits of randomness)
  Gaussian_noise(Boundary1D*, Real, int, Real, Real, size_t); // Seeded
  int generate();
  int add_noise(vector<Real>&);
  vector<Real> noise;

private:
  seed_seq seed;
  mt19937 prng;
  normal_distribution<Real> dist;
  unique_ptr<Boundary1D> boundary;
};

class Boundary1D : protected Lattice_Access {
public:

  enum boundary {
    MIRROR,
    PERIODIC,
  };

  Boundary1D(Lattice*, boundary, boundary); //1D
  Boundary1D(Lattice*, vector<Real>&, boundary, boundary); //1D
  virtual ~Boundary1D();

  virtual int update_boundaries(vector<Real>&);

private:
  function<void(vector<Real>&)> bX0;
  function<void(vector<Real>&)> bXm;

  int set_x_boundaries(boundary, boundary, Real = 0);
  void bX0Mirror(vector<Real>&);
  void bXmMirror(vector<Real>&);
  void bXPeriodic(vector<Real>&);

};

class Boundary2D: public Boundary1D {
public:
  Boundary2D(Lattice*, boundary, boundary, boundary, boundary);
  virtual ~Boundary2D();

  virtual int update_boundaries(vector<Real>&) override;

private:
  function<void(vector<Real>&)> bY0;
  function<void(vector<Real>&)> bYm;

  int set_y_boundaries(boundary, boundary, Real = 0);
  void bY0Mirror(vector<Real>&);
  void bYmMirror(vector<Real>&);
  void bYPeriodic(vector<Real>&);

};

class Boundary3D: public Boundary2D {
public:
  Boundary3D(Lattice*, boundary, boundary, boundary, boundary, boundary, boundary);
  ~Boundary3D();

  virtual int update_boundaries(vector<Real>&) override;

private:
  function<void(vector<Real>&)> bZ0;
  function<void(vector<Real>&)> bZm;

  int set_z_boundaries(boundary, boundary, Real = 0);
  void bZ0Mirror(vector<Real>&);
  void bZmMirror(vector<Real>&);
  void bZPeriodic(vector<Real>&);
};

class Component : protected Lattice_Access {
public:
  Component(Lattice*, Boundary1D*, vector<Real>&); //1D
  ~Component();

  vector<Real> rho;
  vector<Real> alpha;

  Gaussian_noise* gaussian;

  Real rho_at(int, int, int);
  Real alpha_at(int, int, int);
  int update_density(vector<Real>&, int = 1);     //Explicit scheme
  int update_density(vector<Real>&, vector<Real>&, Real ratio, int = 1); //Implicit scheme
  int load_alpha(vector<Real>&);
  int load_rho(vector<Real>&);
  int update_boundaries();
  Real theta();

private:
  unique_ptr<Boundary1D> boundary;
};

class Flux1D : protected Lattice_Access {
public:
  Flux1D(Lattice*, Gaussian_noise*, Real, vector<int>&, Component*, Component*);
  virtual ~Flux1D();

  virtual int langevin_flux();

  Real J_at(int, int, int);
  Real L_at(int, int, int);
  Real mu_at(int, int, int);

  enum error {
    ERROR_SIZE_INCOMPATIBLE,
    ERROR_NOT_IMPLEMENTED,
  };

  vector<Real> J_plus;
  vector<Real> J_minus;
  vector<Real> J;

  shared_ptr<Gaussian_noise> gaussian;


protected:
  int onsager_coefficient(vector<Real>&, vector<Real>&);
  int potential_difference(vector<Real>&, vector<Real>&);
  int langevin_flux(vector<int>&, vector<int>&, int);
  int mask(vector<int>&);
  unique_ptr<Component> A;
  unique_ptr<Component> B;

  vector<Real> L;
  vector<Real> mu;
  const Real D;
  const int JX;
  vector<int> Mask_plus_x;
  vector<int> Mask_minus_x;
};

class Flux2D : public Flux1D {
public:
  Flux2D(Lattice*, Gaussian_noise*, Real, vector<int>&, Component*, Component*);
  virtual ~Flux2D();

  virtual int langevin_flux() override;


protected:
  int mask(vector<int>&);
  const int JY;
  vector<int> Mask_plus_y;
  vector<int> Mask_minus_y;

};

class Flux3D : public Flux2D {
public:
  Flux3D(Lattice*, Gaussian_noise*, Real, vector<int>&, Component*, Component*);
  ~Flux3D();

  virtual int langevin_flux() override;

private:
  int mask(vector<int>&);

protected:
  const int JZ;
  vector<int> Mask_plus_z;
  vector<int> Mask_minus_z;
};

class Mesodyn : private Lattice_Access {

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

  /* Read from file */
  Real D; // diffusionconstant
  Real dt;
  Real mean; // mean of gaussian noise (should be 0)
  Real stddev; // stdev of gaussian noise (should be 1*D)
  Real seed;  // seed of gaussian noise
  bool seed_specified;
  int timesteps; // length of the time evolution
  int timebetweensaves; // how many timesteps before mesodyn writes the current variables to file
  int initialization_mode;
  const size_t component_no; // number of components in the system, read from SysMonMolList


  /* Flow control */
  int RC;
  Real* solve_explicit();
  void get_mask();
  void explicit_start();
  int noise_flux();
  Real* solve_crank_nicolson();
  void load_alpha(vector<Real>&, size_t);
  int sanity_check();
  Real calculate_order_parameter();
  Real cn_ratio; // how much of the old J gets mixed in the crank-nicolson scheme

  /* Initialization*/
  enum init {
    INIT_HOMOGENEOUS,
    INIT_FROMPRO,
    INIT_FROMVTK,
  };
  Real system_volume;
  vector<Real> rho;
  vector<string> tokenize(string, char);
  string read_filename;
  int initial_conditions();
  vector<Real>&  flux_callback(int);
  int init_rho_homogeneous(vector< vector<Real> >&, vector<int>&);
  int norm_density(vector<Real>& rho, Real theta);
  void set_update_lists();
  vector<vector<int>> update_plus;
  vector<vector<int>> update_minus;
  Real boundaryless_volume;

  /* Helper class instances */
  Boundary1D* boundary;
  vector<Component*> component;
  vector<Component*> solver_component;
  vector<Flux1D*> flux;
  vector<Flux1D*> solver_flux;

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
  Mesodyn(vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<State*>, vector<Reaction*>, vector<Molecule*>, vector<System*>, vector<Solve_scf*>, string);
  ~Mesodyn();

  bool mesodyn();

  int norm_theta(vector<Component*>&);

  Real lost;


  /* Inputs / output class interface functions */
  vector<string> ints;
  vector<string> Reals;
  vector<string> bools;
  vector<string> strings;
  vector<Real> Reals_value;
  vector<int> ints_value;
  vector<bool> bools_value;
  vector<string> strings_value;
  int GetValue(string, int&, Real&, string&);
  int write_output();

  std::vector<string> KEYS;
  std::vector<string> PARAMETERS;
  std::vector<string> VALUES;
  bool CheckInput(int);
  string GetValue(string);
};
#endif
