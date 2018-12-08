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

#ifdef PAR_MESODYN
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/device_ptr.h>
namespace stl = thrust;
#else
namespace stl = std;
namespace std {
  template<typename T> using host_vector = vector<T>;   
  template<typename T> using device_vector = vector<T>;   
}
#endif

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

#ifdef PAR_MESODYN
    template <typename T>
  inline T val(stl::device_vector<T>& v, int x, int y, int z) {
      return v[x * JX + y * JY + z * JZ];
  }

  template <typename T>
  inline T* val_ptr(stl::device_vector<T>& v, int x, int y, int z) {
      return thrust::raw_pointer_cast(&v[x * JX + y * JY + z * JZ]);
  }
    template <typename T>
  inline T val(stl::host_vector<T>& v, int x, int y, int z) {
      return v[x * JX + y * JY + z * JZ];
  }

  template <typename T>
  inline T* val_ptr(stl::host_vector<T>& v, int x, int y, int z) {
      return &v[x * JX + y * JY + z * JZ];
  }
#endif

      template <typename T>
  inline T val(std::vector<T>& v, int x, int y, int z) {
      return v[x * JX + y * JY + z * JZ];
  }

  template <typename T>
  inline T* val_ptr(std::vector<T>& v, int x, int y, int z) {
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
  Lattice* Lat;

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
  Gaussian_noise(shared_ptr<Boundary1D>, Real, int, Real, Real); // Not seeded (32 bits of randomness)
  Gaussian_noise(shared_ptr<Boundary1D>, Real, int, Real, Real, size_t); // Seeded
  int generate(size_t);
  int add_noise(stl::device_vector<Real>&);
  stl::device_vector<Real> noise;

private:
  seed_seq seed;
  mt19937 prng;
  normal_distribution<Real> dist;
  shared_ptr<Boundary1D> boundary;
};

class Boundary1D : public Lattice_Access {
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
  Component(Lattice*, shared_ptr<Boundary1D>, stl::host_vector<Real>&); //1D
  ~Component();

  stl::device_vector<Real> rho;
  stl::device_vector<Real> alpha;

  Real* rho_ptr;
  Real* alpha_ptr;

  Real rho_at(int, int, int);
  Real alpha_at(int, int, int);
  int update_density(stl::device_vector<Real>&, int = 1);     //Explicit scheme
  int update_density(stl::device_vector<Real>&, stl::device_vector<Real>&, Real ratio, int = 1); //Implicit scheme
  int load_alpha(stl::device_vector<Real>);
  int load_alpha(Real*);
  int load_rho(stl::device_vector<Real>);
  int load_rho(Real*);
  int update_boundaries();
  Real theta();

private:
  shared_ptr<Boundary1D> boundary;
};

class Flux1D : protected Lattice_Access {
public:
  Flux1D(Lattice*, Real, stl::host_vector<int>&, shared_ptr<Component>, shared_ptr<Component>);
  virtual ~Flux1D();

  virtual int langevin_flux();

  Real J_at(int, int, int);
  Real L_at(int, int, int);
  Real mu_at(int, int, int);

  enum error {
    ERROR_SIZE_INCOMPATIBLE,
    ERROR_NOT_IMPLEMENTED,
  };

  stl::device_vector<Real> J_plus;
  stl::device_vector<Real> J_minus;
  stl::device_vector<Real> J;


protected:
  int onsager_coefficient(stl::device_vector<Real>&, stl::device_vector<Real>&);
  int potential_difference(stl::device_vector<Real>&, stl::device_vector<Real>&);
  int langevin_flux(stl::host_vector<int>&, stl::host_vector<int>&, int);
  int mask(stl::host_vector<int>&);
  shared_ptr<Component> A;
  shared_ptr<Component> B;

  stl::device_vector<Real> L;
  stl::device_vector<Real> mu;
  const Real D;
  const int JX;
  stl::host_vector<int> Mask_plus_x;
  stl::host_vector<int> Mask_minus_x;
};

class Flux2D : public Flux1D {
public:
  Flux2D(Lattice*, Real, stl::host_vector<int>&, shared_ptr<Component>, shared_ptr<Component>);
  virtual ~Flux2D();

  virtual int langevin_flux() override;


protected:
  int mask(stl::host_vector<int>&);
  const int JY;
  stl::host_vector<int> Mask_plus_y;
  stl::host_vector<int> Mask_minus_y;

};

class Flux3D : public Flux2D {
public:
  Flux3D(Lattice*, Real, stl::host_vector<int>&, shared_ptr<Component>, shared_ptr<Component>);
  ~Flux3D();

  virtual int langevin_flux() override;

private:
  int mask(stl::host_vector<int>&);

protected:
  const int JZ;
  stl::host_vector<int> Mask_plus_z;
  stl::host_vector<int> Mask_minus_z;
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
  void explicit_start();
  int noise_flux();
  Real* solve_crank_nicolson();
  void load_alpha(Real*, size_t);
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
  stl::device_vector<Real> rho;
  vector<string> tokenize(string, char);
  string read_filename;
  int initial_conditions();
  int init_rho_homogeneous(stl::host_vector<stl::host_vector<Real>>&,stl::host_vector<int>&);
  int norm_density(vector<Real>& rho, Real theta);
  void set_update_lists();
  vector<vector<int>> update_plus;
  vector<vector<int>> update_minus;
  Real boundaryless_volume;

  /* Helper class instances */
  unique_ptr<Gaussian_noise> gaussian;
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
  Mesodyn(vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<State*>, vector<Reaction*>, vector<Molecule*>, vector<System*>, vector<Solve_scf*>, string);
  ~Mesodyn();

  bool mesodyn();

  int norm_theta(vector< shared_ptr<Component> >&);

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
