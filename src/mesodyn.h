#ifndef MESODYNxH
#define MESODYNxH
#include "input.h"
#include "namics.h"
#include "solve_scf.h"
#include "system.h"
#include "newton.h"   // for...
#include <random>     // noise generation
#include <ctime>      // output labling
#include <cassert>    // making sure all underlying assumptions (e.g. lattice geometry) are satisfied
#include <functional> // boundary conditions
#include <algorithm>  // transform, copy, find functions
#include <limits.h>   // output
#include <unistd.h>   // output

class Newton;

class Gaussian_noise {
  //Makes sure that we keep generating new numbers, instead of the same over and over.

public:
  Gaussian_noise(Real); // Not seeded (32 bits of randomness)
  Gaussian_noise(Real, size_t); // Seeded
  Real noise();

private:
  seed_seq seed;
  mt19937 prng;
  normal_distribution<Real> dist;
};

class Lattice_Access {

public:

  Lattice_Access(Lattice*);
  ~Lattice_Access();

  //TODO: Template these:
  inline Real val(vector<Real>&, int, int, int);
  inline int val(vector<int>&, int, int, int);
  inline Real* valPtr(vector<Real>&, int, int, int);
  inline int index(int x, int y, int z);

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

class Component1D : protected Lattice_Access {
public:

  enum boundary {
    MIRROR,
    PERIODIC,
    BULK,
  };

  enum error {
    ERROR_PERIODIC_BOUNDARY,
    ERROR_SIZE_INCOMPATIBLE,
    ERROR_NOT_IMPLEMENTED,
  };

  Component1D(Lattice*, vector<Real>&, boundary, boundary); //1D
  ~Component1D();

  vector<Real> rho;
  vector<Real> alpha;

  Real rho_at(int, int, int);
  Real alpha_at(int, int, int);
  int update_density(vector<Real>&, int = 1);     //Explicit scheme
  int update_density(vector<Real>&, vector<Real>&, vector<Real>&, int = 1); //Implicit scheme
  int load_alpha(vector<Real>&);
  int load_rho(vector<Real>&);
  int update_boundaries();

private:
  function<void()> bX0;
  function<void()> bXm;

  int set_x_boundaries(boundary, boundary);
  void bX0Mirror(int, int);
  void bXmMirror(int, int, int);
  void bXPeriodic(int, int, int);
  void bX0Bulk(int, int, Real);
  void bXmBulk(int, int, int, Real);
};

class Component2D : public Component1D {
public:
  Component2D(Lattice*, vector<Real>&, boundary, boundary, boundary, boundary); //1D
  ~Component2D();

  int update_boundaries();

private:
  function<void()> bY0;
  function<void()> bYm;

  int set_y_boundaries(boundary, boundary);
  void bY0Mirror(int, int);
  void bYmMirror(int, int, int);
  void bYPeriodic(int, int, int);
  void bY0Bulk(int, int, Real);
  void bYmBulk(int, int, int, Real);
};

class Component3D : public Component2D {
public:
  Component3D(Lattice*, vector<Real>&, boundary, boundary, boundary, boundary, boundary, boundary); //1D
  ~Component3D();

  int update_boundaries();

private:
  function<void()> bZ0;
  function<void()> bZm;

  int set_z_boundaries(boundary, boundary);
  void bZ0Mirror(int, int);
  void bZmMirror(int, int, int);
  void bZPeriodic(int, int, int);
  void bZ0Bulk(int, int, Real);
  void bZmBulk(int, int, int, Real);
};

class Flux1D : protected Lattice_Access {
public:
  Flux1D(Lattice*, Gaussian_noise*, Real, vector<int>&, Component1D*, Component1D*);
  ~Flux1D();

  int langevin_flux();

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

private:
  Component1D* A;
  Component1D* B;


protected:
  int onsager_coefficient(vector<Real>&, vector<Real>&);
  int potential_difference(vector<Real>&, vector<Real>&);
  int langevin_flux(vector<int>&, vector<int>&, int);
  int mask(vector<int>&);

  Gaussian_noise* gaussian;

  vector<Real> L;
  vector<Real> mu;
  const Real D;
  const int JX;
  vector<int> Mask_plus_x;
  vector<int> Mask_minus_x;
  vector<Real> gaussian_noise;
};

class Flux2D : public Flux1D {
public:
  Flux2D(Lattice*, Gaussian_noise*, Real, vector<int>&, Component1D*, Component1D*);
  ~Flux2D();

  int langevin_flux();

private:
  Component2D* A;
  Component2D* B;


protected:
  int mask(vector<int>&);
  const int JY;
  vector<int> Mask_plus_y;
  vector<int> Mask_minus_y;

};

class Flux3D : public Flux2D {
public:
  Flux3D(Lattice*, Gaussian_noise*, Real, vector<int>&, Component1D*, Component1D*);
  ~Flux3D();

  int langevin_flux();

private:
Component3D* A;
Component3D* B;
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
  const vector<System*> Sys;
  const vector<Solve_scf*> New;
  const string brand;

  /* Read from file */
  Real D; // diffusionconstant
  Real mean; // mean of gaussian noise (should be 0)
  Real stdev; // stdev of gaussian noise (should be 1*D)
  Real seed;  // seed of gaussian noise
  int timesteps; // length of the time evolution
  int timebetweensaves; // how many timesteps before mesodyn writes the current variables to file
  const int componentNo; // number of components in the system, read from SysMonMolList

  /* Initialization*/
  int initial_conditions();
  vector<Real>&  flux_callback(int);
  int init_rho(vector< vector<Real> >&, vector<int>&);

  /* Helper class instances */
  vector<Component1D*> component;
  vector<Component1D*> solver_component;
  vector<Flux1D*> flux;
  vector<Flux1D*> solver_flux;

  /* Mesodyn specific output */
  ofstream mesFile;
  void prepareOutputFile();
  void writeRho(int);

  /* Mathematics */
  int factorial (int);
  int combinations (int, int);


public:
  Mesodyn(vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<Molecule*>, vector<System*>, vector<Solve_scf*>, string);
  ~Mesodyn();

  bool mesodyn();

  void gaussianNoise(Real, Real, unsigned int);
  Gaussian_noise* gaussian_noise;


  /* Inputs / output class interface functions */
  vector<string> ints;
  vector<string> Reals;
  vector<string> bools;
  vector<string> strings;
  vector<Real> Reals_value;
  vector<int> ints_value;
  vector<bool> bools_value;
  vector<string> strings_value;
  void push(string, Real);
  void push(string, int);
  void push(string, bool);
  void push(string, string);
  void PushOutput();
  int GetValue(string, int&, Real&, string&);

  std::vector<string> KEYS;
  std::vector<string> PARAMETERS;
  std::vector<string> VALUES;
  bool CheckInput(int);
  void PutParameter(string);
  string GetValue(string);
};
#endif
