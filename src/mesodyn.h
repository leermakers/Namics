#ifndef MESODYNxH
#define MESODYNxH
#include "input.h"
#include "namics.h"
#include "solve_scf.h"
#include "system.h"
#include "newton.h"
#include <random>
#include <ctime>
#include <cassert>
#include <functional>
#include <algorithm>
//for output:
#include <limits.h>
#include <unistd.h>

class Newton;

class Access {
public:

  Access(Lattice*);
  ~Access();

  inline Real val(vector<Real>&, int, int, int);
  inline Real* valPtr(vector<Real>&, int, int, int);
  inline int xyz(int, int, int);

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

class Component1D : protected Access {
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
  int update_density(vector<Real>&, int = 1.0);     //Explicit scheme
  int update_density(vector<Real>&, vector<Real>&); //Implicit scheme
  int load_alpha(Real*, int);  // Update when xx becomes a vector
  int load_rho(Real*, int); // Update when xx becomes a vector
  int update_boundaries();

private:
  function<void()> bX0;
  function<void()> bXm;

  int set_x_boundaries(boundary, boundary);
  void bX0Mirror(int, int);
  void bXmMirror(int, int, int);
  void bXPeriodic(int, int, int);
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
};

class Flux1D : private Access {
public:
  Flux1D(Lattice*, Real, Component1D*, Component1D*);
  ~Flux1D();

  int langevin_flux();

  Real J_at(int, int, int);
  Real L_at(int, int, int);
  Real mu_at(int, int, int);

  enum error {
    ERROR_SIZE_INCOMPATIBLE,
    ERROR_NOT_IMPLEMENTED,
  };

  vector<Real> J;

private:
  Component1D* A;
  Component1D* B;

protected:
  int onsager_coefficient(vector<Real>&, vector<Real>&);
  int potential_difference(vector<Real>&, vector<Real>&);
  int langevin_flux(int);

  vector<Real> L;
  vector<Real> mu;
  const Real D;
  const int JX;
};

class Flux2D : public Flux1D {
public:
  Flux2D(Lattice*, Real, Component1D*, Component1D*);
  ~Flux2D();

  int langevin_flux();

private:
Component2D* A;
Component2D* B;


protected:
  const int JY;

};

class Flux3D : public Flux2D {
public:
  Flux3D(Lattice*, Real, Component1D*, Component1D*);
  ~Flux3D();

  int langevin_flux();

private:
Component3D* A;
Component3D* B;

protected:
  const int JZ;
};


class Mesodyn : private Access {

private:
  const string name;
  const vector<Input*> In;
  const vector<Lattice*> Lat;
  const vector<Molecule*> Mol;
  const vector<Segment*> Seg;
  const vector<System*> Sys;
  const vector<Solve_scf*> New;
  const string brand;
  Real D; //diffusionconstant
  Real mean;
  Real stdev;
  Real seed;
  int timesteps;
  int timebetweensaves;
  const int componentNo;
  vector<Real> noise;

  vector<Component1D*> component;
  vector<Flux1D*> flux;

  ofstream mesFile;
  void prepareOutputFile();
  void writeRho(int);
  int initial_conditions();


public:
  Mesodyn(vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<Molecule*>, vector<System*>, vector<Solve_scf*>, string);
  ~Mesodyn();

  void gaussianNoise(Real, Real, unsigned int);
  bool mesodyn();
  int factorial (int);
  int combinations (int, int);
  int findDimensions();
  vector<Real> alpha;

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
