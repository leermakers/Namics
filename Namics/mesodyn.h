#ifndef MESODYNxH
#define MESODYNxH
#include "input.h"
#include "namics.h"
#include "newton.h"
#include "system.h"
#include <random>

class Mesodyn {
private:
  const string name;
  const vector<Input*> In;
  const vector<Lattice*> Lat;
  const vector<Molecule*> Mol;
  const vector<Segment*> Seg;
  const vector<System*> Sys;
  const vector<Newton*> New;
  const string brand;
  Real D; //diffusionconstant
  void updateDensity();
  Real mean;
  Real stdev;
  Real seed;
  int timesteps;
  int timebetweensaves;
  const int zNeighbor;
  const int yNeighbor;
  const int xNeighbor;
  const int cNeighbor;
  const int componentNo;
  const int size;
  Real initRho;
  const int dimensions;
  bool success {true};

public:
  Mesodyn(vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<Molecule*>, vector<System*>, vector<Newton*>, string);
  ~Mesodyn();
  void AllocateMemory();

  void gaussianNoise(Real, Real, unsigned long);;
  void abort();
  void langevinFlux();
  bool mesodyn();
  void fillRho(Real);
  int factorial (int);
  void onsagerCoefficient();
  int combinations (int, int);
  inline Real at(int, int, int, int);

  vector<Real*> ptrComponentStart;    //densities used in langevinFlux
  vector<Real> noise;
  vector<Real> J;
  vector<Real> rho;
  vector<Real> L;

  vector<Real> dummyVector {1,2,3,4,5};

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
  Real* GetPointer(string);
  int GetValue(string, int&, Real&, string&);

  std::vector<string> KEYS;
  std::vector<string> PARAMETERS;
  std::vector<string> VALUES;
  bool CheckInput(int);
  void PutParameter(string);
  string GetValue(string);
};
#endif
