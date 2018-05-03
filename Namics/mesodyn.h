#ifndef MESODYNxH
#define MESODYNxH
#include "input.h"
#include "namics.h"
#include "newton.h"
#include "system.h"
#include <random>

class Mesodyn {
public:
  Mesodyn(vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<Molecule*>, vector<System*>, vector<Newton*>, string);

  ~Mesodyn();
  void AllocateMemory();

  bool success {true};
  string location;

  void gaussianNoise(Real, Real, unsigned long);
  void pepareForCalculations();
  int findComponentNo();
  void setNeighborIndices(vector<int>&, vector<int>&, vector<int>&);
  void abort();
  void langevinFlux(vector<Real>&, vector<Real>&, vector<Real>&, vector<Real>&);
  void updateDensity();
  bool mesodyn();
  int findComponentVectorSize();
  void setComponentStartIndices (vector<int>&);

  vector<Real> noise;
  int componentNo;
  int size;
  vector<Real> J;
  vector<int> xNeighbors;
  vector<int> yNeighbors;
  vector<int> zNeighbors;
  vector<int> component;

  //dummy variables
  //real variables
  Real D {0.5}; //diffusion constant

  //delete when done
  Real dummyMean {1};
  Real dummyStdev {1};
  vector<Real> dummyVector {1,2,3,4,5};
  int latticeX;
  int latticeY;
  int latticeZ;

  string name;
  vector<Input*> In;
  vector<Lattice*> Lat;
  vector<Segment*> Seg;
  vector<Molecule*> Mol;
  vector<System*> Sys;
  vector<Newton*> New;
  string brand;

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
  int timesteps;
  int timebetweensaves;

  std::vector<string> KEYS;
  std::vector<string> PARAMETERS;
  std::vector<string> VALUES;
  bool CheckInput(int);
  void PutParameter(string);
  string GetValue(string);
};
#endif
