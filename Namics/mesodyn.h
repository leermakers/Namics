#ifndef MESODYNxH
#define MESODYNxH
#include "input.h"
#include "namics.h"
#include "newton.h"
#include "system.h"
#include <random> //For noise in langevinFlux()
//#include <ctime> //To append timestamp to filename.
#include <map>  //For generating the output of settings
#include <ctime>
#include <cassert>


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
  const int JZ;
  const int JY;
  const int JX;
  const int MZ;
  const int MY;
  const int MX;
  const int M;
  const int componentNo;
  const int dimensions;
  vector<Real> J;
  vector<Real> L;
  vector<Real> rho;
  vector<Real*> ptrComponentStart;    //densities used in langevinFlux
  vector<Real> U;

  ofstream mesFile;
  void prepareOutputFile();
  void writeRho(int);
  map<string,Real> settings;
  void boundaryConditions();

public:
  Mesodyn(vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<Molecule*>, vector<System*>, vector<Newton*>, string);
  ~Mesodyn();

  void gaussianNoise(Real, Real, unsigned int);
  void quit();
  void langevinFlux();
  bool mesodyn();
  void initRho();
  int factorial (int);
  void onsagerCoefficient();
  void potentialDifference();
  int combinations (int, int);
  inline Real val(vector<Real>&, int, int, int, int);
  inline Real* valPtr(vector<Real>&, int, int, int, int);
  int findDimensions();

  vector<Real> noise;

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
