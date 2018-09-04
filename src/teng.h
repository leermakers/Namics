#ifndef TENGxH
#define TENGxH
#include "input.h"
#include "namics.h"
#include "solve_scf.h"
#include "system.h"
#include "output.h"

#include <random> //For noise in langevinFlux()
//#include <ctime> //To append timestamp to filename.
#include <map>  //For generating the output of settings
#include <ctime>
#include <cassert>


class Teng {
private:
  const string name;
  const vector<Input*> In;
  const vector<Lattice*> Lat;
  const vector<Molecule*> Mol;
  const vector<Segment*> Seg;
  const vector<System*> Sys;
  const vector<Solve_scf*> New;
  const string brand;

  Real seed;
  void prepareOutputFile();


public:
  Teng(vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<Molecule*>, vector<System*>, vector<Solve_scf*>, string);
  ~Teng();

  int tag_seg;
  int tag_mol;
  int n_particles;
  int n_out;
  int sub_box_size;  
  int MCS;
  int n_p;  
  int t; 
  int save_interval;
  string save_filename;
    vector<Output*> Out;

int PX;
int PY;
int PZ;
vector<int> P;
vector<int> Sx;
vector<int> Sy;
vector<int> Sz; 
vector<int> X;
vector<int> Y;
vector<int> Z;
vector<int> X_bm;
vector<int> Y_bm;
vector<int> Z_bm;

  bool MonteCarlo();
  bool ChangeMode();
  bool IsLegal();

  void push(string, Real);
  void push(string, int);
  void push(string, bool);
  void push(string, string);
  bool CP(transfer);
  void WriteOutput(int);
  void PushOutput();
  int GetValue(string, int&, Real&, string&);
  int GetRandom(int);
  Real GetRandom(Real);

  std::vector<string> KEYS;
  std::vector<string> PARAMETERS;
  std::vector<string> VALUES;
  bool CheckInput(int);
  void PutParameter(string);
  string GetValue(string);

};
#endif
