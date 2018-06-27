#ifndef CLENGxH
#define CLENGxH
#include "input.h"
#include "namics.h"
#include "solve_scf.h"
#include "system.h"
#include "output.h"


//#include <map>  //For generating the output of settings
//#include <cassert>
//random
#include <random>
#include <math.h>       /* exp */
#include <cstdlib>
#include <ctime>
#include <iostream>

using std::setprecision;


class Cleng {
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
    Cleng(vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<Molecule*>, vector<System*>, vector<Solve_scf*>, string);
    ~Cleng();

    int clamp_seg;
    int clp_mol;
    int n_boxes;
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

    vector<int> shift_XYZ;


    bool MonteCarlo();

    void push(string, Real);
    void push(string, int);
    void push(string, bool);
    void push(string, string);
    bool CP(transfer);
    bool MakeShift();
    void WriteOutput(int);
    void PushOutput();
    int GetValue(string, int&, Real&, string&);
    int GetIntRandomValueExclude(int, int, int, bool);
//  int GetIntRandomValue(int, int);
    Real GetRealRandomValue(int, int);

    std::vector<string> KEYS;
    std::vector<string> PARAMETERS;
    std::vector<string> VALUES;
    bool CheckInput(int);
    void PutParameter(string);
    string GetValue(string);
};
#endif
