#ifndef CLENGxH
#define CLENGxH
#include "../input.h"
#include "nodes/simple_node.h"
#include "../namics.h"
#include "../solve_scf.h"
#include "../system.h"
#include "../output.h"
#include <random>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
#include <ostream>

using std::setprecision;


class Cleng {
private:
    const string name;
    const vector<Input*> In;
    const vector<Lattice*> Lat;
    const vector<Segment*> Seg;
    const vector<State*> Sta;
    const vector<Reaction*> Rea;    
    const vector<Molecule*> Mol;
    const vector<System*> Sys;
    const vector<Solve_scf*> New;
    const string brand;

    Real seed;
    void prepareOutputFile();
    void fillXYZ();

public:
    Cleng(vector<Input*>, vector<Lattice*>, vector<Segment*>, vector<State*>, vector<Reaction*>, vector<Molecule*>, vector<System*>, vector<Solve_scf*>, string);
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
    Point BC;

    vector<Output*> Out;


    vector<int> P;
    vector<SimpleNode> simpleNodeList;
    vector<std::shared_ptr<Node>> nodes;
    // temporary arrays for keep nodes coordinates for output
    int* xs = nullptr;
    int* ys = nullptr;
    int* zs = nullptr;

    ofstream out;
    Point shift;
    int id_node_for_move;
    Real free_energy_current;
    Real free_energy_trial;


    bool MonteCarlo();

    void push(string, Real);
    void push(string, int);
    void push(string, bool);
    void push(string, string);
    bool CP(transfer);
    bool MakeShift(bool back);
    bool Checks();
    bool InSubBoxRange();
    bool NotCollapsing();
    bool InRange();
    void WriteOutput(int, Real);
    void PushOutput(int, Real);
    void make_BC();
    int GetValue(string, int&, Real&, string&);
    int GetIntRandomValueExclude(int, int, int, bool);
    Real GetRealRandomValue(int, int);
    Real GetN_times_mu();
    Point PrepareStep();

    std::vector<string> KEYS;
    std::vector<string> PARAMETERS;
    std::vector<string> VALUES;
    bool CheckInput(int);
    void PutParameter(string);
    string GetValue(string);

};
#endif
