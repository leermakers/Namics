#ifndef CLENGxH
#define CLENGxH

#include "../input.h"
#include "../namics.h"
#include "../solve_scf.h"
#include "../system.h"
#include "../output.h"
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
#include <ostream>
#include <map>
#include <cassert>
#include <fstream>
#include <utility>
// features
#include "nodes/simple_node.h"
#include "nodes/monolit.h"
#include "nodes/point.h"
#include "checkpoint/checkpoint.h"
#include "random/random.h"
// tests

using std::setprecision;

class Cleng {
private:
    const string name;
    const vector<Input *> In;
    const vector<Lattice *> Lat;
    const vector<Segment *> Seg;
    const vector<State *> Sta;
    const vector<Reaction *> Rea;
    const vector<Molecule *> Mol;
    const vector<System *> Sys;
    const vector<Solve_scf *> New;
    const string brand;

    int cmajor_version = 1;
    int cminor_version = 1;
    int cpatch = 4;
    int cversion = 99;
    int pseed;

    void fillXYZ();

public:
    Cleng(
            vector<Input *>,
            vector<Lattice *>,
            vector<Segment *>,
            vector<State *>,
            vector<Reaction *>,
            vector<Molecule *>,
            vector<System *>,
            vector<Solve_scf *>,
            string
    );

    ~Cleng();

    vector<string> KEYS;
    vector<string> PARAMETERS;
    vector<string> VALUES;

    int clamp_seg;
    int clp_mol;
    int n_boxes;
    int n_out;
    int sub_box_size;
    int MCS;
    int delta_step;
    int t;
    int delta_save;
    bool checkpoint_save;
    bool checkpoint_load;
    bool simultaneous;
    int  axis;
    string sign_move;
    bool cleng_pos;
    vector<int> ids_node4move;
    int user_node_id;
    int MCS_checkpoint;
    string save_filename;
    Point BC;
    Random rand;

    vector<Output *> Out;
    vector<int> P;
    vector<shared_ptr<SimpleNode>> simpleNodeList;
    vector<std::shared_ptr<Node>> nodes;

    // temporary arrays for keep nodes coordinates for output
    int *xs = nullptr;
    int *ys = nullptr;
    int *zs = nullptr;

    ofstream out;
    Point clamped_move;
    int id_node_for_move;
    Real free_energy_current;
    Real free_energy_trial;

    Real accepted;
    Real rejected;
    int MC_attempt;

    vector<Real> test_vector;

    bool MonteCarlo(bool save_vector);

    bool CP(transfer);

    bool MakeMove(bool back);

    bool Checks();

    bool InSubBoxRange();

    bool NotCollapsing();

    bool InRange();

    void PushOutput(int);

    void WriteOutput(int);

    void WriteClampedNodeDistance(int);

    void make_BC();

    int getLastMCS();

    Real GetN_times_mu();

    Point prepareMove();

    bool CheckInput(int start, bool save_vector);

    string GetValue(string);

};

#endif
