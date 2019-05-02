#include "cleng.h"
#include "cleng_tools.h"
//#include <unistd.h>
#include "iterator/EnumerateIterator.h"

using namespace std;

Cleng::Cleng(
        vector<Input *> In_,
        vector<Lattice *> Lat_,
        vector<Segment *> Seg_,
        vector<State *> Sta_,
        vector<Reaction *> Rea_,
        vector<Molecule *> Mol_,
        vector<System *> Sys_,
        vector<Solve_scf *> New_,
        string name_
) : name(std::move(name_)),
    In(std::move(In_)),
    Lat(std::move(Lat_)),
    Seg(std::move(Seg_)),
    Sta(std::move(Sta_)),
    Rea(std::move(Rea_)),
    Mol(std::move(Mol_)),
    Sys(std::move(Sys_)),
    New(std::move(New_)) {

    if (debug) cout << "Cleng initialized" << endl;
    KEYS.emplace_back("MCS");
    KEYS.emplace_back("delta_step");
    KEYS.emplace_back("delta_save");
    KEYS.emplace_back("save_filename");
    KEYS.emplace_back("seed");
    KEYS.emplace_back("checkpoint_save");
    KEYS.emplace_back("checkpoint_load");
    KEYS.emplace_back("cleng_pos");
    KEYS.emplace_back("simultaneous");
    KEYS.emplace_back("movement_alone");
    KEYS.emplace_back("F_dependency");

    // Debug.log
    //out.open("debug.out", ios_base::out);
}

Cleng::~Cleng() {
    delete[] xs;
    delete[] ys;
    delete[] zs;
    // Debug.log closing
    //out.close();
}

string Cleng::GetValue(string parameter) {
    int i = 0;
    int length = (int) PARAMETERS.size();
    while (i < length) {
        if (parameter == PARAMETERS[i]) {
            return VALUES[i];
        }
        i++;
    }
    return "";
}

bool Cleng::CheckInput(int start, bool save_vector) {
    if (debug) cout << "CheckInput in Cleng" << endl;
    bool success;

    success = In[0]->CheckParameters("cleng", name, start, KEYS, PARAMETERS, VALUES);
    if (success) {

        // MCS
        if (!GetValue("MCS").empty()) {
            success = In[0]->Get_int(GetValue("MCS"), MCS, 1, 10000, "The number of Monte Carlo steps should be between 1 and 10000");
            if (!success) { cout << "MCS will be equal to 5" << endl; MCS = 1; }
        } else MCS = 5;
        if (debug) cout << "MCS is " << MCS << endl;

        // seed
        if (!GetValue("seed").empty()) {
            success = In[0]->Get_int(GetValue("seed"), pseed, 1, 1000, "The seed should be between 1 and 1000");
            if (!success) { cout << "The seed will be equal 1" << endl; pseed = 1; }
        } else pseed = 0;
        rand = pseed==0 ? Random() : Random(pseed);
        if (debug) cout << "seed is " << pseed << endl;

        // delta_step
        if (!GetValue("delta_step").empty()) {
            success = In[0]->Get_int(GetValue("delta_step"), delta_step, 1, 5, "The number of delta_step should be between 1 and 5");
            if (!success) { cout << "The delta_step will be equal 1" << endl; delta_step = 1;}
        } else delta_step = 1;
        if (debug) cout << "delta_step is " << delta_step << endl;

        // delta_save
        if (!GetValue("delta_save").empty()) {
            success = In[0]->Get_int(GetValue("delta_save"), delta_save, 1, MCS, "The delta_save interval should be between 1 and " + to_string(MCS));
        } else delta_save = 1;
        if (debug) cout << "delta_save_interval " << delta_save << endl;

        // Cleng molecules
        if (Sys[0]->SysClampList.empty()) {
            cout << "Cleng needs to have clamped molecules in the system" << endl;
            success = false;
        } else {
            clamp_seg = Sys[0]->SysClampList[0];
            if (Sys[0]->SysClampList.size() > 1) {
                success = false;
                cout << "Currently the clamping is limited to one molecule per system. " << endl;
            }
        }

        // checkpoint save
        if (!GetValue("checkpoint_save").empty()) {
            checkpoint_save = In[0]->Get_bool(GetValue("checkpoint_save"), false);
        } else checkpoint_save = false;
        if (debug) cout << "checkpoint_save " << checkpoint_save << endl;

        // checkpoint load
        if (!GetValue("checkpoint_load").empty()) {
            checkpoint_load = In[0]->Get_bool(GetValue("checkpoint_load"), false);
        } else checkpoint_load = false;
        if (debug) cout << "checkpoint_load " << checkpoint_load << endl;

        // saving cleng position of nodes
        if (!GetValue("cleng_pos").empty()) {
            cleng_pos = In[0]->Get_bool(GetValue("cleng_pos"), false);
        } else cleng_pos = false;
        if (debug) cout << "cleng_pos " << cleng_pos << endl;

        // simultaneous
        if (!GetValue("simultaneous").empty()) simultaneous = In[0]->Get_bool(GetValue("simultaneous"), false);
        else simultaneous = false;
        if (debug) cout << "simultaneous move " << simultaneous << endl;

        // movement_alone
        if (!GetValue("movement_alone").empty()) {
            success = In[0]->Get_int(GetValue("movement_alone"), axis, 1, 3, "The number of delta_step should be between 1 and 3");
            if (!success) {
                cout << "Sorry, you provide incorrect axis number. The movement_alone will be disable" << endl;
                axis = -1;
                success = true;
            }
        } else axis = -1;
        if (debug) cout << "movement_alone axis " << axis << endl;

        // F_dependent
        if (!GetValue("F_dependency").empty()) F_dependency = In[0]->Get_bool(GetValue("F_dependency"), false);
        else F_dependency = false;
        if (debug) cout << "F_dependent move " << F_dependency << endl;

        // ???
        if (success) { n_boxes = Seg[clamp_seg]->n_box; sub_box_size = Seg[clamp_seg]->mx; }
        clp_mol = -1;
        int length = (int) In[0]->MolList.size();
        for (int i = 0; i < length; i++) if (Mol[i]->freedom == "clamped") clp_mol = i;
    }

    if (success) {
        n_out = (int) In[0]->OutputList.size();
        if (n_out == 0) cout << "Warning: no output defined!" << endl;

        for (int i = 0; i < n_out; i++) {
            Out.push_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[i], i, n_out));
            if (!Out[i]->CheckInput(start)) {
                cout << "input_error in output " << endl;
                success = false;
            }
        }
        cout << cmajor_version << "." << cminor_version << "." << cpatch << " v." << cversion << endl;
        success = MonteCarlo(save_vector);
    }
    return success;
}

bool Cleng::CP(transfer tofrom) {
    if (debug) cout << "CP in Cleng" << endl;

    bool success = true;
    Point box{Lat[0]->MX, Lat[0]->MY, Lat[0]->MZ};

    int JX = Lat[0]->JX;
    int JY = Lat[0]->JY;
    int M  = Lat[0]->M;

    Segment *clamped = Seg[clamp_seg];
    map<int, Point> system_points;
    switch (tofrom) {
        case to_cleng:
            simpleNodeList.clear();
            for (int i = 0; i < n_boxes; i++) {
                auto first_node  = fromSystemToNode(clamped->px1[i], clamped->py1[i], clamped->pz1[i], 2 * i, box);
                auto second_node = fromSystemToNode(clamped->px2[i], clamped->py2[i], clamped->pz2[i], 2 * i + 1, box);
                first_node ->set_cnode(second_node);
                second_node->set_cnode(first_node);
                //
                simpleNodeList.push_back(first_node);
                simpleNodeList.push_back(second_node);
            }; nodes = createNodes(simpleNodeList);
            break;

        case to_segment:
            Zero(clamped->H_MASK, M);

            for (auto &&SN : Enumerate(simpleNodeList)) {
                size_t index = SN.first;    //
                if (index % 2 == 1) continue;
                size_t n_box = index / 2;   // n_boxes
                Point p1 = SN.second->point();                // all time in box [0:box_l]
                Point p2 = SN.second->get_cnode()->point();   // all time in box [0:box_l]

                Point p1sys  = SN.second->get_system_point();
                Point p2sys  = SN.second->get_cnode()->get_system_point();

//                cout << " n_box: " << n_box << endl;
//                cout << " index: " << index << " v: " << p1.x << " " << "("<<p1sys.x<<")" << " " << p1.y << " "  << "("<<p1sys.y<<")" << " " << p1.z << " "  << "("<<p1sys.z<<")" << endl;
//                cout << " index: " << index << " v: " << p2.x << " " << "("<<p2sys.x<<")" << " " << p2.y << " "  << "("<<p2sys.y<<")" << " " << p2.z << " "  << "("<<p2sys.z<<")" << endl;

                if ((p1.x - p2.x) != (p1sys.x - p2sys.x)) {if (p1.x < p2.x) p1.x += box.x; else p2.x += box.x;}
                if ((p1.y - p2.y) != (p1sys.y - p2sys.y)) {if (p1.y < p2.y) p1.y += box.y; else p2.y += box.y;}
                if ((p1.z - p2.z) != (p1sys.z - p2sys.z)) {if (p1.z < p2.z) p1.z += box.z; else p2.z += box.z;}

                clamped->px1[n_box] = p1.x;
                clamped->py1[n_box] = p1.y;
                clamped->pz1[n_box] = p1.z;

                clamped->px2[n_box] = p2.x;
                clamped->py2[n_box] = p2.y;
                clamped->pz2[n_box] = p2.z;
// box
                clamped->bx[n_box] = (p2.x + p1.x - sub_box_size) / 2;
                clamped->by[n_box] = (p2.y + p1.y - sub_box_size) / 2;
                clamped->bz[n_box] = (p2.z + p1.z - sub_box_size) / 2;

                if (clamped->bx[n_box] < 1) { clamped->bx[n_box] += box.x; clamped->px1[n_box] += box.x; clamped->px2[n_box] += box.x;}
                if (clamped->by[n_box] < 1) { clamped->by[n_box] += box.y; clamped->py1[n_box] += box.y; clamped->py2[n_box] += box.y;}
                if (clamped->bz[n_box] < 1) { clamped->bz[n_box] += box.z; clamped->pz1[n_box] += box.z; clamped->pz2[n_box] += box.z;}
// clearing
                auto hp1x = ((clamped->px1[n_box] - 1) % box.x + 1) * JX;
                auto hp1y = ((clamped->py1[n_box] - 1) % box.y + 1) * JY;
                auto hp1z =  (clamped->pz1[n_box] - 1) % box.z + 1;

                auto hp2x = ((clamped->px2[n_box] - 1) % box.x + 1) * JX;
                auto hp2y = ((clamped->py2[n_box] - 1) % box.y + 1) * JY;
                auto hp2z = ( clamped->pz2[n_box] - 1) % box.z + 1;

                clamped->H_MASK[hp1x + hp1y + hp1z] = 1;
                clamped->H_MASK[hp2x + hp2y + hp2z] = 1;
            }
            break;

        default:
            success = false;
            cout << "error in transfer" << endl;
            break;
    }
    return success;
}

bool Cleng::Checks() {
    bool result;
    bool not_collapsing;
    bool in_range;
    bool in_subbox_range;

    // check subbox_ranges
    in_subbox_range = InSubBoxRange();

    // check distances between all nodes => preventing collapsing
    not_collapsing = NotCollapsing();

    // check distance between all nodes and constrains (walls)
    // BC.x, BC.y, BC.z = mirror
    if (BC.x and BC.y and BC.z) in_range = InRange();
    // BC.x and/or BC.y and/or BC.z != mirror
    else in_range = true;

    result = not_collapsing and in_range and in_subbox_range;
    return result;
}

bool Cleng::MakeMove(bool back) {

    if (debug) cout << "MakeMove in Cleng" << endl;
    bool success = true;

    if (back) {
        cout << "MakeMove back" << endl;
        Point _clamped_move = clamped_move.negate();
        if (!simultaneous) nodes[id_node_for_move]->shift(_clamped_move);
        else for (auto &node : nodes) node->shift(_clamped_move);

    } else {
        clamped_move = prepareMove();
        if (!simultaneous) {
            if (!F_dependency) {
                id_node_for_move = rand.getInt(0, (int) nodes.size() - 1);
                cout << "Prepared id: " << id_node_for_move << " clamped_move: " << clamped_move.to_string() << endl;
                while (!Checks()) {
                    cout << "Prepared MC step for a node does not pass checks. Rejected." << endl;
                    clamped_move = {0, 0, 0};
                    rejected++;
                    success = false;
                }
                nodes[id_node_for_move]->shift(clamped_move);
                cout << "Moved: \n" << "node: " << id_node_for_move << ", " << "MC step: " << clamped_move.to_string()
                     << endl;
            } else {
                id_node_for_move = 1;
                cout << "Prepared id: " << id_node_for_move << " clamped_move: " << clamped_move.to_string() << endl;
                if (!Checks()) {
                    cout << "Prepared MC step for a node does not pass checks. Rejected." << endl;
                    clamped_move = {0, 0, 0};
                    rejected++;
                    success = false;
                }
                nodes[id_node_for_move]->shift(clamped_move);
                cout << "Moved: \n" << "node: " << id_node_for_move << ", " << "MC step: " << clamped_move.to_string()
                     << endl;
            }
        } else {
            for (auto &node : nodes) node->shift(clamped_move);
            cout << "Moved: \n" << "*All* " << "MC step: " << clamped_move.to_string() << endl;
        }
    }
    return success;
}

bool Cleng::MonteCarlo(bool save_vector) {
    if (debug) cout << "Monte Carlo in Cleng" << endl;
    bool success = true;

    MCS_checkpoint = 0;
    Checkpoint checkpoint;
    if (checkpoint_load) {
        Point box{Lat[0]->MX, Lat[0]->MY, Lat[0]->MZ};
        CP(to_cleng);
        nodes = checkpoint.loadCheckpoint(nodes, box);
        CP(to_segment);
        if (getLastMCS() != 0) MCS_checkpoint = getLastMCS() + 1;
    }

// Analysis MC
    accepted = 0.0;
    rejected = 0.0;
    MC_attempt = 0;
    make_BC();


// init system outlook
    New[0]->Solve(true);
    free_energy_current = Sys[0]->GetFreeEnergy() - GetN_times_mu();

// init save
    WriteOutput(MC_attempt + MCS_checkpoint);

    cout << "Initialization done.\n" << endl;
    CP(to_cleng);
    cout << "Here we go..." << endl;
    for (MC_attempt = 1; MC_attempt <= MCS; MC_attempt++) { // loop for trials
        bool success_;

        success_ = MakeMove(false);
        CP(to_segment);

        cout << endl;
        cout << "[Cleng] System for calculation: " << endl;
        for (auto &&n : nodes) cout << n->to_string() << endl;
        cout << endl;

        if (success_) {
            New[0]->Solve(true);
            free_energy_trial = Sys[0]->GetFreeEnergy() - GetN_times_mu();
            if (save_vector) test_vector.push_back(free_energy_trial);

            cout << "Free Energy (c): " << free_energy_current    << endl;
            cout << "            (t): " << free_energy_trial << endl;
//            cout << "trial - current = " << std::to_string(free_energy_trial - free_energy_current) << endl;
            if (is_ieee754_nan( free_energy_trial )) {cout << " Sorry, Free Energy is NaN. Termination..." << endl; break;}

            if (free_energy_trial - free_energy_current <= 0.0) {
                cout << "Accepted" << endl;
                free_energy_current = free_energy_trial;
                accepted++;
            } else {
                Real acceptance = rand.getReal(0, 1);

                if (acceptance < exp(-1.0 * (free_energy_trial - free_energy_current))) {
                    cout << "Accepted with probability" << endl;
                    free_energy_current = free_energy_trial;
                    accepted++;
                } else {
                    cout << "Rejected" << endl;
                    MakeMove(true);
                    CP(to_segment);
                    rejected++;
                }
            }
        }
        cout << "Monte Carlo attempt: " << MC_attempt << endl;
        cout << "Accepted: # " << accepted << " | " << 100 * (accepted / MC_attempt) << "%" << endl;
        cout << "Rejected: # " << rejected << " | " << 100 * (rejected / MC_attempt) << "%" << endl;

        if ((MC_attempt % delta_save) == 0) WriteOutput(MC_attempt + MCS_checkpoint);
        if (checkpoint_save) checkpoint.updateCheckpoint(simpleNodeList);
    }

    cout << endl;
    cout << "Finally:" << endl;
    cout << "Monte Carlo attempt: " << MC_attempt-1 << endl;
    cout << "Accepted: # " << accepted << " | " << 100 * (accepted / (MC_attempt-1)) << "%" << endl;
    cout << "Rejected: # " << rejected << " | " << 100 * (rejected / (MC_attempt-1)) << "%" << endl;

    if (checkpoint_save) { cout << "Saving checkpoint..." << endl; checkpoint.saveCheckpoint(simpleNodeList);}
    cout << "Have a fun. " << endl;
    return success;
}
