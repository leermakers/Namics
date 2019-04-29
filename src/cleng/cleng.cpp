#include "cleng.h"
#include "check_is_nan.h"
//#include <unistd.h>

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
            success = In[0]->Get_int(GetValue("delta_save"), delta_save, 1, MCS, "The delta_save interval should be between 1 and " + std::to_string(MCS));
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

shared_ptr<SimpleNode> fromSystemToNode(int x, int y, int z, int id, const Point &box) {
    return make_shared<SimpleNode>(Point(x, y, z), id, box);
}

vector<shared_ptr<Node>> createNodes(const vector<shared_ptr<SimpleNode>> &simple_nodes) {
    vector<shared_ptr<Node>> result;
    map<SimpleNode, vector<shared_ptr<SimpleNode>>> m;
    for (auto &&n  : simple_nodes) {
        m[*n].push_back(n);
    }
    for (auto &&entry : m) {
        if (entry.second.size() == 1) {
            result.push_back(entry.second[0]);
        } else {
            result.push_back(make_shared<Monolit>(entry.second));
        }
    }

    return result;
}

bool Cleng::CP(transfer tofrom) {
    if (debug) cout << "CP in Cleng" << endl;

    bool success = true;
    Point box{Lat[0]->MX, Lat[0]->MY, Lat[0]->MZ};

    int JX = Lat[0]->JX;
    int JY = Lat[0]->JY;
    int M = Lat[0]->M;

    Segment *clamped = Seg[clamp_seg];
    map<int, Point> system_points;
    switch (tofrom) {
        case to_cleng:
            simpleNodeList.clear();
            for (int i = 0; i < n_boxes; i++) {
                auto first_node = fromSystemToNode(clamped->px1[i], clamped->py1[i], clamped->pz1[i], 2 * i, box);
                auto second_node = fromSystemToNode(clamped->px2[i], clamped->py2[i], clamped->pz2[i], 2 * i + 1, box);
                first_node->set_cnode(second_node);
                second_node->set_cnode(first_node);
                //
                simpleNodeList.push_back(first_node);
                simpleNodeList.push_back(second_node);
            }

            nodes = createNodes(simpleNodeList);
//            for (auto &&n : nodes) {
//                cout << n->to_string() << endl;
//            }
//            assert(nodes.size() != 3);
            break;

        case to_segment:
            Zero(clamped->H_MASK, M);

            for (auto &&n : nodes) n->pushSystemPoints(system_points);

            for (int index=0; index < (int) system_points.size(); index=index+2) {
                auto first_node = system_points[index];
                auto second_node = system_points[index+1];

                if (!first_node.point_in_range(box) and !second_node.point_in_range(box)) {
                    cout << "1* box" << endl;
                    system_points[index]   = first_node % box;
                    system_points[index+1] = second_node % box;
                }
                if (!first_node.point_in_range(box+box) and !second_node.point_in_range(box+box)) {
                    cout << "2* box" << endl;
                    system_points[index]   = first_node % box;
                    system_points[index+1] = second_node % box;
                }
            }

            for (auto &&entry : system_points) {
                int i = entry.first / 2;
                Point &p = entry.second;
                if (entry.first % 2 == 0) {
                    clamped->px1[i] = p.x;
                    clamped->py1[i] = p.y;
                    clamped->pz1[i] = p.z;
                } else {
                    clamped->px2[i] = p.x;
                    clamped->py2[i] = p.y;
                    clamped->pz2[i] = p.z;
                }
            }

            for (int i = 0; i < n_boxes; i++) {
                clamped->bx[i] = (clamped->px2[i] + clamped->px1[i] - sub_box_size) / 2;
                clamped->by[i] = (clamped->py2[i] + clamped->py1[i] - sub_box_size) / 2;
                clamped->bz[i] = (clamped->pz2[i] + clamped->pz1[i] - sub_box_size) / 2;

                if (clamped->bx[i] < 1) {
                    clamped->bx[i] += box.x;
                    clamped->px1[i] += box.x;
                    clamped->px2[i] += box.x;
                }
                if (clamped->by[i] < 1) {
                    clamped->by[i] += box.y;
                    clamped->py1[i] += box.y;
                    clamped->py2[i] += box.y;
                }
                if (clamped->bz[i] < 1) {
                    clamped->bz[i] += box.z;
                    clamped->pz1[i] += box.z;
                    clamped->pz2[i] += box.z;
                }

                clamped->H_MASK[((clamped->px1[i] - 1) % box.x + 1) * JX + ((clamped->py1[i] - 1) % box.y + 1) * JY +
                                (clamped->pz1[i] - 1) % box.z + 1] = 1;
                clamped->H_MASK[((clamped->px2[i] - 1) % box.x + 1) * JX + ((clamped->py2[i] - 1) % box.y + 1) * JY +
                                (clamped->pz2[i] - 1) % box.z + 1] = 1;

            }
            break;

        default:
            success = false;
            cout << "error in transfer" << endl;
            break;
    }
    return success;
}

Real Cleng::GetN_times_mu() {
    int n_mol = (int) In[0]->MolList.size();
    Real n_times_mu = 0;
    for (int i = 0; i < n_mol; i++) {
        Real Mu = Mol[i]->Mu;
        Real n = Mol[i]->n;
        if (Mol[i]->IsClamped()) n = Mol[i]->n_box;
        n_times_mu += n * Mu;
    }
    return n_times_mu;
}

bool Cleng::InSubBoxRange() {
//    cout << "InSubBoxRange ... " << endl;
    int not_in_subbox_range = 0;
    Point sub_box = {sub_box_size-2, sub_box_size-2, sub_box_size-2};  // change it in future
//    cout << "sub_box: " << sub_box.to_string() << endl;
    if (!nodes[id_node_for_move]->inSubBoxRange(sub_box, clamped_move)) not_in_subbox_range++;

    if (not_in_subbox_range != 0) {
        cout << "There are nodes not in sub-box range!" << endl;
        return false;
    }
    return true;
}

bool Cleng::NotCollapsing() {
    bool not_collapsing = false;
    Point MP(nodes[id_node_for_move]->point());
    Point MPs(nodes[id_node_for_move]->point() + clamped_move);

//    cout << "MP: "  << MP.to_string() << endl;
//    cout << "MPs: " << MPs.to_string() << endl;

    double min_dist = 2; // minimal distance between nodes
    int i = 0;
    for (const auto &n : nodes) {
        if (MP != n->point()) {
            Real distance = MPs.distance(n->point());
            if (distance < min_dist) {
                cout << "Nodes are too close to each other." << endl;
                i++;
            }
        }
    }
    if (i == 0) {
        not_collapsing = true;
    }
    return not_collapsing;
}

bool Cleng::InRange() {
    bool in_range = false;
    Point box = {Lat[0]->MX, Lat[0]->MY, Lat[0]->MZ};
    Point down_boundary = {1, 1, 1};
    Point up_boundary = box - down_boundary;
    Point MPs(nodes[id_node_for_move]->point() + clamped_move);

//    cout << "MPs" << MPs.to_string() << endl;
//    cout << "down bound" << down_boundary.to_string() << endl;
//    cout << "up boundary bound" << up_boundary.to_string() << endl;

    if ((down_boundary.all_elements_less_than(MPs)) and (MPs.all_elements_less_than(up_boundary))) in_range = true;
    return in_range;
}

void Cleng::make_BC() {
//    make boundary condition point => BC = {1,1,1} => minor; BC = {0,0,0} => periodic
    BC = {
            Lat[0]->BC[0] == "mirror",
            Lat[0]->BC[2] == "mirror",
            Lat[0]->BC[4] == "mirror",
    };
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

vector<int> makeExcludedArray(int step) {
    vector<int> result;
    for (int i = -step + 1; i < step; i++) result.push_back(i);
    return result;
}

Point Cleng::prepareMove() {

    if (debug) cout << "prepareMove in Cleng" << endl;

    if (axis != -1) {
        if (axis == 1) clamped_move = {-2, 0, 0};
        if (axis == 2) clamped_move = {0, -2, 0};
        if (axis == 3) clamped_move = {0, 0, -2};
    } else {

        clamped_move = {
                rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
                rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
                rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
        };

        if ((delta_step % 2) == 1) {  // 1 3 5
            int id4rm = rand.getInt(0, 2);
            if (id4rm == 0) clamped_move.x = 0;
            else {
                if (id4rm == 1) clamped_move.y = 0;
                else clamped_move.z = 0;
            }
        } else {  // 2 4 6
            int id4rm1 = rand.getInt(0, 2);
            int id4rm2 = rand.getInt(0, 2);

            if (id4rm1 == 0) clamped_move.x = 0;
            else {
                if (id4rm1 == 1) clamped_move.y = 0;
                else clamped_move.z = 0;
            }

            if (id4rm2 == 0) clamped_move.x = 0;
            else {
                if (id4rm2 == 1) clamped_move.y = 0;
                else clamped_move.z = 0;
            }
        }
    }

    return clamped_move;
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


int Cleng::getLastMCS() {

    string filename;
    vector<string> sub;
    int MS_step = 0;

    string infilename = In[0]->name;
    In[0]->split(infilename, '.', sub);
    filename = sub[0];
    filename = In[0]->output_info.getOutputPath() + filename;
    //read kal file
    ifstream infile(filename + ".kal");

    if (infile) {
        string line;
        while (infile >> std::ws && std::getline(infile, line)); // skip empty lines

        std::istringstream iss(line);
        std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                         std::istream_iterator<std::string>());

        MS_step = stoi(results[0]);
    } else cout << "Unable to open kal file.\n";

    return MS_step;
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
        cout << "System for calculation: " << endl;
        for (auto &&n : nodes) cout << n->to_string() << endl;
        cout << endl;

        if (success_) {

            New[0]->Solve(true);
            free_energy_trial = Sys[0]->GetFreeEnergy() - GetN_times_mu();
            if (save_vector) test_vector.push_back(free_energy_trial);

            cout << "Free Energy (c): " << free_energy_current    << endl;
            cout << "            (t): " << free_energy_trial << endl;
//            cout << "trial - current = " << std::to_string(free_energy_trial - free_energy_current) << endl;
            if (is_ieee754_nan( free_energy_trial )) {
                cout << " Sorry, Free Energy is NaN. Termination..." << endl;
                break;
            }

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
                    CP(to_cleng);
                    MakeMove(true);
                    CP(to_segment);
                    rejected++;
                }
            }
        }
        cout << "Monte Carlo attempt: " << MC_attempt << endl;
        cout << "Accepted: # " << accepted << " | " << 100 * (accepted / MC_attempt) << "%" << endl;
        cout << "Rejected: # " << rejected << " | " << 100 * (rejected / MC_attempt) << "%" << endl;

        if ((MC_attempt % delta_save) == 0) {
            WriteOutput(MC_attempt + MCS_checkpoint);
        }
        if (checkpoint_save) checkpoint.updateCheckpoint(simpleNodeList);
    }

    cout << endl;
    cout << "Finally:" << endl;
    cout << "Monte Carlo attempt: " << MC_attempt-1 << endl;
    cout << "Accepted: # " << accepted << " | " << 100 * (accepted / (MC_attempt-1)) << "%" << endl;
    cout << "Rejected: # " << rejected << " | " << 100 * (rejected / (MC_attempt-1)) << "%" << endl;

    if (checkpoint_save) {
        cout << "Saving checkpoint..." << endl;
        checkpoint.saveCheckpoint(simpleNodeList);
    }

    cout << "Have a fun. " << endl;

    return success;
}

void Cleng::WriteOutput(int attempt) {
    if (debug) cout << "WriteOutput in Cleng" << endl;
    PushOutput(attempt);
    New[0]->PushOutput();
    for (int i = 0; i < n_out; i++) Out[i]->WriteOutput(attempt);
    if (cleng_pos) WriteClampedNodeDistance(attempt);
}

void Cleng::PushOutput(int attempt) {
//	int* point;
    for (int i = 0; i < n_out; i++) {
        Out[i]->PointerVectorInt.clear();
        Out[i]->PointerVectorReal.clear();
        Out[i]->SizeVectorInt.clear();
        Out[i]->SizeVectorReal.clear();
        Out[i]->strings.clear();
        Out[i]->strings_value.clear();
        Out[i]->bools.clear();
        Out[i]->bools_value.clear();
        Out[i]->Reals.clear();
        Out[i]->Reals_value.clear();
        Out[i]->ints.clear();
        Out[i]->ints_value.clear();

        if (Out[i]->name == "ana" || Out[i]->name == "kal") {
            Out[i]->push("MC_attempt", attempt);
            Out[i]->push("MCS", MCS);
            Out[i]->push("free_energy", free_energy_current);

        }

        if (Out[i]->name == "ana" ||
            Out[i]->name == "vec") { //example for putting an array of Reals of arbitrary length to output
            string s = "vector;0"; //the keyword 'vector' is used for Reals; the value 0 is the first vector, use 1 for the next etc,
            Out[i]->push("gn", s); //"gn" is the name that will appear in output file
            Out[i]->PointerVectorReal.push_back(
                    Mol[clp_mol]->gn); //this is the pointer to the start of the 'vector' that is reported to output.
            Out[i]->SizeVectorReal.push_back(
                    sizeof(Mol[clp_mol]->gn)); //this is the size of the 'vector' that is reported to output
        }

        if (Out[i]->name == "ana" ||
            Out[i]->name == "pos") { //example for putting 3 array's of Integers of arbitrary length to output
            string s = "array;0";
            Out[i]->push("X", s);

//			 TODO using normal point in output and remove xs, ys, zs
            fillXYZ();
//			point=X.data();
            Out[i]->PointerVectorInt.push_back(xs);
            Out[i]->SizeVectorInt.push_back((int) nodes.size());
            s = "array;1";
            Out[i]->push("Y", s);
//			point = Y.data();
            Out[i]->PointerVectorInt.push_back(ys);
            Out[i]->SizeVectorInt.push_back((int) nodes.size());
            s = "array;2";
            Out[i]->push("Z", s);
//			point=Z.data();
            Out[i]->PointerVectorInt.push_back(zs);
            Out[i]->SizeVectorInt.push_back((int) nodes.size());
        }
    }
}

void Cleng::fillXYZ() {
    delete[] xs;
    delete[] ys;
    delete[] zs;

    int n = nodes.size();
    xs = new int[n];
    ys = new int[n];
    zs = new int[n];
    for (int i = 0; i < n; i++) {
        xs[i] = nodes[i]->point().x;
        ys[i] = nodes[i]->point().y;
        zs[i] = nodes[i]->point().z;
    }
}

void Cleng::WriteClampedNodeDistance(int MS_step) {
    vector<Real> distPerMC;
    int i = 0;
    for (auto &&SN :  simpleNodeList) {
        if (!(i % 2)) distPerMC.push_back(SN->distance(SN->get_cnode()->get_system_point()));
        i++;
    }

    string filename;
    vector<string> sub;

    string infilename = In[0]->name;
    In[0]->split(infilename, '.', sub);
    filename = sub[0];
    filename = In[0]->output_info.getOutputPath() + filename;

    ofstream outfile;
    outfile.open(filename + ".cpos", std::ios_base::app);

    outfile << MS_step << " ";
    for (auto n : distPerMC)outfile << n << " ";
    outfile << endl;
}
