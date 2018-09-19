#include <utility>

#include "cleng.h"
#include "nodes/point.h"
#include "nodes/monolit.h"
#include <map>
#include <cassert>

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
    KEYS.emplace_back("save_interval");
    KEYS.emplace_back("save_filename");
    KEYS.emplace_back("seed");

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


bool Cleng::CheckInput(int start) {
    if (debug) cout << "CheckInput in Cleng" << endl;
    bool success;

    success = In[0]->CheckParameters("cleng", name, start, KEYS, PARAMETERS, VALUES);
    if (success) {
        vector<string> options;
        if (!GetValue("MCS").empty()) {
            success = In[0]->Get_int(GetValue("MCS"), MCS, 1, 10000,
                                     "The number of timesteps should be between 1 and 10000");
        }
        if (debug) cout << "MCS is " << MCS << endl;

        if (!GetValue("save_interval").empty()) {
            success = In[0]->Get_int(GetValue("save_interval"), save_interval, 1, MCS,
                                     "The save interval nr should be between 1 and 100");
        }
        if (debug) cout << "Save_interval " << save_interval << endl;

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

        if (success) {
            n_boxes = Seg[clamp_seg]->n_box;
            sub_box_size = Seg[clamp_seg]->mx;
        }
        clp_mol = -1;
        int length = (int)In[0]->MolList.size();
        for (int i = 0; i < length; i++) if (Mol[i]->freedom == "clamped") clp_mol = i;
    }
    if (success) {
        n_out = (int)In[0]->OutputList.size();
        if (n_out == 0) cout << "Warning: no output defined!" << endl;

        for (int i = 0; i < n_out; i++) {
            Out.push_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[i], i, n_out));
            if (!Out[i]->CheckInput(start)) {
                cout << "input_error in output " << endl;
                success = false;
            }
        }
        MonteCarlo();
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
            for (auto &&n : nodes) {
                cout << n->to_string() << endl;
            }
//            assert(nodes.size() == 3);
            break;

        case to_segment:
            Zero(Seg[clamp_seg]->H_MASK, M);

            for (auto &&n : nodes) {
                n->pushSystemPoints(system_points);
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
                Seg[clamp_seg]->bx[i] = (Seg[clamp_seg]->px2[i] + Seg[clamp_seg]->px1[i] - sub_box_size) / 2;
                Seg[clamp_seg]->by[i] = (Seg[clamp_seg]->py2[i] + Seg[clamp_seg]->py1[i] - sub_box_size) / 2;
                Seg[clamp_seg]->bz[i] = (Seg[clamp_seg]->pz2[i] + Seg[clamp_seg]->pz1[i] - sub_box_size) / 2;

                if (Seg[clamp_seg]->bx[i] < 1) {
                    Seg[clamp_seg]->bx[i] += box.x;
                    Seg[clamp_seg]->px1[i] += box.x;
                    Seg[clamp_seg]->px2[i] += box.x;
                }
                if (Seg[clamp_seg]->by[i] < 1) {
                    Seg[clamp_seg]->by[i] += box.y;
                    Seg[clamp_seg]->py1[i] += box.y;
                    Seg[clamp_seg]->py2[i] += box.y;
                }
                if (Seg[clamp_seg]->bz[i] < 1) {
                    Seg[clamp_seg]->bz[i] += box.z;
                    Seg[clamp_seg]->pz1[i] += box.z;
                    Seg[clamp_seg]->pz2[i] += box.z;
                }

                Seg[clamp_seg]->H_MASK[((Seg[clamp_seg]->px1[i] - 1) % box.x + 1) * JX +
                                       ((Seg[clamp_seg]->py1[i] - 1) % box.y + 1) * JY +
                                       (Seg[clamp_seg]->pz1[i] - 1) % box.z + 1] = 1;
                Seg[clamp_seg]->H_MASK[((Seg[clamp_seg]->px2[i] - 1) % box.x + 1) * JX +
                                       ((Seg[clamp_seg]->py2[i] - 1) % box.y + 1) * JY +
                                       (Seg[clamp_seg]->pz2[i] - 1) % box.z + 1] = 1;

            }

            // debug.log
            //out << "Nodes:" << endl;
            //for (auto &&n : nodes) {
            //    out << n->to_string() << endl;
            //}
            //for (int i = 0; i < n_boxes; i++) {
            //    out << "Box " << i << ": ";
            //    out << "{" << clamped->bx[i] << ", " << clamped->by[i] << ", " << clamped->bz[i] << "}" << endl;
            //}
            //out << endl;
//            nodes.clear();
            break;

        default:
            success = false;
            cout << "error in transfer" << endl;
            break;
    }
    return success;
}

void Cleng::WriteOutput(int MS_step, Real exp_diff) {
    if (debug) cout << "WriteOutput in Cleng" << endl;
    PushOutput(MS_step, exp_diff);
    New[0]->PushOutput();
    for (int i = 0; i < n_out; i++) {
        Out[i]->WriteOutput(MS_step);
    }
}

int Cleng::GetIntRandomValueExclude(int min_value, int max_value, int exclude_value, bool need_exclude) {
    if (debug) cout << "Int GetRandomValue in Cleng" << endl;
    int out;
    random_device rd;
    default_random_engine gen(rd());
    uniform_int_distribution<> dist(min_value, max_value);
    out = dist(gen);
    if (need_exclude) {
        while (out == exclude_value) {
            out = dist(gen);
        }
    }
    return out;
}

Real Cleng::GetRealRandomValue(int min_value, int max_value) {
    if (debug) cout << "Real GetRandomValue in Cleng" << endl;
    random_device rd;
    default_random_engine gen(rd());
    uniform_real_distribution<> dist(min_value, max_value);
    return dist(gen);
}

Real Cleng::GetN_times_mu() {
    int n_mol=(int)In[0]->MolList.size();
    Real n_times_mu=0;
    for (int i=0; i<n_mol; i++) {
        Real Mu=Mol[i]->Mu;
        Real n=Mol[i]->n;
        if (Mol[i]->IsClamped()) n=Mol[i]->n_box;
        n_times_mu +=  n*Mu;
    }
    return n_times_mu;
}

void Cleng::makeInSubBoxRange() {
    int not_in_subbox_range = 0;
    Point sub_box = {sub_box_size, sub_box_size, sub_box_size};  // change it in future

    for (auto &&n : nodes) {
        if(!n->inSubBoxRange(sub_box)) not_in_subbox_range ++;
    }
    cout << "not_in_subbox_range: " << not_in_subbox_range << endl;
    if (not_in_subbox_range != 0) {
        cout << "MakeShift back!" << endl;
        MakeShift(true);
        for (auto &&n : nodes) {
            cout << n->to_string() << endl;
        }
    }
}

bool Cleng::NotCollapsing() {
    bool not_collapsing = false;
    Point MP (nodes[id_node_for_move]->point());
    Point MPs (nodes[id_node_for_move]->point() + shift);
    double min_dist = 4; // minimal distance between nodes
    int i = 0;
    for (const auto &n : nodes) {
        if (MP != n->point()) {
            Real distance = MPs.distance(n->point());
            if (distance < min_dist) i++;
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
    Point down_boundary = { 3, 3, 3};
    Point up_boundary   = box - down_boundary;
    Point MPs (nodes[id_node_for_move]->point() + shift);

    cout << "MPs" << MPs.to_string() << endl;
    cout << "down bound" << down_boundary.to_string() << endl;
    cout << "up boundary bound" << up_boundary.to_string() << endl;
    if ((down_boundary.less_all_elements_than(MPs)) and (MPs.less_all_elements_than(up_boundary))) in_range = true;
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
    bool not_collapsing;
    bool in_range;

    // check distances between all nodes => preventing collapsing
    not_collapsing = NotCollapsing();

    // check distance between all nodes and constrains (walls)
    if (BC.x and BC.y and BC.z) { // BC.x, BC.y, BC.z = mirror
        in_range = InRange();
    } else { // BC.x and/or BC.y and/or BC.z != mirror
        in_range = true;
    }
    cout << "not_collapsing: " << not_collapsing << " in_range: " << in_range << endl;
    return not_collapsing and in_range;
}

Point Cleng::PrepareStep() {
    shift = {GetIntRandomValueExclude(-1, 1, 0, true), GetIntRandomValueExclude(-1, 1, 0, true),
             GetIntRandomValueExclude(-1, 1, 0, true)};

    int id_pos_eq0 = GetIntRandomValueExclude(0, 2, 0, false);
    if (id_pos_eq0 == 0) {
        shift.x = 0;
    } else {
        if (id_pos_eq0 == 1) {
            shift.y = 0;
        } else {
            shift.z = 0;
        }
    }
    return shift;
}

bool Cleng::MakeShift(bool back) {
    if (debug) cout << "MakeShift in Cleng" << endl;
    bool success = true;

    if (!back) {
        shift = PrepareStep();
        id_node_for_move = GetIntRandomValueExclude(0, (int) nodes.size() - 1, 0, false);
        cout << "Trying: \n node_id: " << id_node_for_move
             << ", MC step: { " << shift.x << ", " << shift.y << ", " << shift.z << " }" << endl;
        if (Checks()) {
            nodes[id_node_for_move]->shift(shift);
            makeInSubBoxRange();
        } else {
            while (!Checks()) {
                cout << "Choose another particle id and step..." << endl;
                id_node_for_move = GetIntRandomValueExclude(0, (int) nodes.size() - 1, 0, false);
                shift = PrepareStep();
                cout << "Trying: \n node_id: " << id_node_for_move
                     << ", MC step: { " << shift.x << ", " << shift.y << ", " << shift.z << " }" << endl;
            }
            nodes[id_node_for_move]->shift(shift);
            makeInSubBoxRange();
        }
        cout << "Finally: \n node_id: " << id_node_for_move
        << ", MC step: { " << shift.x << ", " << shift.y << ", " << shift.z << " }" << endl;
    } else {
        cout << "MakeShift back" << endl;
        Point neg_shift = shift.negate();
        nodes[id_node_for_move]->shift(neg_shift);
    }

    return success;
}

bool Cleng::MonteCarlo() {
    if (debug) cout << "Monte Carlo in Cleng" << endl;
    bool success = true;
// Analysis MC
    Real accepted=0.0;
    Real rejected=0.0;
    Real exp_diff=0.0;
    make_BC();

// init system outlook
    New[0]->Solve(true);
    free_energy_current = Sys[0]->GetFreeEnergy() - GetN_times_mu();
    WriteOutput(0,exp_diff);

    for (int MS_step = 1; MS_step < MCS; MS_step++) { // loop for trials
        CP(to_cleng);
        MakeShift(false);
        CP(to_segment);
        New[0]->Solve(true);
        free_energy_trial = Sys[0]->GetFreeEnergy() - GetN_times_mu();

        assert(!std::isnan(free_energy_trial));
        for (auto &&n : nodes) {
            cout << "n: " << n->point().to_string() << endl;
            assert((n->point().x > 3) and (n->point().x < 28));
            assert((n->point().y > 3) and (n->point().y < 28));
            assert((n->point().z > 3) and (n->point().z < 28));
        }


        cout << "free_energy_current: " << free_energy_current << endl;
        cout << "free_energy_trial: " << free_energy_trial << endl;

        if (free_energy_trial - free_energy_current <= 0.0) {
            cout << "Accepted" << endl;
            free_energy_current = free_energy_trial;
            accepted ++;
        } else {
            Real acceptance = GetRealRandomValue(0, 1);
            cout << "acceptance: " << acceptance << endl;
            cout << "diff: " << free_energy_trial - free_energy_current << endl;
            exp_diff = exp(-1.0 * (free_energy_trial - free_energy_current));
            cout << "exp(diff): " << exp_diff << endl;

            if (acceptance < exp(-1.0 * (free_energy_trial - free_energy_current))) {
                cout << "Accepted" << endl;
                free_energy_current = free_energy_trial;
                accepted ++;
            } else {
                cout << "Deny" << endl;
                MakeShift(true);
                cout << "Done. No saving. \n" << endl;
                rejected ++;
                continue;
            }
        }
        cout << "MonteCarlo steps: " << MS_step << endl;
        cout << "Accepted %: " << 100* (accepted / MS_step) << endl;
        cout << "Rejected %: " << 100* (rejected / MS_step) << endl;

        if ((MS_step % save_interval) == 0) {
            WriteOutput(MS_step, exp_diff);
        }
        cout << "Done. \n" << endl;
    }

    cout << "Finally:" << endl;
    cout << "Accepted %: " << 100* (accepted / (MCS-1)) << endl;
    cout << "Rejected %: " << 100* (rejected / (MCS-1)) << endl;

    return success;
}

void Cleng::PutParameter(string new_param) {
    KEYS.push_back(new_param);
}

string Cleng::GetValue(string parameter) {
    int i = 0;
    int length = (int)PARAMETERS.size();
    while (i < length) {
        if (parameter == PARAMETERS[i]) {
            return VALUES[i];
        }
        i++;
    }
    return "";
}

void Cleng::PushOutput(int MS_step, Real exp_diff) {
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
            Out[i]->push("MC_step", MS_step);
            Out[i]->push("MCS", MCS);
            Out[i]->push("free_energy_current", free_energy_current);
            Out[i]->push("exp_diff", exp_diff);
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
            Out[i]->SizeVectorInt.push_back((int)nodes.size());
            s = "array;1";
            Out[i]->push("Y", s);
//			point = Y.data();
            Out[i]->PointerVectorInt.push_back(ys);
            Out[i]->SizeVectorInt.push_back((int)nodes.size());
            s = "array;2";
            Out[i]->push("Z", s);
//			point=Z.data();
            Out[i]->PointerVectorInt.push_back(zs);
            Out[i]->SizeVectorInt.push_back((int)nodes.size());
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
