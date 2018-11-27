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
) : name(name_),
    In(In_),
    Lat(Lat_),
    Seg(Seg_),
    Sta(Sta_),
    Rea(Rea_),
    Mol(Mol_),
    Sys(Sys_),
    New(New_) {
    if (debug) cout << "Cleng initialized" << endl;
    KEYS.push_back("MCS");
    KEYS.push_back("save_interval");
    KEYS.push_back("save_filename");
    KEYS.push_back("seed");

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
    bool success = true;

    success = In[0]->CheckParameters("cleng", name, start, KEYS, PARAMETERS, VALUES);
    if (success) {
        vector<string> options;
        if (GetValue("MCS").size() > 0) {
            success = In[0]->Get_int(GetValue("MCS"), MCS, 1, 10000,
                                     "The number of timesteps should be between 1 and 10000");
        }
        if (debug) cout << "MCS is " << MCS << endl;

        if (GetValue("save_interval").size() > 0) {
            success = In[0]->Get_int(GetValue("save_interval"), save_interval, 1, MCS,
                                     "The save interval nr should be between 1 and 100");
        }
        if (debug) cout << "Save_interval " << save_interval << endl;

        if (Sys[0]->SysClampList.size() < 1) {
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
        int length = In[0]->MolList.size();
        for (int i = 0; i < length; i++) if (Mol[i]->freedom == "clamped") clp_mol = i;
    }
    if (success) {
        n_out = In[0]->OutputList.size();
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

SimpleNode fromSystemToNode(int x, int y, int z, int id, const Point &box) {
    return SimpleNode({x, y, z}, id, box);
}

vector<shared_ptr<Node>> createNodes(const vector<SimpleNode> &simple_nodes) {
    vector<shared_ptr<Node>> result;
    map<SimpleNode, vector<SimpleNode>> m;
    for (auto &&n  : simple_nodes) {
        m[n].push_back(n);
    }
    for (auto &&entry : m) {
        if (entry.second.size() == 1) {
            result.push_back(make_shared<SimpleNode>(entry.first));
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

    vector<SimpleNode> sn;
    Segment *clamped = Seg[clamp_seg];
    map<int, Point> system_points;
    switch (tofrom) {
        case to_cleng:
            for (int i = 0; i < n_boxes; i++) {
                sn.push_back(fromSystemToNode(clamped->px1[i], clamped->py1[i], clamped->pz1[i], 2 * i, box));
                sn.push_back(fromSystemToNode(clamped->px2[i], clamped->py2[i], clamped->pz2[i], 2 * i + 1, box));
            }

            nodes = createNodes(sn);
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

void Cleng::WriteOutput(int MS_step) {
    if (debug) cout << "WriteOutput in Cleng" << endl;

    PushOutput(MS_step);
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

bool Cleng::CheckRange_and_distances() {
    double min_dist = 9;
    bool result_distances = true;
    bool result_range = false;
    int up_boundary = 3; int down_boundary = 3;
    Point box_wr { Lat[0]->MX - up_boundary, Lat[0]->MY - up_boundary, Lat[0]->MZ - up_boundary };
    Point MP {
        nodes[id_part_for_move]->point().x,
        nodes[id_part_for_move]->point().y,
        nodes[id_part_for_move]->point().z
    };
    Point MPs {
        nodes[id_part_for_move]->point().x + shift.x,
        nodes[id_part_for_move]->point().y + shift.y,
        nodes[id_part_for_move]->point().z + shift.z
    };

//    cout << "Try MP : " << MP.to_string() << endl;

    // Distances
    int i = 0;
    for (const auto &n : nodes) {
        if (MP != n->point()) {
            double dx = pow(MPs.x - n->point().x, 2);
            double dy = pow(MPs.y - n->point().y, 2);
            double dz = pow(MPs.z - n->point().z, 2);
            double dr = sqrt(dx+dy+dz);
            if (dr < min_dist) i++;
        }
    }
    if (i > 0) {
        result_distances = false;
    }
//    result ? cout<<endl : cout << "Distances are not okay...\n" << endl;

    // Range
    if (
            (down_boundary < MPs.x) and (MPs.x < box_wr.x) and
            (down_boundary < MPs.y) and (MPs.y < box_wr.y) and
            (down_boundary < MPs.z) and (MPs.z < box_wr.z)
            )
    {
        result_range = true;
    }
    return result_distances and result_range;
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
        id_part_for_move = GetIntRandomValueExclude(0, (int) nodes.size() - 1, 0, false);
        cout << "Trying: \n part_id: " << id_part_for_move
             << ", MC step: { " << shift.x << ", " << shift.y << ", " << shift.z << " }" << endl;

        if (CheckRange_and_distances()) {
            nodes[id_part_for_move]->shift(shift);
        } else {
            while (!CheckRange_and_distances()) {
                cout << "Choose another particle id and step..." << endl;
                id_part_for_move = GetIntRandomValueExclude(0, (int) nodes.size() - 1, 0, false);
                shift = PrepareStep();
            }
            nodes[id_part_for_move]->shift(shift);
        }
        cout << "Finally: \n part_id: " << id_part_for_move
        << ", MC step: { " << shift.x << ", " << shift.y << ", " << shift.z << " }" << endl;
    } else {
        cout << "MakeShift back" << endl;
        Point neg_shift = shift.negate();
        nodes[id_part_for_move]->shift(neg_shift);
    }

    return success;
}

bool Cleng::MonteCarlo() {
    if (debug) cout << "Monte Carlo in Cleng" << endl;
    bool success = true;

    Real free_energy_c;
    Real free_energy_t;

// init system outlook
    New[0]->Solve(true);
    WriteOutput(0);

    for (int MS_step = 1; MS_step < MCS; MS_step++) { // loop for trials

        Real my_rand = GetRealRandomValue(0, 1);
        free_energy_c = Sys[0]->FreeEnergy;
        CP(to_cleng);
        MakeShift(false);
        CP(to_segment);
        New[0]->Solve(true);

        free_energy_t = Sys[0]->FreeEnergy;

        cout << "my_rand:" << my_rand << endl;
        cout << "free_energy_c:" << free_energy_c << endl;
        cout << "free_energy_t:" << free_energy_t << endl;

        if (std::isnan(free_energy_t)) {
            cout << "Free Energy is nan!!" << endl;
            MakeShift(true);
            break;
        }


        if (free_energy_t <= free_energy_c) {
            cout << "Accepted" << endl;
        } else {
            if (my_rand < exp(free_energy_t - free_energy_c)) {
                cout << "Accepted" << endl;
            } else {
                cout << "Deny" << endl;
                MakeShift(true);
                cout << "Done. No saving. \n" << endl;
                continue;
            }
        }

        if ((MS_step % save_interval) == 0) {
            WriteOutput(MS_step);
        }
        cout << "Done. \n" << endl;
    }
    return success;
}

void Cleng::PutParameter(string new_param) {
    KEYS.push_back(new_param);
}

string Cleng::GetValue(string parameter) {
    int i = 0;
    int length = PARAMETERS.size();
    while (i < length) {
        if (parameter == PARAMETERS[i]) {
            return VALUES[i];
        }
        i++;
    }
    return "";
}

void Cleng::PushOutput(int MS_step) {
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
            Out[i]->SizeVectorInt.push_back(nodes.size());
            s = "array;1";
            Out[i]->push("Y", s);
//			point = Y.data();
            Out[i]->PointerVectorInt.push_back(ys);
            Out[i]->SizeVectorInt.push_back(nodes.size());
            s = "array;2";
            Out[i]->push("Z", s);
//			point=Z.data();
            Out[i]->PointerVectorInt.push_back(zs);
            Out[i]->SizeVectorInt.push_back(nodes.size());
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
