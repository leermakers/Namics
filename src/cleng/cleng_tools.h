#include <utility>
#include "cleng.h"
#include "iterator/EnumerateIterator.h"

volatile sig_atomic_t cleng_flag_termination = 0;

using Fp_info = numeric_limits<double>;

inline auto is_ieee754_nan(double const x) -> bool {
    static constexpr bool is_claimed_ieee754 = Fp_info::is_iec559;
    static constexpr int n_bits_per_byte = CHAR_BIT;
    using Byte = unsigned char;

    static_assert(is_claimed_ieee754, "!");
    static_assert(n_bits_per_byte == 8, "!");
    static_assert(sizeof(x) == sizeof(uint64_t), "!");

#ifdef _MSC_VER
    uint64_t const bits = reinterpret_cast<uint64_t const&>( x );
#else
    Byte bytes[sizeof(x)];
    memcpy(bytes, &x, sizeof(x));
    uint64_t int_value;
    memcpy(&int_value, bytes, sizeof(x));
    uint64_t const &bits = int_value;
#endif

    static constexpr uint64_t sign_mask = 0x8000000000000000;
    static constexpr uint64_t exp_mask = 0x7FF0000000000000;
    static constexpr uint64_t mantissa_mask = 0x000FFFFFFFFFFFFF;

    (void) sign_mask;
    return (bits & exp_mask) == exp_mask and (bits & mantissa_mask) != 0;
}

shared_ptr<SimpleNode> fromSystemToNode(int x, int y, int z, int id, const Point &box) {
    return make_shared<SimpleNode>(Point(x, y, z), id, box);
}

map<int, vector<shared_ptr<Node>>> createNodes(const vector<shared_ptr<SimpleNode>> &simple_nodes, map<int, vector<int>> &pivot_arm_nodes, int &arms) {
    map<int, vector<shared_ptr<Node>>> result_map;
    int index = 0;
    arms = 0;
    map<SimpleNode, vector<shared_ptr<SimpleNode>>> m;
    for (auto &&n : simple_nodes) m[*n].push_back(n);
    for (auto &&entry : m) {
        if (entry.second.size() == 1) result_map[index] = {entry.second[0]};
        else result_map[index] = {make_shared<Monolit>(entry.second)};
        index++;
    }
    // pivot_arm_nodes
    Point P_first = simple_nodes[0]->point();
    for (auto &&n : simple_nodes) {
        int ID_test = n->get_ID();
//        cout << "ID_test: " << ID_test << endl;
        if (n->point() == P_first) arms++;
        for (auto &&node : result_map) {
            if (node.second.data()->get()->isIdInside(ID_test)) {pivot_arm_nodes[arms].push_back(node.first);}
        }
    }
    // clean up the repeating elements
    for (auto &&pair : pivot_arm_nodes){
        int arm_ = pair.first;
        auto last = std::unique(pivot_arm_nodes[arm_].begin(), pivot_arm_nodes[arm_].end());
        pivot_arm_nodes[arm_].erase(last, pivot_arm_nodes[arm_].end());
    }

    return result_map;
}

vector<int> makeExcludedArray(int step) {
    vector<int> result;
    for (int i = -step + 1; i < step; i++) result.push_back(i);
    return result;
}

bool Cleng::InSubBoxRange(int id_node_for_move) {
    bool success = true;
    vector<int> id_nodes;
    Point sub_box = {sub_box_size.x, sub_box_size.y, sub_box_size.z};  // change it in future

    if (!nodes_map[id_node_for_move].data()->get()->inSubBoxRange(sub_box)) {
        id_nodes.push_back(id_node_for_move);
        success = false;
        cout << endl;
        cout << "[WARNING] Subbox check is failed..." << endl;
    }
    return success;
}

bool Cleng::IsCommensuratable() {
    bool success = true;
    int length = (int) In[0]->MolList.size();

    // TODO: procedure goes throw all Molecules in the System...
    // TODO: Need to take only Polymer molecules...
    for (int i = 0; i < length; i++) {
        int chain_length = Mol[i]->chainlength;
        if (chain_length < 10) continue;   // hack

        for (auto &&SN : Enumerate(simpleNodeList)) {
            auto p3 = SN.second->get_cnode()->get_system_point() - SN.second->get_system_point();
            int path_length = abs(p3.x) + abs(p3.y) + abs(p3.z);

            int path_length_even = path_length % 2;
            int chain_length_even = chain_length % 2;

//            cout << "path_length:  "  << path_length
//            << " | chain_length: " << chain_length << endl;

//            if (path_length > 20) cout << "Ooowww!! path_len > 20! "<< endl;
//            if (chain_length > 20) cout << "Ooowww!! chain_len > 20! "<< endl;

            if (path_length_even == chain_length_even) success = false;
            if (path_length >= chain_length) success = false;
        }
    }
    if (!success) {
        cout << endl;
        cout << "[WARNING] paths between clamps is not commensurate with the chain length!" << endl;
    }
    return success;
}

//  Idea: check if shifted_point (candidate) will occupy position with *0* value of KSAM
//
//  (BUG) At very first check KSAM probably zero array. Need to be initialized... TODO
//  a coordinate x,y,z is put on element j:  KSAM[j]=KSAM[JX*x+JY*y+z] with JX=(MX+2)*(MY+2) and JY=(MY+2)
//std::cout << "is it 0/1?" << std::endl; // if it's *0* --> frozen state --> bad move
//std::cout << Sys[0] -> KSAM[JX*shifted_point.x + JY*shifted_point.y + shifted_point.z] << std::endl;
bool Cleng::NotViolatedFrozenStates(int id_node_for_move) {
if (debug) cout <<"NotViolatedFrozenStates " << endl;
    bool not_violated_frozen_states = true;

    //int sum = 0;
    //Sum(sum, Sys[0]->KSAM, Lat[0]->MX*Lat[0]->MY*Lat[0]->MZ);
    //std::cout << sum << std::endl;
    //if (sum == 0) return true; // hack;

    // Apparently these quantities are available
    int JX = Lat[0]->JX;
    int JY = Lat[0]->JY;

    //cout << "JX:" << JX << endl;
    //cout << "JY:" << JY << endl;

    Point shifted_point(nodes_map[id_node_for_move].data()->get()->point());
    //cout << "id node for move: " << id_node_for_move << " | coord: " << shifted_point.to_string() << endl;

    // additionally check the hypothetical nodes between node and clamped node
    int i = 0;
    for (auto &&SN :  simpleNodeList) {
        i++;
        if  (i%2 == 0) continue;
        //cout << i << " node "<< SN->get_system_point().to_string() << endl;
        //cout << i << " cnode"<< SN->get_cnode()->get_system_point().to_string() << endl;

        //int s_x = SN->get_system_point().x < SN->get_cnode()->get_system_point().x ?SN->get_system_point().x:SN->get_cnode()->get_system_point().x;
        //int b_x = SN->get_system_point().x > SN->get_cnode()->get_system_point().x ?SN->get_system_point().x:SN->get_cnode()->get_system_point().x;
        //int s_y = SN->get_system_point().y < SN->get_cnode()->get_system_point().y ?SN->get_system_point().y:SN->get_cnode()->get_system_point().y;
        //int b_y = SN->get_system_point().y > SN->get_cnode()->get_system_point().y ?SN->get_system_point().y:SN->get_cnode()->get_system_point().y;
        //int s_z = SN->get_system_point().z < SN->get_cnode()->get_system_point().z ?SN->get_system_point().z:SN->get_cnode()->get_system_point().z;
        //int b_z = SN->get_system_point().z > SN->get_cnode()->get_system_point().z ?SN->get_system_point().z:SN->get_cnode()->get_system_point().z;

		int s_x = SN->get_system_point().x;
        int b_x = SN->get_cnode()->get_system_point().x;
        int s_y = SN->get_system_point().y;
        int b_y = SN->get_cnode()->get_system_point().y;
        int s_z = SN->get_system_point().z ;
        int b_z = SN->get_cnode()->get_system_point().z;

        //cout << "s_x:" << s_x << " b_x:" << b_x << endl;
        //cout << "s_y:" << s_y << " b_y:" << b_y << endl;
        //cout << "s_z:" << s_z << " b_z:" << b_z << endl;
        // it is possible to extend the range of checked hypothetical nodes
		int MX=Lat[0]->MX;
		int MY=Lat[0]->MY;
		int MZ=Lat[0]->MZ;
		//bool MirrorX=false;
		//bool MirrorY=false;
		//bool MirrorZ=false;
		//if (Lat[0]->BC[0]=="mirror") {MirrorX=true; cout <<"mirrorx"<< endl;}
		//if (Lat[0]->BC[1]=="mirror") {MirrorY=true; cout <<"mirrory"<< endl;}
		//if (Lat[0]->BC[2]=="mirror") {MirrorZ=true; cout <<"mirrorz"<< endl;}
		if (s_x>MX ) s_x-=MX;
		if (b_x>MX ) b_x-=MX;
		if (s_y>MY ) s_y-=MY;
		if (b_y>MY ) b_y-=MY;
		if (s_z>MZ ) s_z-=MZ;
		if (b_z>MZ ) b_z-=MZ;

        //cout << " s_x:" << s_x << " b_x:" << b_x << endl;
        //cout << " s_y:" << s_y << " b_y:" << b_y << endl;
        //cout << " s_z:" << s_z << " b_z:" << b_z << endl;

        //for (int _x = s_x; _x < b_x+1; _x++) {
            //for (int _y = s_y; _y < b_y+1; _y++)  {
                //for (int _z = s_z; _z < b_z+1; _z++)  {
                    //cout << "{ " << _x << " " << _y << " " << _z << " }" << endl;
                    //cout << "KSAM[JX*x+JY*y+z] JX:" << JX << "| JY:" << JY << "|:" << Sys[0] -> KSAM[JX*_x + JY*_y + _z]  << endl;

                    int length = In[0]->MonList.size();
                    for (int i = 0; i < length; i++)
                    {

                        if (Seg[i]->freedom == "frozen") {
                            //cout << "Seg -> MASK[JX*x+JY*y+z]: " << Seg[i]->MASK[JX*_x + JY*_y + _z] << endl;
                            if (Seg[i] -> MASK[JX*s_x + JY*s_y + s_z] == 1 || Seg[i] -> MASK[JX*b_x + JY*b_y + b_z] == 1) {
                                //cout << "node on  FROZEN site at";
                                //cout <<" x,y,z = " << s_x << "," << s_y << "," << s_z << endl;

                                not_violated_frozen_states = false;
                                return not_violated_frozen_states;
                            }
                            //if (Seg[i] -> MASK[JX*b_x + JY*b_y + b_z] == 1) {
                            //    cout << "node on  FROZEN Site at";
                            //    cout <<" x,y,z = " << b_x << "," << b_y << "," << b_z << endl;

                            //    not_violated_frozen_states = false;
                            //    return not_violated_frozen_states;
                            //}
                        }
                    }

                    // Apparently KSAM could be *0* for non frozen state... TODO: need to check...
                    //
                    //if (Sys[0] -> KSAM[JX*_x + JY*_y + _z] == 0) {
                    //    // either it can be frozen state or another node -> checking for another node
                    //    Point p = Point(_x, _y, _z);

                    //    bool is_it_node = false;
                    //    int _i = 0;
                    //    for (auto &&_SN: simpleNodeList) {
                    //        _i++;
                    //        if (_i%2 == 0) continue;

                    //        Point p1 = _SN->point();
                    //        Point p2 = _SN->get_cnode()->point();
                    //        if ((p1 == p) || (p2 == p))
                    //            is_it_node = true;
                    //    }

                    //    cout << "is it node: " << is_it_node << endl;
                    //    //cout << "POINT _x _y _z " << p.to_string() << endl;

                    //    if (!is_it_node) {
                    //        cout << endl;
                    //        cout << "Virtual node "<< "{ "<< _x << " " << _y << " " << _z << " }"
                    //             << " is crossing frozen state!"
                    //             << endl;
                    //        not_violated_frozen_states = false;
                    //        return not_violated_frozen_states;
                    //    }
                    //}
                //}
            //}
        //}
        //bool less = SN->get_system_point() < SN->get_cnode()->get_system_point();
        //cout << i << " node < cnode: " << less << endl;
    }
    return not_violated_frozen_states;
}

bool Cleng::NotCollapsing(int id_node_for_move) {
    bool not_collapsing = true;

    Point shifted_point(nodes_map[id_node_for_move].data()->get()->point());
    //cout << "id node for move: " << id_node_for_move << " | coord: " << shifted_point.to_string() << endl;

    double min_dist = 0.3; // minimal distance between nodes_map. at least 1. --> too dense --> crash state TODO: write test!

    for (auto &&n : nodes_map) {
        Point test_point  = n.second.data()->get()->point();
        int test_point_id = n.first;

        if ( (shifted_point != test_point) ) { // not the same as shifted point
            Real distance = shifted_point.distance(n.second.data()->get()->point());
            //cout << "Id test [check]: " << test_point_id << " | coord: " << test_point.to_string() << " <-- Checking"  << endl;

            if (distance <= min_dist) {
                cout << endl;
                cout << "Nodes are too close to each other." << endl;
                cout << "Shifted point from nodes_map: "     << shifted_point.to_string() << endl;
                cout << "Test point from nodes_map:    "     << test_point.to_string()    << endl;
                not_collapsing = false;
            }
        } else {
            if (id_node_for_move != test_point_id) {
                Real distance = shifted_point.distance(n.second.data()->get()->point());
                if (distance <= min_dist) {
                    cout << endl;
                    cout << "Nodes are too close to each other." << endl;
                    cout << "Shifted point from nodes_map: "     << shifted_point.to_string() << endl;
                    cout << "Test point from nodes_map:    "     << test_point.to_string()    << endl;
                    not_collapsing = false;
                }
            }
           // else {
           //     cout << "Same point: " << "Id test [check]: " << test_point_id << " | coord: " << test_point.to_string() << endl;
           // }
        }
    }
    return not_collapsing;
}

bool Cleng::InRange(int id_node_for_move) {
    bool in_range = false;
    Point down_boundary = {1, 1, 1};
    Point up_boundary = box - down_boundary;
    Point shifted_point(nodes_map[id_node_for_move].data()->get()->point());
    if ((down_boundary.all_elements_less_than(shifted_point)) and
        (shifted_point.all_elements_less_than(up_boundary)))
        in_range = true;
    return in_range;
}

void signalHandler(int signum) {
    cout << "---> Termination..." << endl;
    cleng_flag_termination++;
    if (cleng_flag_termination > 1) exit(signum);
}

void Cleng::make_BC() {
//    make boundary condition point => BC = {1,1,1} => minor; BC = {0,0,0} => periodic
    BC = {
            Lat[0]->BC[0] == "mirror",
            Lat[0]->BC[2] == "mirror",
            Lat[0]->BC[4] == "mirror",
    };
}

vector<Real> Cleng::prepare_vtk(string key, string name, string prop) {
    int Size = 0;
    vector<Real> vtk;
    vtk.clear();
    Real *X = Out[0]->GetPointer(std::move(key), std::move(name), std::move(prop), Size);
    if (X != nullptr) {
        for (int i = 1; i < box.x + 1; i++)
        for (int j = 1; j < box.y + 1; j++)
        for (int k = 1; k < box.z + 1; k++)
        vtk.push_back(X[i * J.x + j * J.y + k]);
    } else {
        vtk.push_back(0.0);
    }

    return vtk;
}

int Cleng::prepareIdNode() {
    if (debug) cout << "prepareIdNode in Cleng" << endl;
    int id_node_for_move = 0;
    if (ids_node4move[0] != -1) { id_node_for_move = ids_node4move[rand.getInt(0, (int) ids_node4move.size() - 1)]; }
    else { id_node_for_move = rand.getInt(0, (int) nodes_map.size() - 1); }

    if (id_node_for_move > (int) nodes_map.size() - 1) {
        cout << "###" << endl;
        cout << "Node id is too high. Trying to move node id: " << id_node_for_move << "." << endl;
        cout << "Available nodes_map: " << endl;
        for (auto &&n: nodes_map)
            cout << "id: " << n.first << " " << n.second.data()->get()->to_string() << endl;
        cout << "Termination..." << endl;
        exit(1);
    }

    return id_node_for_move;
}

void Cleng::prepareIdsNode() {
    if (debug) cout << "prepareIdsNode in Cleng" << endl;
    pivot_node_ids.clear();
    int _arm = rand.getInt(1, pivot_arms);
    int _start_index = rand.getInt(0, pivot_arm_nodes[_arm].size()-2);
    for (size_t i = _start_index; i<pivot_arm_nodes[_arm].size(); i++)
        pivot_node_ids.push_back(pivot_arm_nodes[_arm][i]);
//    for (auto &&iD: pivot_node_ids) cout << "id: " << iD << endl;
}

Point Cleng::prepareMove(const string& type_move) {
    if (debug) cout << "prepareMove in Cleng" << endl;
    Point clamped_move;

    if (type_move == "pivot_move") {
        prepareIdsNode();
        prepareRotationMatrix<int>();
//        prepareScalingMatrix<Real>();

    } else {
        if (type_move == "pivot_one_bond_move") {prepareIdsNode();}
        if (axis) {
            int c1 = 1;  // sing of movement +1 or -1
            int c2 = 2;  // MC step by default 2; could be 1 in two_ends_extension mode

            if (two_ends_extension) c2 = 1;
            if (sign_move == "-") c1 = -1;
            // currently it is implemented for step
            if (axis == 1) clamped_move = {c1 * c2, 0, 0};
            if (axis == 2) clamped_move = {0, c1 * c2, 0};
            if (axis == 3) clamped_move = {0, 0, c1 * c2};
        } else {
            clamped_move = {
                    rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
                    rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
                    rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
            };

            if ((delta_step % 2) == 1) {  // 1 3 5
                int id4rm = rand.getInt(0, 2);
                if (id4rm == 0) clamped_move.x = 0;
                else { if (id4rm == 1) clamped_move.y = 0; else clamped_move.z = 0; }
            } else {  // 2 4 6
                int id4rm1 = rand.getInt(0, 2);
                int id4rm2 = rand.getInt(0, 2);

                if (id4rm1 == 0) clamped_move.x = 0;
                else { if (id4rm1 == 1) clamped_move.y = 0; else clamped_move.z = 0; }

                if (id4rm2 == 0) clamped_move.x = 0;
                else { if (id4rm2 == 1) clamped_move.y = 0; else clamped_move.z = 0; }
            }
        }
    }
    return clamped_move;
}

template<class T>
ClMatrix<T> Cleng::prepareRotationMatrix() {
    if (debug) cout << "prepareRotationMatrix in Cleng" << endl;
    int pivot_coef = rand.getIntExcludeArray(-1, 1, {0});
    if (pivot_axis == -1) {
        int pivot_axis_current = rand.getInt(1, 3);
        rotation_matrix = _create_rotational_matrix<T>(pivot_axis_current, pivot_move*pivot_coef);
    } else rotation_matrix = _create_rotational_matrix<T>(pivot_axis, pivot_move*pivot_coef);
    return rotation_matrix;
}

template<class T>
ClMatrix<T> Cleng::prepareScalingMatrix() {
    if (debug) cout << "prepareScalingMatrix in Cleng" << endl;
    vector<Real> scaling_coeff = {0.5, 1.0, 2.0};
    int scaling_index = rand.getInt(0, 2);
    scaling_matrix = _create_scaling_matrix<Real>(scaling_coeff[scaling_index]);

    cout << "Scaling matrix" << endl;
    cout << scaling_matrix << endl;
    return scaling_matrix;
}

template<class T>
ClMatrix<T> Cleng::_create_rotational_matrix(int axis_rotation, int grad) {
    ClMatrix<T> rotation_matrix_(3, 3);
    switch (axis_rotation) {
        case 1:
            rotation_matrix_.put(0, 0, 1);
            rotation_matrix_.put(0, 1, 0);
            rotation_matrix_.put(0, 2, 0);

            rotation_matrix_.put(1, 0, 0);
            rotation_matrix_.put(1, 1, cos(grad * PIE / 180.0));
            rotation_matrix_.put(1, 2, -sin(grad * PIE / 180.0));

            rotation_matrix_.put(2, 0, 0);
            rotation_matrix_.put(2, 1, sin(grad * PIE / 180.0));
            rotation_matrix_.put(2, 2, cos(grad * PIE / 180.0));

            break;
        case 2:
            rotation_matrix_.put(0, 0, cos(grad * PIE / 180.0));
            rotation_matrix_.put(0, 1, 0);
            rotation_matrix_.put(0, 2, sin(grad * PIE / 180.0));

            rotation_matrix_.put(1, 0, 0);
            rotation_matrix_.put(1, 1, 1);
            rotation_matrix_.put(1, 2, 0);

            rotation_matrix_.put(2, 0, -sin(grad * PIE / 180.0));
            rotation_matrix_.put(2, 1, 0);
            rotation_matrix_.put(2, 2, cos(grad * PIE / 180.0));

            break;
        case 3:
            rotation_matrix_.put(0, 0, cos(grad * PIE / 180.0));
            rotation_matrix_.put(0, 1, -sin(grad * PIE / 180.0));
            rotation_matrix_.put(0, 2, 0);

            rotation_matrix_.put(1, 0, sin(grad * PIE / 180.0));
            rotation_matrix_.put(1, 1, cos(grad * PIE / 180.0));
            rotation_matrix_.put(1, 2, 0);

            rotation_matrix_.put(2, 0, 0);
            rotation_matrix_.put(2, 1, 0);
            rotation_matrix_.put(2, 2, 1);
            break;

        default:
            cout << "[Warning] Strange axis: " << axis_rotation << endl;
    }
    return rotation_matrix_;
}

template<class T>
ClMatrix<T> Cleng::_create_scaling_matrix(Real scaling_coef) {
    ClMatrix<T> scaling_matrix_(3, 3);
    scaling_matrix_.put(0, 0, 1*scaling_coef);
    scaling_matrix_.put(0, 1, 0);
    scaling_matrix_.put(0, 2, 0);

    scaling_matrix_.put(1, 0, 0);
    scaling_matrix_.put(1, 1, 1*scaling_coef);
    scaling_matrix_.put(1, 2, 0);

    scaling_matrix_.put(2, 0, 0);
    scaling_matrix_.put(2, 1, 0);
    scaling_matrix_.put(2, 2, 1*scaling_coef);
    return scaling_matrix_;
}

int Cleng::getLastMCS() const {

    int MC_attemp = 0;

    //read kal file
    ifstream infile(filename + ".kal");

    if (infile) {
        string line;
        while (infile >> ws && getline(infile, line)); // skip empty lines

        istringstream iss(line);
        vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());
        MC_attemp = stoi(results[0]);
    } else cout << "Unable to open kal file.\n";

    return MC_attemp;
}

void Cleng::WriteOutput(int num) {
    if (debug) cout << "WriteOutput in Cleng" << endl;
    PushOutput(num);
    New[0]->PushOutput();
    for (int i = 0; i < n_out; i++) Out[i]->WriteOutput(num);
}

void Cleng::PushOutput(int num) {

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
            Out[i]->push("MC_attempt", num);
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
            Out[i]->SizeVectorInt.push_back((int) nodes_map.size());
            s = "array;1";
            Out[i]->push("Y", s);
//			point = Y.data();
            Out[i]->PointerVectorInt.push_back(ys);
            Out[i]->SizeVectorInt.push_back((int) nodes_map.size());
            s = "array;2";
            Out[i]->push("Z", s);
//			point=Z.data();
            Out[i]->PointerVectorInt.push_back(zs);
            Out[i]->SizeVectorInt.push_back((int) nodes_map.size());
        }
    }
}

void Cleng::fillXYZ() {
    delete[] xs;
    delete[] ys;
    delete[] zs;

    int n = nodes_map.size();
    xs = new int[n];
    ys = new int[n];
    zs = new int[n];
    for (int i = 0; i < n; i++) {
        xs[i] = nodes_map[i].data()->get()->point().x;
        ys[i] = nodes_map[i].data()->get()->point().y;
        zs[i] = nodes_map[i].data()->get()->point().z;
    }
}

void Cleng::WriteClampedNodeDistance(int num) {
    vector<Real> distPerMC;
    int i = 0;
    for (auto &&SN :  simpleNodeList) {
        if (!(i % 2)) distPerMC.push_back(SN->distance(SN->get_cnode()->get_system_point()));
        i++;
    }

    ofstream outfile;
    outfile.open(filename + ".cdis", std::ios_base::app);
    outfile << num << " ";
    for (auto n : distPerMC) outfile << n << " ";
    outfile << endl;
}

void Cleng::WriteClampedNodeDistanceWithCenter(int num) {
    vector<Real> distPerMC;
    int central_node_id = pivot_arm_nodes[1].begin()[0];
    Point central_node  = nodes_map[central_node_id].data()->get()->point();
    for (auto &&pair_arm_nodes : pivot_arm_nodes) {
        for (auto &&node_id : pair_arm_nodes.second) {
            if ( central_node_id != node_id ) {
                Point node = nodes_map[node_id].data()->get()->point();
                distPerMC.push_back(central_node.distance(node));
            }
        }
    }
    // all distances are collected
    ofstream outfile;
    outfile.open(filename + ".cdis2c", std::ios_base::app);
    outfile << num << " " << tpure_simulation << " ";
    for (auto &&value : distPerMC) outfile << value << " ";
    outfile << endl;
}

void Cleng::WriteClampedNodePosition(int num) {

    ofstream outfile;
    outfile.open(filename + ".cpos", std::ios_base::app);

//    cout
    outfile
            << "#step " << num << " {X, Y, Z} #" << endl;

    int i = 0;
    for (auto &&SN :  simpleNodeList) {
        if (!(i % 2)) {
//            cout
            outfile
                    << SN->get_system_point().to_string() << SN->get_cnode()->get_system_point().to_string() << endl;
//            << SN->get_system_point().to_string() << endl;
        }
        i++;
    }
}

void Cleng::Write2File(int step, const string& what, Real value) const {

    ofstream outfile;
    outfile.open(filename + "." +what, std::ios_base::app);

//    cout
    outfile << step << " " << value << endl;
}

void Cleng::Write2File(int step, const string& what, const vector<Real>& values, bool same_line) const {

    ofstream outfile;
    outfile.open(filename + "." +what, std::ios_base::app);

    string delim = (same_line ? " " : "\n");
//    cout
    outfile << step << delim;
    for (auto const& v : values) {
        outfile << v << delim;
    }
    outfile << endl;
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

string Cleng::GetValue(const string& parameter) {
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

void Cleng::update_ids_node4move() {
    if (!ids_node4fix.empty()) {
        ids_node4move.clear();
        for (auto const &n_: nodes_map) {
            bool exists = any_of(begin(ids_node4fix), end(ids_node4fix), [&](int i) { return i == n_.first; });
            if (!exists) ids_node4move.push_back(n_.first);
        }
    }
}


Real Cleng::getFreeEnergyBox() {
    Real F_total  = 0.0;
    Real dist     = 0.0;
    Real chi      = 0.0;
    //
    int mon_pol_id    = 0;
    int mon_sol_id    = 0;
    //int mon_clamp_id  = 0;
    //
    int mol_pol_id = 0;
    //int mol_sol_id = 0;
    int chain_len = 0;
    // through arms
    for (auto &&arm : pivot_arm_nodes) {
 //       cout << "arm: " << arm.first << endl;

        // through arm nodes
        for (int idx=1; idx < (int)arm.second.size(); ++idx)
        {
            int id1 = arm.second[idx-1];
            int id2 = arm.second[idx];

            //
            dist = nodes_map[id1].data()->get()->point().distance(nodes_map[id2].data()->get()->point());

            // find polymer segment [typically it A] // TODO: implement interface for providing polymer mon
            for (size_t seg_idx=0; seg_idx < Seg.size(); seg_idx++)
            {
                if (Seg[seg_idx]->name == "W") mon_sol_id   = seg_idx;
                if (Seg[seg_idx]->name == "A") mon_pol_id   = seg_idx;
                //if (Seg[seg_idx]->name == "X") mon_clamp_id = seg_idx;
            }
            chi = Seg[mon_pol_id]->chi[mon_sol_id];

            // len
            for (size_t mol_idx=0; mol_idx < Mol.size(); mol_idx++)
            {
                if (Mol[mol_idx]->name == "pol")   mol_pol_id = mol_idx;
                //if (Mol[mol_idx]->name == "water") mol_sol_id = mol_idx;
            }
            chain_len = Mol[mol_pol_id]->chainlength;
            //cout << "chain_len:" << chain_len << endl;
            Real proposed_free_energy = calcFreeEnergyBox(chain_len-2,dist,chi);
            F_total += proposed_free_energy;
        }
    }
    return F_total;
}

Real Cleng::calcFreeEnergyBox(const Real& N, const Real& R, const Real& chi) {
    Real F      = 0.0;
    Real F_conf = 0.0;
    Real F_int  = 0.0;
    Real kT = 1;
    Real nu = 1;
    if (chi==0.0) nu = 0.9;
    if (chi==0.5) nu = 0.9;
    if (chi>0.5)  nu = 0.9;  // 0.3
    //F_conf = N * ( (pow(R,2) / (2.0*pow(N,2)) ) - log( 1.0 - (pow(R,2)/pow(N,2)) ) ) / kT;
    F_conf = 3.0/2.0 * pow(R,2) / N;
//    cout << "F_conf:" << F_conf << endl;

    Real V = pow(R * pow(N, nu), 3);
    //cout << "V:" << V << endl;

    F_int  = V * ( ( (1-(N/V))*log(1.0 - (N/V)) ) + chi * N/V * (1.0 - (N/V) )) / kT;
//    cout << "F_int:" << F_int << endl;
    F = F_conf + F_int;
    return F;
}

bool Cleng::solveAndCheckFreeEnergy() {
    bool success = true;
    bool success_iteration = New[0]->Solve(true);
    // breakpoint of free energy value
//    //// Simulation without rescue procedure --->
//    if (is_ieee754_nan(Sys[0]->GetFreeEnergy())) {
//        cout << "#?# Sorry, Free Energy is NaN. " << endl;
//        cout << "#?# Here is result from solver: " << success_iteration << endl;
//
//        cout << "#?# The step will be rejected! "
//                "Probably your system is too dense! "
//                "Simulation will stop... " << endl;
//        cout << internal_name << "[CRASH STATE] " << "returning back the system configuration... " << endl;
//        success = false;
//    }
//    //// Simulation without rescue procedure <---

    //// Simulation with rescue procedure --->
    if (is_ieee754_nan(Sys[0]->GetFreeEnergy())) {
        cout << "%?% Sorry, Free Energy is NaN.  " << endl;
        cout << "%?% Here is result from solver: " << success_iteration << endl;

        // New[0]->attempt_DIIS_rescue("none");  // before the function expects char[]
        New[0]->attempt_DIIS_rescue();
        rescue_times ++;

        cout << "%?% Restarting iteration." << endl;
        success_iteration = New[0]->Solve(true);

        if (is_ieee754_nan(Sys[0]->GetFreeEnergy())) {
            cout << "%?% Sorry, Free Energy is still NaN. " << endl;
            cout << "%?% Here is result from solver: " << success_iteration << endl;

            cout << "%?% The step will be rejected! "
                    "Probably your system is too dense! "
                    "Simulation will continue... " << endl;
            cout << internal_name << "[CRASH STATE] " << "returning back the system configuration... " << endl;
            cleng_rejected++;

            success = false;
        }
    }
    //// Simulation with rescue procedure <---

    return success;
}

bool Cleng::initSystemOutlook() {
    bool success;

    auto t0_noanalysis_simulation = std::chrono::steady_clock::now();
    if (!loaded) {
        CP(to_cleng);
        if (pivot_move) { // star HACK
            ids_node4fix.clear();
            ids_node4fix.push_back(pivot_arm_nodes[1][0]);
        }
    }
    // checks some node
    if (!Checks(0)) {cout << internal_name << "Some checks are not passed. Termination..." << endl; exit(1); }

    success = solveAndCheckFreeEnergy();
    if (!success) exit(1);

    auto t1_noanalysis_simulation = std::chrono::steady_clock::now();
    tpure_simulation = std::chrono::duration_cast<std::chrono::seconds> (t1_noanalysis_simulation - t0_noanalysis_simulation).count();

    return success;
}

void Cleng::notification() {
    cout << internal_name << "System for calculation: " << endl;
    for (auto &&n : nodes_map) cout << "Node id: " << n.first << " " << n.second.data()->get()->to_string() << endl;
    cout << endl;
    if (pivot_move) {
        cout << internal_name << "System [pivot -> stars only]: " << endl;
        for (auto &&pair_pivot : pivot_arm_nodes) {
            cout << "--> arm: " << pair_pivot.first << " ids: ";
            for (auto &&ids: pair_pivot.second) cout << ids << " ";
            cout << endl;
        }
        cout << endl;
    }
}

void Cleng::save(int num, Analyzer& analyzer) {

    cout << internal_name << "Saving... " << num << endl;

    if ((int(num) % delta_save) == 0) WriteOutput(num);
    if (cleng_dis) WriteClampedNodeDistanceWithCenter(num);
    if (cleng_dis) WriteClampedNodeDistance(num);
    if (cleng_pos) WriteClampedNodePosition(num);

    // STARS
    // vector<Real> vtk;
    // vtk.clear();
    // vtk = prepare_vtk("mol", "pol", "phi");
    // analyzer.updateVtk2PointsRepresentation(vtk, box);
    // // Re
    // Real Re_value = Analyzer::calculateRe(pivot_arm_nodes, nodes_map);
    // // Rg
    // Real Rg2_value = analyzer.calculateRg();
    // // phi
    // Real nr_check_sum = 0, phi_check_sum = 0;
    // map<int, vector<Real>> phi = analyzer.calculatePhi();
    // for (auto const& pair : phi) {
    //     //cout << "[phi] Layer: "<< pair.first << " | value[nr, phi]:";
    //     nr_check_sum  += pair.second[0];
    //     phi_check_sum += pair.second[1];
    //     //for (auto const& value : pair.second) {
    //     //    cout << value << " ";
    //     //}
    //     //cout << endl;
    // }
    // cout << "Total [nr]: " << nr_check_sum << " | [phi]:" << phi_check_sum << endl;

    // Write2File(num, "re",  Re_value);
    // Write2File(num, "rg2", Rg2_value);

    // vector<Real> phi_vector; vector<Real> nr_vector;
    // phi_vector.clear(); nr_vector.clear();
    // for (auto const& pair : phi) {nr_vector.push_back(pair.second[0]); phi_vector.push_back(pair.second[1]);}
    // Write2File(num, "phi", phi_vector);
    // Write2File(num, "nr",  nr_vector);

#ifdef CLENG_EXPERIMENTAL
    save2h5vtk();
    // free energy calculation
    n_times_mu = GetN_times_mu();
    vector<Real> MC_free_energy = {static_cast<Real>(num), free_energy_current, free_energy_current-n_times_mu};
    save2h5("free_energy", dims_3, MC_free_energy);
    // Not working yet.  TODO: FIX IT
    //// ReRg2
    //vector<Real> MC_ReRg2 = {static_cast<Real>(MC_attempt+MCS_checkpoint), Re_value, Rg2_value};
    //save2h5("ReRg2", dims_3, MC_ReRg2);
    //// phi
    //vector<Real> phi_vector; vector<Real> nr_vector;
    //for (auto const& pair : phi) {nr_vector.push_back(pair.second[0]); phi_vector.push_back(pair.second[1]);}
    //save2h5("phi/phi"+to_string(MC_attempt+MCS_checkpoint), dims_phi, phi_vector);
    //save2h5("nr/nr"+to_string(MC_attempt+MCS_checkpoint), dims_phi, nr_vector);
#endif

}

void Cleng::_save_differences(int mcs_, Real SCF_free_energy_trial, Real F_proposed) const {
//    Real F_proposed = getFreeEnergyBox();
    vector<Real> mc_energy_vector;
    mc_energy_vector.clear();
    mc_energy_vector.push_back(SCF_free_energy_trial);
    mc_energy_vector.push_back(F_proposed);
    Write2File(static_cast<int>(MC_attempt+MCS_checkpoint+mcs_), "ESCF_Ebox", mc_energy_vector, true);
}

#ifdef CLENG_EXPERIMENTAL
    void Cleng::save2h5vtk() {
        // vtk profiles per line of h5
        vector<Real> vtk;
        for (auto && line : out_per_line) {
            vtk.clear();
            vtk = prepare_vtk(line[0], line[1], line[2]);
            if (vtk.size() == (unsigned int)dims_vtk[0]) cleng_writer.write("/VTK_data", "vtk_"+line[0]+"|"+line[1]+"|"+line[2]+to_string(MC_attempt+MCS_checkpoint), dims_vtk, vtk);
            else cout << "vtk file was not saved because 'profile' was not found for " << line[0] << "|" << line[1] << "|" << line[2] << endl;
        }
    }

    void Cleng::save2h5(string what, vector<int> dims, vector<Real> value) {
        cleng_writer.append("/system_info/", what, dims, value);
    }

#endif
