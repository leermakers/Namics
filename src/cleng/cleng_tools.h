#include "cleng.h"
#include "iterator/EnumerateIterator.h"

volatile sig_atomic_t cleng_flag_termination = 0;

using Fp_info = numeric_limits<double>;
inline auto is_ieee754_nan( double const x ) -> bool {
    static constexpr bool   is_claimed_ieee754  = Fp_info::is_iec559;
    static constexpr int    n_bits_per_byte     = CHAR_BIT;
    using Byte = unsigned char;

    static_assert( is_claimed_ieee754, "!" );
    static_assert( n_bits_per_byte == 8, "!" );
    static_assert( sizeof( x ) == sizeof( uint64_t ), "!" );

#ifdef _MSC_VER
    uint64_t const bits = reinterpret_cast<uint64_t const&>( x );
#else
    Byte bytes[sizeof(x)];
    memcpy( bytes, &x, sizeof( x ) );
    uint64_t int_value;
    memcpy( &int_value, bytes, sizeof( x ) );
    uint64_t const& bits = int_value;
#endif

    static constexpr uint64_t   sign_mask       = 0x8000000000000000;
    static constexpr uint64_t   exp_mask        = 0x7FF0000000000000;
    static constexpr uint64_t   mantissa_mask   = 0x000FFFFFFFFFFFFF;

    (void) sign_mask;
    return (bits & exp_mask) == exp_mask and (bits & mantissa_mask) != 0;
}

shared_ptr<SimpleNode> fromSystemToNode(int x, int y, int z, int id, const Point &box) {
    return make_shared<SimpleNode>(Point(x, y, z), id, box);
}

map<int, vector<shared_ptr<Node>>> createNodes(const vector<shared_ptr<SimpleNode>> &simple_nodes) {
    map<int, vector<shared_ptr<Node>>> result_map;
    int index=0;
    map<SimpleNode, vector<shared_ptr<SimpleNode>>> m;
    for (auto &&n  : simple_nodes) { m[*n].push_back(n); }
    for (auto &&entry : m) {
        if (entry.second.size() == 1) result_map[index] = {entry.second[0]};
        else result_map[index] = {make_shared<Monolit>(entry.second)};
        index ++;
    }
    return result_map;
}

vector<int> makeExcludedArray(int step) {
    vector<int> result;
    for (int i = -step + 1; i < step; i++) result.push_back(i);
    return result;
}

bool Cleng::InSubBoxRange() {
    bool success = true;
    vector<int> id_nodes;
    Point sub_box = {sub_box_size.x - 2, sub_box_size.y - 2, sub_box_size.z - 2};  // change it in future

    if (!nodes_map[id_node_for_move].data()->get()->inSubBoxRange(sub_box, clamped_move)) {
        id_nodes.push_back(id_node_for_move); success = false;
    }
    if (!id_nodes.empty()) {
        cout << "These nodes_map[id] are not in sub-box range:" << endl;
        for (auto &&nid: id_nodes) cout << nid << endl;
    }
    return success;
}

bool Cleng::IsCommensuratable() {
    // TODO name think about changing the name.
    bool success = true;
    // trial movement for simpleNodeList
    nodes_map[id_node_for_move].data()->get()->shift(clamped_move);

    int length = (int) In[0]->MolList.size();

    for (int i=0; i<length;i++) {
    int chain_length = Mol[i]->chainlength-2;
    // TODO it will not work for multicomponent systems (gel + salt for example...)!
    // TODO to next commit
    if (chain_length < 10) continue;


    for (auto &&SN : Enumerate(simpleNodeList)) {
        auto p1 = SN.second->get_system_point();
        auto p2 = SN.second->get_cnode()->get_system_point();
        auto p3 = p2 - p1;
        int path_length = abs(p3.x) + abs(p3.y) + abs(p3.z);

        int path_length_even = path_length % 2;
        int chain_length_even = chain_length % 2;

        // cout << "path: " << path_length << endl;
        // cout << "chain: " << chain_length << endl;

        if (path_length_even == chain_length_even) success =false;
        if (path_length >= chain_length) success =false;
    }
    }
    if (!success) cout << "Warning, the paths between clamps is not commensurate with the chain length!" << endl;
    // put back
    nodes_map[id_node_for_move].data()->get()->shift(clamped_move.negate());
    return success;
}

bool Cleng::NotCollapsing() {
    bool not_collapsing = true;

    Point P_not_shifted (nodes_map[id_node_for_move].data()->get()->point());
    Point P_shifted(nodes_map[id_node_for_move].data()->get()->point() + clamped_move);

    double min_dist = 0; // minimal distance between nodes_map. 0 for a while...
    for (auto &&n : nodes_map) {
        Point P_test = n.second.data()->get()->point();

        if (P_not_shifted != P_test) {
            Real distance = P_shifted.distance(n.second.data()->get()->point());
            if (distance <= min_dist) {
                cout << "Nodes are too close to each other."           << endl;
                cout << "Shifted nodes_map: " << P_shifted.to_string() << endl;
                cout << "Nodes id:          " << P_test.to_string()    << endl;
                not_collapsing =false;
            }
        }
    }
    return not_collapsing;
}

void signalHandler(int signum) {
    cout << "Termination..." << endl;
    cleng_flag_termination++;
    if (cleng_flag_termination > 1) exit(signum);
}

bool Cleng::InRange() {
    bool in_range = false;
    Point down_boundary = {1, 1, 1};
    Point up_boundary = box - down_boundary;
    Point MPs(nodes_map[id_node_for_move].data()->get()->point() + clamped_move);

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

vector<Real> Cleng::prepare_vtk() {
    int Size = 0;
    vector<Real> vtk;
    Real *X = Out[0]->GetPointer(Out[0]->OUT_key[0], Out[0]->OUT_name[0], Out[0]->OUT_prop[0], Size);
    for (int i = 1; i < box.x + 1; i++) {
        for (int j = 1; j < box.y + 1; j++) {
            for (int k = 1; k < box.z + 1; k++) {
                vtk.push_back(X[i * J.x + j * J.y + k]);
            }
        }
    }
    return vtk;
}

Point Cleng::prepareMove() {
    if (debug) cout << "prepareMove in Cleng" << endl;

    if (axis) {
        int c1=1;  // sing of movement +1 or -1
        int c2=2;  // MC step by default 2; could be 1 in two_ends_extension mode

        if (two_ends_extension) c2 = 1;
        if (sign_move == "-") c1=-1;
        // currently it is implemented for step
        if (axis == 1) clamped_move = {c1*c2, 0, 0};
        if (axis == 2) clamped_move = {0, c1*c2, 0};
        if (axis == 3) clamped_move = {0, 0, c1*c2};
    } else {

        clamped_move = {
                rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
                rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
                rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
        };

        if ((delta_step % 2) == 1) {  // 1 3 5
            int id4rm = rand.getInt(0, 2);
            if (id4rm == 0) clamped_move.x = 0;
            else { if (id4rm == 1) clamped_move.y = 0; else clamped_move.z = 0;}
        } else {  // 2 4 6
            int id4rm1 = rand.getInt(0, 2);
            int id4rm2 = rand.getInt(0, 2);

            if (id4rm1 == 0) clamped_move.x = 0;
            else {if (id4rm1 == 1) clamped_move.y = 0; else clamped_move.z = 0;}

            if (id4rm2 == 0) clamped_move.x = 0;
            else { if (id4rm2 == 1) clamped_move.y = 0; else clamped_move.z = 0;}
        }
    }
    return clamped_move;
}

int Cleng::getLastMCS() {

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

void Cleng::WriteOutput() {
    if (debug) cout << "WriteOutput in Cleng" << endl;
    PushOutput();
    New[0]->PushOutput();
    for (int i = 0; i < n_out; i++) Out[i]->WriteOutput(MC_attempt + MCS_checkpoint);
}

void Cleng::PushOutput() {

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
            Out[i]->push("MC_attempt", MC_attempt + MCS_checkpoint);
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

void Cleng::WriteClampedNodeDistance() {
    vector<Real> distPerMC;
    int i = 0;
    for (auto &&SN :  simpleNodeList) {
        if (!(i % 2)) distPerMC.push_back(SN->distance(SN->get_cnode()->get_system_point()));
        i++;
    }

    ofstream outfile;
    outfile.open(filename + ".cdis", std::ios_base::app);
    outfile << MC_attempt + MCS_checkpoint << " ";
    for (auto n : distPerMC)outfile << n << " ";
    outfile << endl;
}

void Cleng::WriteClampedNodePosition() {

    ofstream outfile;
    outfile.open(filename + ".cpos", std::ios_base::app);

//    cout
    outfile
    << "#step " << MC_attempt + MCS_checkpoint << " {X, Y, Z} #" << endl;

    int i = 0;
    for (auto &&SN :  simpleNodeList) {
        if (!(i % 2)) {
//            cout
            outfile
//            << SN->get_system_point().to_string() << SN->get_cnode()->get_system_point().to_string() << endl;
            << SN->get_system_point().to_string() << endl;
        }
        i++;
    }
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


