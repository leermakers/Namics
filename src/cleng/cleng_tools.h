#include "cleng.h"

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

vector<shared_ptr<Node>> createNodes(const vector<shared_ptr<SimpleNode>> &simple_nodes) {
    vector<shared_ptr<Node>> result;
    map<SimpleNode, vector<shared_ptr<SimpleNode>>> m;
    for (auto &&n  : simple_nodes) { m[*n].push_back(n); }
    for (auto &&entry : m) {
        if (entry.second.size() == 1) result.push_back(entry.second[0]);
        else result.push_back(make_shared<Monolit>(entry.second));
    }

    return result;
}

vector<int> makeExcludedArray(int step) {
    vector<int> result;
    for (int i = -step + 1; i < step; i++) result.push_back(i);
    return result;
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
    bool success = true;
    vector<int> id_nodes;
    Point sub_box = {sub_box_size-2, sub_box_size-2, sub_box_size-2};  // change it in future

    if (!nodes[id_node_for_move]->inSubBoxRange(sub_box, clamped_move)) {
        id_nodes.push_back(id_node_for_move); success = false;
    }
    if (!id_nodes.empty()) {
        cout << "These nodes[id] are not in sub-box range:" << endl;
        for (auto &&nid: id_nodes) cout << nid << endl;
    }
    return success;
}

bool Cleng::NotCollapsing() {
    bool not_collapsing = true;

    Point P_not_shifted (nodes[id_node_for_move]->point());
    Point P_shifed(nodes[id_node_for_move]->point() + clamped_move);

    double min_dist = 0; // minimal distance between nodes
    for (auto &&n : nodes) {
        Point P_test = n->point();

        if (P_not_shifted != P_test) {
            Real distance = P_shifed.distance(n->point());
            if (distance <= min_dist) {
                cout << "Nodes are too close to each other." << endl;
                cout << "Shifted nodes: " << P_shifed.to_string() << endl;
                cout << "Nodes id:      " << P_test.to_string() << endl;
                not_collapsing =false;
            }
        }
    }
    return not_collapsing;
}

bool Cleng::InRange() {
    bool in_range = false;
    Point box = {Lat[0]->MX, Lat[0]->MY, Lat[0]->MZ};
    Point down_boundary = {1, 1, 1};
    Point up_boundary = box - down_boundary;
    Point MPs(nodes[id_node_for_move]->point() + clamped_move);

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

Point Cleng::prepareMove() {
    if (debug) cout << "prepareMove in Cleng" << endl;

    if (axis) {
        int c1=1; if (sign_move == "-") c1=-1;
        // currently it is implemented for step
        if (axis == 1) clamped_move = {c1*2, 0, 0};
        if (axis == 2) clamped_move = {0, c1*2, 0};
        if (axis == 3) clamped_move = {0, 0, c1*2};
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
        while (infile >> ws && getline(infile, line)); // skip empty lines

        istringstream iss(line);
        vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());
        MS_step = stoi(results[0]);
    } else cout << "Unable to open kal file.\n";

    return MS_step;
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