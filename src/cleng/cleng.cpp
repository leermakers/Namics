#include "cleng.h"
#include "point.h"
#include "monolit.h"
#include <map>
#include <cassert>

using namespace std;

Cleng::Cleng(
        vector<Input*> In_,
        vector<Lattice*> Lat_,
        vector<Segment*> Seg_,
        vector<Molecule*> Mol_,
        vector<System*> Sys_,
        vector<Solve_scf*> New_,
        string name_
) : name(name_),
    In(In_),
    Lat(Lat_),
    Mol(Mol_),
    Seg(Seg_),
    Sys(Sys_),
    New(New_)

{
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
            success = In[0]->Get_int(GetValue("MCS"), MCS, 1, 10000, "The number of timesteps should be between 1 and 10000");
        }
        if (debug) cout << "MCS is " << MCS << endl;
        
        if (GetValue("save_interval").size() > 0) {
            success = In[0]->Get_int(GetValue("save_interval"), save_interval,1,MCS,"The save interval nr should be between 1 and 100");
        }
        if (debug) cout << "Save_interval " << save_interval << endl;
        
        if (Sys[0]->SysClampList.size() <1) {
            cout <<"Cleng needs to have clamped molecules in the system" << endl; success=false;
        }
        else {
            clamp_seg=Sys[0]->SysClampList[0]; 
            if (Sys[0]->SysClampList.size()>1) {
                success=false; cout <<"Currently the clamping is limited to one molecule per system. " << endl; 
            }
        }
        
        if (success) {
            n_boxes = Seg[clamp_seg]->n_box;
            sub_box_size=Seg[clamp_seg]->mx;
        }
        clp_mol=-1;
        int length = In[0]->MolList.size();
        for (int i=0; i<length; i++) if (Mol[i]->freedom =="clamped") clp_mol=i; 
    }
    if (success) {
        n_out = In[0]->OutputList.size();
        if (n_out == 0) cout << "Warning: no output defined!" << endl;
        
        for (int i =  0; i < n_out; i++) {
            Out.push_back(new Output(In, Lat, Seg, Mol, Sys, New, In[0]->OutputList[i], i, n_out));
            if (!Out[i]->CheckInput(start)) {
                cout << "input_error in output " << endl;
                success=false;
            }
        }
        MonteCarlo();
    }
    return success;
}

SimpleNode fromSystemToNode(int x, int y, int z, int id, const Point& box) {
    return SimpleNode({x, y, z}, id, box);
}

vector<shared_ptr<Node>> createNodes(const vector<SimpleNode>& simple_nodes) {
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

    bool success=true;
    Point box {Lat[0]->MX, Lat[0]->MY, Lat[0]->MZ};
	
    int JX=Lat[0]->JX;
	int JY=Lat[0]->JY;
	int M =Lat[0]->M;

    // Not used
	//int j;
    //int length;
	//bool found;

    vector<SimpleNode> sn;
    Segment* clamped = Seg[clamp_seg];
    map<int, Point> system_points;
	// TODO think about simplify
	switch(tofrom) {
        case to_cleng:
            for (int i = 0; i < n_boxes; i++) {
                sn.push_back(fromSystemToNode(clamped->px1[i], clamped->py1[i], clamped->pz1[i], 2 * i, box));
                sn.push_back(fromSystemToNode(clamped->px2[i], clamped->py2[i], clamped->pz2[i], 2 * i + 1, box));
            }

            nodes = createNodes(sn);
            for (auto && n : nodes) {
                cout << n->to_string() << endl;
            }
//            assert(nodes.size() == 3);
            break;

        case to_segment:
			Zero(Seg[clamp_seg]->H_MASK,M);

            for (auto &&n : nodes) {
                n->pushSystemPoints(system_points);
            }

            for (auto &&entry : system_points) {
                int i = entry.first / 2;
                Point& p = entry.second;
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

				if (Seg[clamp_seg]->bx[i] < 1) {Seg[clamp_seg]->bx[i] += box.x; Seg[clamp_seg]->px1[i] +=box.x; Seg[clamp_seg]->px2[i] +=box.x;}
				if (Seg[clamp_seg]->by[i] < 1) {Seg[clamp_seg]->by[i] += box.y; Seg[clamp_seg]->py1[i] +=box.y; Seg[clamp_seg]->py2[i] +=box.y;}
				if (Seg[clamp_seg]->bz[i] < 1) {Seg[clamp_seg]->bz[i] += box.z; Seg[clamp_seg]->pz1[i] +=box.z; Seg[clamp_seg]->pz2[i] +=box.z;}

				Seg[clamp_seg]->H_MASK[((Seg[clamp_seg]->px1[i]-1)%box.x+1)*JX + ((Seg[clamp_seg]->py1[i]-1)%box.y+1)*JY + (Seg[clamp_seg]->pz1[i]-1)%box.z+1]=1;
				Seg[clamp_seg]->H_MASK[((Seg[clamp_seg]->px2[i]-1)%box.x+1)*JX + ((Seg[clamp_seg]->py2[i]-1)%box.y+1)*JY + (Seg[clamp_seg]->pz2[i]-1)%box.z+1]=1;

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
			success=false; 
			cout <<"error in transfer" << endl;
		break;
	}
	return success;
}

void Cleng::WriteOutput(int subloop){
    if (debug) cout << "WriteOutput in Cleng" << endl;
    int mon_length;
    int mol_length;
    int molal_length;

    PushOutput();
    Sys[0]->PushOutput(); // needs to be after pushing output for seg.
    Lat[0]->PushOutput();
    New[0]->PushOutput();
    mon_length = (int)In[0]->MonList.size();
    for (int i = 0; i < mon_length; i++) {
        Seg[i]->PushOutput();
    }
    mol_length = (int)In[0]->MolList.size();
    for (int i = 0; i < mol_length; i++) {
        molal_length = (int) Mol[i]->MolAlList.size();
        for (int k = 0; k < molal_length; k++) {
            Mol[i]->Al[k]->PushOutput();
        }
        Mol[i]->PushOutput();
    }
    // length = In[0]->AliasList.size();
    // for (int i=0; i<length; i++) Al[i]->PushOutput();

    for (int i = 0; i < n_out; i++) {
        Out[i]->WriteOutput(subloop);
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

bool Cleng::InBoxRange() {
    int up_bondary = 3; int down_bondary = 3;
    // TODO simplify this
    return    (down_bondary < nodes[rand_part_index]->point().x + shift.x) && (nodes[rand_part_index]->point().x + shift.x < (int)Lat[0]->MX - up_bondary)
           && (down_bondary < nodes[rand_part_index]->point().y + shift.y) && (nodes[rand_part_index]->point().y + shift.y < (int)Lat[0]->MY - up_bondary)
           && (down_bondary < nodes[rand_part_index]->point().z + shift.z) && (nodes[rand_part_index]->point().z + shift.z < (int)Lat[0]->MZ - up_bondary);
}

bool Cleng::NotTooClose() {
    
    Point MP{nodes[rand_part_index]->point().x + shift.x, nodes[rand_part_index]->point().y+ shift.y, nodes[rand_part_index]->point().z+ shift.z};

    for (auto n : nodes) {
        //TODO is this equal refactoring (maybe will be problems with monolit)?
        if (MP == n->point()) {
            return false;
        }
    }
    return true;
}

bool Cleng::MakeShift1(bool back) {
    if (debug) cout << "MakeShift1 in Cleng" << endl;
    bool success = true;

    cout << "clamped seg:" << clamp_seg << endl;
    for (int i=0; i < n_boxes; i++) {
        cout << "clemped seg pos_x1: " << Seg[clamp_seg]-> px1[i] << " " << Seg[clamp_seg]-> py1[i] << " " << Seg[clamp_seg]-> pz1[i] << endl;
        cout << "clemped seg pos_x2: " << Seg[clamp_seg]-> px2[i] << " " << Seg[clamp_seg]-> py2[i] << " " << Seg[clamp_seg]-> pz2[i] << endl;
    }

    cout << "len X:" << nodes.size() << endl;
    cout << "len Y:" << nodes.size() << endl;
    cout << "len Z:" << nodes.size() << endl;

    for (int i=0; i < (int)nodes.size(); i++) {
        cout << "i:" << i << " X:" << nodes[i]->point().x << " Y:" << nodes[i]->point().y << " Z:" << nodes[i]->point().z << endl;
    }

    return success;
}

bool Cleng::MakeShift(bool back) {
    if (debug) cout << "MakeShift in Cleng" << endl;
    bool success = true;

    if (!back) {
        shift = {0, 0, 0};
        rand_part_index = GetIntRandomValueExclude(0, (int)nodes.size() - 1, 0, false);

        cout << "rand_part_index:" << rand_part_index << endl;
        cout << "X:" << nodes[rand_part_index]->point().x << endl;
        cout << "Y:" << nodes[rand_part_index]->point().y << endl;
        cout << "Z:" << nodes[rand_part_index]->point().z << endl;

//    TODO: rethink physics
//    pos_array = GetIntRandomValueExclude(0, 2, 0, false);
//    cout << "pos_array:" << pos_array << endl;
//    //choosing what direction should I change (x,y or z)
//    for (int i=0; i<3; i++) {
//        shift[i] = GetIntRandomValueExclude(-1, 1, 0, true);
//    }
//    shift[pos_array] = 0;

        int pos_array = GetIntRandomValueExclude(0, 1, 0, false);
        cout << "pos_array:" << pos_array << endl;

        if (pos_array == 0) {
//      z-direction
            shift = {shift.x, shift.y, GetIntRandomValueExclude(-1, 1, 0, true)};
        } else {
//      xy-direction
            shift = {GetIntRandomValueExclude(-1, 1, 0, true), GetIntRandomValueExclude(-1, 1, 0, true), shift.z};
        }

        if (InBoxRange() && NotTooClose()) {
            nodes[rand_part_index]->shift(shift);

        } else {
            if (pos_array == 0) {
//      z-direction
                shift = {shift.x, shift.y, GetIntRandomValueExclude(-1, 1, 0, true)};
            } else {
//      xy-direction
                shift = {GetIntRandomValueExclude(-1, 1, 0, true), GetIntRandomValueExclude(-1, 1, 0, true), shift.z};
            }
        }

        cout << "len X:" << nodes.size() << endl;
        cout << "len Y:" << nodes.size() << endl;
        cout << "len Z:" << nodes.size() << endl;


//    cout << "changed:" << changed << endl;
        cout << "Shift:" << shift.x << " " << shift.y << " " << shift.z << endl;
//    cout << "len P:" << P.size() << endl;
//    cout << "len Sx:" << Sx.size() << endl;
//    cout << "len Sy:" << Sy.size() << endl;
//    cout << "len Sz:" << Sz.size() << endl;
//    cout << "len X:" << nodes.size() << endl;
//    cout << "len Y:" << nodes.size() << endl;
//    cout << "len Z:" << nodes.size() << endl;
    cout << endl;
    }
    else {
        cout << "MakeShift back" << endl;

        cout << "rand_part_index:" << rand_part_index << endl;
        cout << "X:" << nodes[rand_part_index]->point().x << endl;
        cout << "Y:" << nodes[rand_part_index]->point().y << endl;
        cout << "Z:" << nodes[rand_part_index]->point().z << endl;

        cout << -shift.x << " " << -shift.y << " " << -shift.z << endl;
        cout << "last particle:" << rand_part_index << endl;
        // TODO simplify
        Point neg_shift = shift.negate();
        nodes[rand_part_index]->shift(neg_shift);

        cout << "rand_part_index:" << rand_part_index << endl;
        cout << "X:" << nodes[rand_part_index]->point().x << endl;
        cout << "Y:" << nodes[rand_part_index]->point().y << endl;
        cout << "Z:" << nodes[rand_part_index]->point().z << endl;

        cout << endl;

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
    for (int i = 1; i < MCS; i++) { // loop for trials

        //cout << "Step i:" << i << endl;

        //for (int i = 0; i < n_boxes; i++) {
        //    cout << "Subbox_x" << Seg[clamp_seg]->bx[i] << " " << Seg[clamp_seg]->by[i] << " " << Seg[clamp_seg]->bz[i] << endl;
        //}
//            sn.push_back(fromSystemToNode(clamped->px1[i], clamped->py1[i], clamped->pz1[i], 2 * i, box));
//            sn.push_back(fromSystemToNode(clamped->px2[i], clamped->py2[i], clamped->pz2[i], 2 * i + 1, box));
//        }

        Real my_rand = GetRealRandomValue(0, 1);
        free_energy_c = Sys[0]-> FreeEnergy;
        success=CP(to_cleng);
        MakeShift(false);
        success=CP(to_segment);
        New[0]->Solve(true);

//        for (int i=0; i < n_boxes; i++) {

//            cout << "################ START" << endl;
//            cout << "clemped seg pos_x1: " << Seg[clamp_seg]-> px1[i] << " pos_y1:" << Seg[clamp_seg]-> py1[i] << " pos_z1: " << Seg[clamp_seg]-> pz1[i] << endl;
//            cout << "clemped seg pos_x2: " << Seg[clamp_seg]-> px2[i] << " pos_y2:" << Seg[clamp_seg]-> py2[i]  << " pos_z2: " << Seg[clamp_seg]-> pz2[i] << endl;
//        }
        free_energy_t = Sys[0]-> FreeEnergy;

        cout << "my_rand:" << my_rand << endl;
        cout << "free_energy_c:" << free_energy_c << endl;
        cout << "free_energy_t:" << free_energy_t << endl;

        if (std::isnan(free_energy_t)) {
            for (int k = 0; k < (int) nodes.size(); k++) {
                cout << nodes[k]->point().x << " " << nodes[k]->point().y << " " << nodes[k]->point().z << endl;
            }
            break;
        }


        if ( my_rand < exp(free_energy_c-free_energy_t) ) {
            cout << "Accepted" << endl;
        } else {
            cout << "Deny" << endl;
            MakeShift(true);
            continue;
        }
        if ((i % save_interval) == 0) {
            WriteOutput(i);
        }
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

void Cleng::PushOutput() {
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
		if (Out[i]->name=="ana" || Out[i]->name=="kal") Out[i]->push("t",t);
		if (Out[i]->name=="ana" || Out[i]->name=="kal") Out[i]->push("MCS",MCS);
		if (Out[i]->name=="ana" || Out[i]->name=="vec") { //example for putting an array of Reals of arbitrary length to output
			string s="vector;0"; //the keyword 'vector' is used for Reals; the value 0 is the first vector, use 1 for the next etc, 
			Out[i]->push("gn",s); //"gn" is the name that will appear in output file
			Out[i]->PointerVectorReal.push_back(Mol[clp_mol]->gn); //this is the pointer to the start of the 'vector' that is reported to output.
			Out[i]->SizeVectorReal.push_back(sizeof(Mol[clp_mol]->gn)); //this is the size of the 'vector' that is reported to output
		}
		if (Out[i]->name=="ana" || Out[i]->name=="pos") { //example for putting 3 array's of Integers of arbitrary length to output
			string s="array;0";
			Out[i]->push("X",s);

//			 TODO using normal point in output and remove xs, ys, zs
			fillXYZ();
//			point=X.data();
			Out[i]->PointerVectorInt.push_back(xs);
			Out[i]->SizeVectorInt.push_back(nodes.size());
			s="array;1";
			Out[i]->push("Y",s);
//			point = Y.data();
			Out[i]->PointerVectorInt.push_back(ys);
			Out[i]->SizeVectorInt.push_back(nodes.size());
			s="array;2";
			Out[i]->push("Z",s);
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