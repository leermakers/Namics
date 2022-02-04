#include "cleng.h"
#include "cleng_tools.h"
#include "cleng_moves.h"

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
    KEYS.emplace_back("cleng_dis");
    KEYS.emplace_back("simultaneously");
    KEYS.emplace_back("movement_along");
    KEYS.emplace_back("sign_move");
    KEYS.emplace_back("user_node_id_move");
    KEYS.emplace_back("two_ends_extension");
    KEYS.emplace_back("metropolis");
    KEYS.emplace_back("prefactor_kT");
    KEYS.emplace_back("one_node");
    KEYS.emplace_back("pivot_move");
    KEYS.emplace_back("pivot_axis");
    KEYS.emplace_back("pivot+one_node");
    KEYS.emplace_back("pivot+one_bond");
    KEYS.emplace_back("h5");
    KEYS.emplace_back("warming_up_steps");
    KEYS.emplace_back("warming_up_stage");
    KEYS.emplace_back("ESCF_Ebox");
    //
    KEYS.emplace_back("mc_engine");
    //
    KEYS.emplace_back("inner_loop_each");
    KEYS.emplace_back("mcs");
    KEYS.emplace_back("start_inner_loop_from");

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

//bool Cleng::CheckInput(int start, bool save_vector) { // testing
bool Cleng::CheckInput(int start) {
    if (debug) cout << "CheckInput in Cleng" << endl;
    bool success;

    success = In[0]->CheckParameters("cleng", name, start, KEYS, PARAMETERS, VALUES);
    if (success) {

        // MCS
        if (!GetValue("MCS").empty()) {
            success = In[0]->Get_int(GetValue("MCS"), MCS, 0, 10000000, "The number of Monte Carlo steps should be between 0 and 10_000_000; MCS = 0 is special case for calculating using only SCF part and Cleng tools.");
            if (!success) { MCS = 1; cout << "MCS will be equal to " << MCS << endl; }
        } else MCS = 0;
        if (debug) cout << "MCS is " << MCS << endl;
        
        // mcs
        if (!GetValue("mcs").empty()) {
            success = In[0]->Get_int(GetValue("mcs"), mcs, 1, 10000000, "The number of inner Monte Carlo steps should be between 1 and 10_000_000;");
            if (!success) { mcs = 10; cout << "Inner mcs will be equal to " << mcs << endl; }
        } else mcs = 0;
        if ( (MCS == 0) && (mcs != 0) ) {
            cout << "[WARNING] MCS equals to zero, however inner loop is not zero." << endl;
            cout << "   -->>   inner loop will be set to zero too." << endl;
            cout << "   -->>   use MCS value instead." << endl;
            mcs = 0;
        }
        if (debug) cout << "mcs is " << mcs << endl;

        // start_inner_loop_from SILF
        if (!GetValue("start_inner_loop_from").empty()) {
            success = In[0]->Get_int(GetValue("start_inner_loop_from"), SILF, 1, MCS, "Starting point for inner loop should be between 1 and MCS flag.");
            if (!success) { SILF = 1; cout << "Starting point for inner loop will be equal to " << SILF << endl; }
        } else SILF = MCS+5;
        if ( (mcs != 0 ) && (SILF == 0) ) {
            SILF = 1;
            cout << "[WARNING] You should provided 'start_inner_loop_from' flag for inner MC loop [mcs flag]." << endl;
            cout << "   -->>   'start_inner_loop_from' will be equal to " << SILF << endl;
        }
        if (debug) cout << "SILF is " << SILF << endl;

        // inner_loop_each ILE
        if (!GetValue("inner_loop_each").empty()) {
            success = In[0]->Get_int(GetValue("inner_loop_each"), ILE, 1, MCS, "Entrance to inner loop should be between 1 and MCS flag;");
            if (!success) { ILE = 1; cout << "ILE will be equal to " << ILE << endl; }
        } else ILE = 0;
        if ( (ILE == 0) && (mcs != 0) ) {
            ILE = 1;
        }
        if (debug) cout << "ILE is " << ILE << endl;

        // seed
        if (!GetValue("seed").empty()) {
            success = In[0]->Get_int(GetValue("seed"), pseed, 1, 1000, "The seed should be between 1 and 1000");
            if (!success) { pseed = 1; cout << "The seed will be equal to " << pseed << endl; }
        } else pseed = 0;
        rand = pseed==0 ? Random() : Random(pseed);
        if (debug) cout << "seed is " << pseed << endl;

        // delta_step
        if (!GetValue("delta_step").empty()) {
            success = In[0]->Get_int(GetValue("delta_step"), delta_step, 1, 5, "The number of delta_step should be between 1 and 5");
            if (!success) { delta_step = 1; cout << "The delta_step will be equal to " << delta_step << endl; }
        } else delta_step = 0;
        if (debug) cout << "delta_step is " << delta_step << endl;
        //

        // pivot_move
        if (!GetValue("pivot_move").empty()) {
            success = In[0]->Get_int(GetValue("pivot_move"), pivot_move, 1, 360, "The angle of pivot_move should be between 1 and 360");
            if (!success) { cout << "The pivot_move will be disable." << endl; pivot_move = 0;}
        } else pivot_move = 0;
        if (debug) cout << "pivot_move is " << pivot_move << endl;

        // pivot_axis
        if (pivot_move) {
            if (!GetValue("pivot_axis").empty()) {
                success = In[0]->Get_int(GetValue("pivot_axis"), pivot_axis, 1, 3,
                                         "The axis of pivot_move should be between 1 and 3");
                if (!success) {
                    cout << "The pivot_axis will be:  all axis." << endl;
                    pivot_axis = -1;
                }
            } else pivot_axis = -1;
            if (debug) cout << "pivot_axis is " << pivot_axis << endl;
        } else {pivot_axis = 0;} // pivot_axis = 0 means disable pivot.

        // one_node HACK: clean up after
        if (!GetValue("one_node").empty()) {
            bool one_node_flag = In[0]->Get_bool(GetValue("one_node"), false);
            if (pivot_move) {
                if (one_node_flag) {
                    cout << "Sorry, if you like to use combination pivot+one_node move, please, specify pivot movement "
                            "by parameter 'pivot_move'"
                            "\n Termination..." << endl;
                    exit(0);
                }
            } else {
                if (one_node_flag) {
                    pivot_move = false;
                }
            }
        }
        if (debug) cout << "one_node " << pivot_one_node << endl;

        // pivot+one_node
        if (!GetValue("pivot+one_node").empty()) {
            pivot_one_node = In[0]->Get_bool(GetValue("pivot+one_node"), false);
        } else pivot_one_node = false;
        if (debug) cout << "pivot+one_node " << pivot_one_node << endl;
        if ((pivot_one_node) and (!pivot_move)) {
            cout << "Sorry, if you like to use combination pivot+one_node move, please, specify pivot movement "
                    "by parameter 'pivot_move'"
                    "\n Termination..." << endl;
            exit(0);
        }

        // pivot+one_bond
        if (!GetValue("pivot+one_bond").empty()) {
            pivot_one_bond = In[0]->Get_bool(GetValue("pivot+one_bond"), false);
        } else pivot_one_bond = false;
        if (debug) cout << "pivot+one_node " << pivot_one_bond << endl;
        if ((pivot_one_bond) and (!pivot_move)) {
            cout << "Sorry, if you like to use combination pivot+one_bond move, please, specify pivot movement "
                    "by parameter 'pivot_move'"
                    "\n Termination..." << endl;
            exit(0);
        }

        if (pivot_one_bond) pivot_one_node = false;
        // TODO: think about user
//        if ((delta_step) and (pivot_move) and ()) {
//            cout << "Sorry, but you have to choose either one_node_move by `delta_step` parameter or `pivot_move`"
//                    "\n Termination..." << endl;
//            exit(0);
//        }
//
        // DEBUG Energy SCF and box model 
        if (!GetValue("ESCF_Ebox").empty()) {
            cout << "[WARNING] It is experimental/debug feature. Please, be sure what you do." << endl;
            escf_ebox_flag = In[0]->Get_bool(GetValue("ESCF_Ebox"), false);
        }
        if (debug) cout << "one_node " << pivot_one_node << endl;

        // delta_save
        if (!GetValue("delta_save").empty()) {
            success = In[0]->Get_int(GetValue("delta_save"), delta_save, 1, MCS+1, "The delta_save interval should be between 1 and " + to_string(MCS+1));
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
        if (!GetValue("checkpoint_save").empty()) {checkpoint_save = In[0]->Get_bool(GetValue("checkpoint_save"), false);}
        else checkpoint_save = false;
        if (debug) cout << "checkpoint_save " << checkpoint_save << endl;

        // checkpoint load
        if (!GetValue("checkpoint_load").empty()) {checkpoint_load = In[0]->Get_bool(GetValue("checkpoint_load"), false);}
        else checkpoint_load = false;
        if (debug) cout << "checkpoint_load " << checkpoint_load << endl;

        // saving cleng position of nodes_map
        if (!GetValue("cleng_pos").empty()) {cleng_pos = In[0]->Get_bool(GetValue("cleng_pos"), false);}
        else cleng_pos = false;
        if (debug) cout << "cleng_pos " << cleng_pos << endl;

        // saving distance between of nodes_map
        if (!GetValue("cleng_dis").empty()) {cleng_dis = In[0]->Get_bool(GetValue("cleng_dis"), false);}
        else cleng_dis = false;
        if (debug) cout << "cleng_dis " << cleng_dis << endl;

        // simultaneously
        if (!GetValue("simultaneously").empty()) simultaneously = In[0]->Get_bool(GetValue("simultaneously"), false);
        else simultaneously = false;
        if (debug) cout << "simultaneously move " << simultaneously << endl;

        // movement_along
        if (!GetValue("movement_along").empty()) {
            cout << "Warning!!! In movement_along mode delta step will be ignored! Monte Carlo step will be {sign_move*2} depending on axis " << endl;
            success = In[0]->Get_int(GetValue("movement_along"), axis, 1, 3, "The number of delta_step should be between 1 and 3");
            if (!success) {
                cout << "Sorry, you provide incorrect axis number. The movement_along will be disable" << endl;
                axis = 0;
                success = true;
            }
        } else axis = 0;
        if (debug) cout << "movement_along axis " << axis << endl;

        // warming_up_steps  NOT IMPLEMENTED
        if (!GetValue("warming_up_steps").empty()) {
            success = In[0]->Get_int(GetValue("warming_up_steps"), warming_up_steps, 1, 1000, "The warming_up_steps should be between 1 and 1000");
            if (!success) { warming_up_steps = 10; cout << "The warming_up steps will be equal to " << warming_up_steps << endl; }
        } else warming_up_steps = 10;
        if (debug) cout << "warming_up_steps is " << warming_up_steps << endl;

        // warming_up stage  NOT IMPLEMENTED
        if (!GetValue("warming_up_stage").empty()) do_warming_up_stage = In[0]->Get_bool(GetValue("warming_up_stage"), false);
        else do_warming_up_stage = false;
        if (debug) cout << "do_warming_up_stage " << do_warming_up_stage << endl;

        // sign_move
        vector<string> options {"+", "-"};
        if (axis) {
            if (!GetValue("sign_move").empty()) {
                success = In[0]->Get_string(GetValue("sign_move"), sign_move, options, "The sigh could be ether + or -");
            } else { sign_move = "+";}
        } else {
            if (!GetValue("sign_move").empty()) {
                cout << "Sorry, but you cannot use sign_move without movement_along axis flag" << endl;
                success = false;
                return success;
            }
        }
        if (debug) cout << "sign_move " << sign_move << endl;

        // user_node_move_id
        string struser_node_move_id;
        if (!GetValue("user_node_id_move").empty()) {
            success = In[0]->Get_string(GetValue("user_node_id_move"), struser_node_move_id, "");
            string node_id;
            string _ ;
            stringstream stream(struser_node_move_id);
            while( getline(stream, node_id, ',') ) {
                if (node_id[0] == '~') {
                    _ = node_id.erase(0,1);
                    try { ids_node4fix.push_back(stoi(_)); }
                    catch (invalid_argument &e) {
                        success = false;
                        cout << "Check yours node id structure" << endl;
                        return success;
                    }
                } else {
                    try { ids_node4move.push_back(stoi(node_id)); }
                    catch (invalid_argument &e) {
                        success = false;
                        cout << "Check yours node id structure" << endl;
                        return success;
                    }
                }
            }
        }
        else ids_node4move = {-1};
        if (debug) for (auto &&id:ids_node4move) cout << "user_node_id_move: " << id << endl;
        if (debug) for (auto &&id:ids_node4fix)  cout << "user_node_id_fix:  "  << id << endl;

        // 2 end extension mode
        if (!GetValue("two_ends_extension").empty()) {two_ends_extension = In[0]->Get_bool(GetValue("two_ends_extension"), false);}
        else two_ends_extension = false;
        if (debug) cout << "two_ends_extension " << two_ends_extension << endl;

        // metropolis enable/disable
        if (!GetValue("metropolis").empty()) {metropolis = In[0]->Get_bool(GetValue("metropolis"), false);}
        else metropolis = true;
        if (debug) cout << "metropolis " << metropolis << endl;

        // prefactor kT constant
        if (!GetValue("prefactor_kT").empty()) {
            success = In[0]->Get_Real(GetValue("prefactor_kT"), prefactor_kT, 0, 10,
                    "Prefactor_kT is a constant in Metropolis algorithm (-1/C1) * (delta Fs/kT), where C1 - prefactor.\n"
                    "Currently available range is from 0 to 10.");
            if (!success) {return success;}
        } else prefactor_kT = 1;
        if (debug) cout << "prefactor_kT is " << prefactor_kT << endl;

        // h5

        vector<string> key_name_props;

        string values2h5;
        if (!GetValue("h5").empty()) {
            success = In[0]->Get_string(GetValue("h5"), values2h5, "");
            // cutting {first, last}
            values2h5.erase(0, 1);
            values2h5.erase(values2h5.length()-1, values2h5.length());

            // parsing through ,
            string value;
            stringstream stream(values2h5);
            while( getline(stream, value, ',') ) {
                key_name_props.push_back(value);
            }
            // parsing through |
            string knp;
            for (auto && one_output : key_name_props) {
                vector<string> _;
                stringstream stream(one_output);
                while( getline(stream, knp, '|') ) {
                    _.push_back(knp);
                }
                out_per_line.push_back(_);
            }
        }
        else values2h5 = "";
        if (debug) cout << "values2h5 is " << values2h5 << endl;

        for (auto && f : out_per_line) {
            if (debug) cout << "h5_out: " << f[0] << " " << f[1] << " " << f[2] << endl;
        }

        // TODO: EXTEND CLENG
        if (success) {
            n_boxes = Seg[clamp_seg]->n_box;
            sub_box_size = {Seg[clamp_seg]->mx, Seg[clamp_seg]->my, Seg[clamp_seg]->mz};
        }
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
                return success;
            }
        }

        vector<string> sub;
        In[0]->split(In[0]->name, '.', sub);
        filename = In[0]->output_info.getOutputPath() + sub[0];

        cout << CLENG_VERSION << endl;
        t0_simulation = std::chrono::steady_clock::now();
        //success = MonteCarlo(save_vector); // testing
        success = MonteCarlo(false); //
        t1_simulation= std::chrono::steady_clock::now();
        std::cout << "# It took (s) = " << std::chrono::duration_cast<std::chrono::seconds> (t1_simulation - t0_simulation).count() <<
                  " or (m) = " << std::chrono::duration_cast<std::chrono::minutes> (t1_simulation - t0_simulation).count() <<
                  " or (h) = " << std::chrono::duration_cast<std::chrono::hours> (t1_simulation - t0_simulation).count() << std::endl;

        cout << "Have a fun. " << endl;
    }
    return success;
}

bool Cleng::CP(transfer tofrom) {
    if (debug) cout << "CP in Cleng" << endl;

    bool success = true;

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
            }
            nodes_map = createNodes(simpleNodeList, pivot_arm_nodes, pivot_arms);
            break;

        case to_segment:
            //Zero(clamped->H_MASK, Lat[0]->M);  [CLENG]
            // merged:
            std::fill(clamped->H_MASK, clamped->H_MASK+Lat[0]->M, 0);
            //std::fill(Seg[clamp_seg]->H_MASK, Seg[clamp_seg]->H_MASK+M, 0); [Master]

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
                clamped->bx[n_box] = (p2.x + p1.x - sub_box_size.x) / 2;
                clamped->by[n_box] = (p2.y + p1.y - sub_box_size.y) / 2;
                clamped->bz[n_box] = (p2.z + p1.z - sub_box_size.z) / 2;

                if (clamped->bx[n_box] < 1) { clamped->bx[n_box] += box.x; clamped->px1[n_box] += box.x; clamped->px2[n_box] += box.x;}
                if (clamped->by[n_box] < 1) { clamped->by[n_box] += box.y; clamped->py1[n_box] += box.y; clamped->py2[n_box] += box.y;}
                if (clamped->bz[n_box] < 1) { clamped->bz[n_box] += box.z; clamped->pz1[n_box] += box.z; clamped->pz2[n_box] += box.z;}
// clearing
                auto hp1x = ((clamped->px1[n_box] - 1) % box.x + 1) * J.x;
                auto hp1y = ((clamped->py1[n_box] - 1) % box.y + 1) * J.y;
                auto hp1z =  (clamped->pz1[n_box] - 1) % box.z + 1;

                auto hp2x = ((clamped->px2[n_box] - 1) % box.x + 1) * J.x;
                auto hp2y = ((clamped->py2[n_box] - 1) % box.y + 1) * J.y;
                auto hp2z = (clamped->pz2[n_box] - 1) % box.z + 1;

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

bool Cleng::MakeMove(bool back, bool inner_loop, bool cleng) {
    if (debug) cout << "MakeMove in Cleng" << endl;
    bool success = true;

    if (back) {
        cout << internal_name << "Moving back...";
        string _nodes;

        if (inner_loop) {
            //cout << "[INNER LOOP]" << endl;
            // reverse order of movements
            reverse(MC_attempt_nodeIDs_clampedMove_info.begin(), MC_attempt_nodeIDs_clampedMove_info.end());
            // go through MC attempts
            for (auto &&nodeIDs_clampedMove_per_MC : MC_attempt_nodeIDs_clampedMove_info) {
                for (auto &&nodeID_clampedMove : nodeIDs_clampedMove_per_MC) {
                    _nodes += " " + to_string(nodeID_clampedMove.first);
                    _moveClampedNode(back, nodeID_clampedMove.first, nodeID_clampedMove.second);
                }
            }
            cout << _nodes << " [Moved back]" << endl;
            MC_attempt_nodeIDs_clampedMove_info.clear();
        } else {
            //cout << "[NOT INNER LOOP]" << endl;
            for (auto &&nodeID_clampedMove : nodeIDs_clampedMove) {
                _nodes +=  " " + to_string(nodeID_clampedMove.first);
                _moveClampedNode(back, nodeID_clampedMove.first, nodeID_clampedMove.second);
            }
            cout << _nodes << " [Moved back]" << endl;
            // clean the last step
            if (!MC_attempt_nodeIDs_clampedMove_info.empty()) MC_attempt_nodeIDs_clampedMove_info.pop_back();
        }
    } else {
        nodeIDs_clampedMove.clear();
        if (inner_loop) {

            int type_move = rand.getInt(0, 2);
            if (type_move==0) success = _pivotMoveClampedNode(false);
            if (type_move==1) success = _oneNodeMoveClampedNode(false);
            if (type_move==2) success = _pivotMoveOneBond(false);

        } else {
            if (pivot_move) {
                if (pivot_one_bond) {
                    int type_move = rand.getInt(0, 1);
                    if (type_move) success = _pivotMoveClampedNode(false);
                    else success = _pivotMoveOneBond(false);
                }
                if (pivot_one_node) {
                    int type_move = rand.getInt(0, 1);
                    if (type_move) success = _pivotMoveClampedNode(false);
                    else success = _oneNodeMoveClampedNode(false);
                }
                if (!pivot_one_node and !pivot_one_bond) {
                    success = _pivotMoveClampedNode(false);
                }
            } else {
                success = _oneNodeMoveClampedNode(false);
            }
        }

        MC_attempt_nodeIDs_clampedMove_info.push_back(nodeIDs_clampedMove);
        //if (inner_loop) MC_attempt_nodeIDs_clampedMove_info.push_back(nodeIDs_clampedMove);
    }
    return success;
}

// Main procedure
bool Cleng::MonteCarlo(bool save_vector) {
    if (debug) cout << "Monte Carlo in Cleng" << endl;
    bool success = true;

    signal(SIGINT, signalHandler);
#ifdef CLENG_EXPERIMENTAL
    cleng_writer.init(filename);
#endif
    MC_attempt = 0; // initial value for loop and cleng_writer in init outlook
    mcs_done   = 0; // initial value for loop and _save_differences

    Real _accepted = 0.0;
    Real _rejected = 0.0;

    Checkpoint checkpoint;
    if (checkpoint_load) {
        if (checkpoint.isLoadable()) {
            simpleNodeList = checkpoint.loadCheckpoint(simpleNodeList, box, MC_attempt, mcs_done, _accepted, _rejected);
            nodes_map = createNodes(simpleNodeList, pivot_arm_nodes, pivot_arms);
            cout << internal_name << "From checkpoint next nodes are available: " << endl;
            for (auto &&n : nodes_map) cout << "id: " << n.first << " " << n.second.data()->get()->to_string() << endl;
            CP(to_segment);
            if (getLastMCS() != 0)  {
                MCS_checkpoint = getLastMCS() + 1;
            }
            loaded = true;
            // for stars
            //for (auto &&pair_pivot : pivot_arm_nodes) {
            //    cout << "---> arm: " << pair_pivot.first << " ids: ";
            //    for (auto &&ids: pair_pivot.second) cout << ids << " ";
            //    cout << endl;
            //}
            //if (pivot_move) {ids_node4fix.clear();ids_node4fix.push_back(pivot_arm_nodes[1][0]);}
        }
    }

    make_BC();      // check boundary conditions

// init system outlook
    cout << "Initial system outlook...\n" << endl;
    success = initSystemOutlook();
    if (!success) exit(1);

    free_energy_current = Sys[0]->GetFreeEnergy();
    if (save_vector) test_vector.push_back(free_energy_current);

    // central node of the star
    Analyzer analyzer = Analyzer(filename, box.x / 2,
                                 nodes_map[pivot_arm_nodes[1].begin()[0]].data()->get()->point(),
                                 _accepted,
                                 _rejected);

    // init save
    save(MC_attempt, analyzer);
    if (checkpoint_save) {checkpoint.saveCheckpoint(simpleNodeList, analyzer.accepted, analyzer.rejected);}
    MC_attempt += 1;
    cout << "Initialization --> done.\n" << endl;

    if (MCS) {
        auto t0_noanalysis_simulation = std::chrono::steady_clock::now();
        update_ids_node4move();
        notification();
        cout << internal_name <<  "Here we go..." << endl;
        bool success_move;

        Real F_proposed_init    = 0.0;  // init step
        Real F_proposed_current = 0.0;  // running step
        Real F_proposed_trial = 0.0;    // running step
        vector<Real> mcs_values;

        Real correction = 0.0;
        bool update_correction = true;

        for (; MC_attempt <= MCS; MC_attempt++) { // main loop for trials

            if (MC_attempt <= SILF) {
                // standard behavior -->
                success_move = MakeMove(false);
                if (success_move) {
                    CP(to_segment);
                    cout << endl;
                    notification();
                    success = solveAndCheckFreeEnergy(); // free energy is not Nan
                    if (!success) break;

                    free_energy_trial = Sys[0]->GetFreeEnergy();
                    // TESTING
                    if (save_vector) test_vector.push_back(Sys[0]->GetFreeEnergy());
                    // notification
                    cout << "Free Energy (c): " << free_energy_current << endl;
                    cout << "            (t): " << free_energy_trial << endl;
                    cout << "   prefactor kT: " << prefactor_kT << endl;
                    cout << endl;

                    Real F_proposed = getFreeEnergyBox();
                    if (escf_ebox_flag) _save_differences(0, free_energy_trial, F_proposed);

                    if (!metropolis) {
                        cout << "Metropolis is disabled. " << endl;
                        analyzer.accepted++;
                    } else {
                        if (MC_attempt <= SILF) { // check inner loop is working
                            cout << "[METROPOLIS]" << endl;
                            success = analyzer.Metropolis(rand, prefactor_kT, free_energy_trial, free_energy_current);
                        }

                        if (success) {
                            MC_attempt_nodeIDs_clampedMove_info.clear(); // clean before the next step
                            free_energy_current = free_energy_trial;
                        }
                        else { // rejection the trial step
                            MakeMove(true);
                            CP(to_segment);
                            success = solveAndCheckFreeEnergy(); //
                            if (!success) break;
                            // TODO check cleng returned to normal numbers... [2nd solution]
                            cout << "... [Done]" << endl;
                            //notification();
                            //cout << "CHECK THIS1" << endl;
                            free_energy_current = Sys[0]->GetFreeEnergy();  // in case of small tolerance could be differences
                        }
                    }

                    // notification
                    if (metropolis) { analyzer.notification(MC_attempt, cleng_rejected); }

                    // Saving data
                    auto t1_noanalysis_simulation = std::chrono::steady_clock::now();
                    tpure_simulation = std::chrono::duration_cast<std::chrono::seconds>(
                            t1_noanalysis_simulation - t0_noanalysis_simulation).count();
                    save(int(MC_attempt), analyzer);
                    if (checkpoint_save) checkpoint.updateCheckpoint(simpleNodeList, MC_attempt, mcs_done, analyzer.accepted, analyzer.rejected);
                } else MakeMove(true);
                // <-- standard behavior
            } else {
                // INNER_LOOP
                // Starting point for inner MC
                if (MC_attempt % ILE == 0) {         // if this just right MC_attempt --> entering to inner loop

                    if ( (MC_attempt % 100 == 0) && (mcs_done > 0) ) update_correction = true;
                    if ( update_correction ) {
                        Real F_proposed_ = getFreeEnergyBox();      // current proposed free energy
                        Real SCF_        = Sys[0]->GetFreeEnergy();
                        correction = abs(abs(F_proposed_) - abs(SCF_));  // absolute diff
                        if (F_proposed_ > SCF_) correction = -correction;
                        correction = 0.0;
                        cout << "Coefficient:" << correction << endl;
                        update_correction = false;
                    }

                    // inner loop
                    if (mcs_done == 0) {F_proposed_init = F_proposed_current = getFreeEnergyBox() + correction;}

                    //cout << "[INNER] F_proposed_init: " << F_proposed_init << endl;
                    F_proposed_current = F_proposed_init;

                    mcs_values.clear();
                    for (int mc_attempt = 0; mc_attempt < mcs; mc_attempt++) {
                        cout << "[" << mc_attempt << "] ";
                        success_move = MakeMove(false, true);  // storing all mc_attempt
                        if (success_move) {
                            cout << "[INNER METROPOLIS]" << endl;
                            F_proposed_trial = getFreeEnergyBox() + correction;
                            //cout << "F_proposed_trial:" << F_proposed_trial << endl;
                            success = analyzer.Metropolis(rand, prefactor_kT, F_proposed_trial, F_proposed_current, false);
                            if (success) {
                                //cout << "ACCEPTED" << endl;
                                F_proposed_current = F_proposed_trial;
                            }
                            else { // rejection the trial step
                                //cout << "REJECTED" << endl;
                                MakeMove(true);
                            }
                            mcs_values.push_back(F_proposed_current); // add value if the move is okay!

                            //cout << "Inner Free Energies: [ ";
                            //for (auto &&v : mcs_values) cout << v << " ";
                            //cout << "];" << endl;
                        }
                        else {
                            MakeMove(true, true); // restoring all mc_attempt back
                            mc_attempt = 0;       // repeat --> leads to an increase cleng_rejects
                            mcs_values.clear();   // clearing
                        }
                    }
                    // after finishing inner part update namics nodes
                    CP(to_segment);
                    cout << endl;

                    // Predicted last step box | SCF
                    F_proposed_trial = getFreeEnergyBox() + correction;
                    notification();
                    success          = solveAndCheckFreeEnergy();
                    if (!success) break;
                    free_energy_trial = Sys[0]->GetFreeEnergy();
                } else {
                    // TODO: do standard behavior;
                    //cout << "Not just right attempt." << endl;
                }
                // END INNER LOOP

                if (!metropolis) {
                    cout << "Metropolis is disabled. " << endl;
                    analyzer.accepted++;
                } else {
                    cout << "[METROPOLIS+INNER]" << endl;
                    Real free_energy        = free_energy_trial - free_energy_current;
                    Real free_energy_primed = F_proposed_trial  - F_proposed_init;

                    // notification
                    cout << "[SCF] Free Energy  (c): " << free_energy_current << endl;
                    cout << "                   (t): " << free_energy_trial << endl;
                    cout << "[Inner] Box        (c): " << F_proposed_init << endl;
                    cout << "                   (t): " << F_proposed_trial << endl;
                    cout << "[Outer/Inner] (energy): "  << free_energy         << endl;
                    cout << "             (energy'): " << free_energy_primed << endl;
                    //cout << "     (energy - energy'): " << free_energy - free_energy_primed << endl;
                    cout << "prefactor kT: " << prefactor_kT << endl;
                    cout << "  correction: " << correction << endl;
                    cout << endl;
                    success = analyzer.OuterInner(rand, prefactor_kT, free_energy, free_energy_primed);

                    //cout << "WHAT WE HAVE BEFORE -->" << endl;
                    //notification();
                    //cout << "<-- WHAT WE HAVE BEFORE" << endl;

                    if (success) {
                        // set of inner step is accepted
                        MC_attempt_nodeIDs_clampedMove_info.clear(); // clean before the next step
                        free_energy_current = free_energy_trial;
                        F_proposed_init     = F_proposed_trial;

                        // saving initial and last before the last one
                        int mcs_done_ = mcs_done;
                        if (escf_ebox_flag) {
                            for (size_t idx = 0; idx < mcs_values.size(); idx++) {
                                _save_differences(static_cast<int>(mcs_done_+idx), nan("1"), mcs_values[idx]);
                            }
                        }
                        mcs_done += mcs; // increase mcs done steps
                        if (escf_ebox_flag) _save_differences(mcs_done, free_energy_trial, F_proposed_trial);
                    }
                    else { // rejection the trial step
                        MakeMove(true, true);
                        CP(to_segment);
                        
                        //cout << "WHAT WE HAVE NOW -->" << endl;
                        //notification();
                        //cout << "<-- WHAT WE HAVE NOW" << endl;
                        
                        success = solveAndCheckFreeEnergy(); //
                        if (!success) break;
                        // TODO check cleng returned to normal numbers... [2nd solution]
                        cout << "... [Done]" << endl;
                        notification();
                        //cout << "CHECK THIS" << endl;
                        free_energy_current = Sys[0]->GetFreeEnergy();  // in case of small tolerance could be differences
                    }
                }
                // notification
                if (metropolis) {analyzer.notification(MC_attempt, cleng_rejected, mcs_done);}

                // Saving data
                auto t1_noanalysis_simulation = std::chrono::steady_clock::now();
                tpure_simulation = std::chrono::duration_cast<std::chrono::seconds> (t1_noanalysis_simulation - t0_noanalysis_simulation).count();
//                save(int(MC_attempt+MCS_checkpoint-cleng_rejected), analyzer);
                save(int(MC_attempt), analyzer);
                if (checkpoint_save) checkpoint.updateCheckpoint(simpleNodeList, MC_attempt, mcs_done, analyzer.accepted, analyzer.rejected);

            } // <-- INNER

            if (cleng_flag_termination) break;
            cout << endl;
        }  // main loop

        cout << endl;
        cout << "[Finally]" << endl;
        analyzer.notification(MC_attempt-1, cleng_rejected, mcs_done);
    }
    return success;
}

