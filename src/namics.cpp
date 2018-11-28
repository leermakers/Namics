#define MAINxH
#include "alias.h"
#include "input.h"
#include "lattice.h"
#include "molecule.h"
#include "namics.h"
#include "mesodyn.h"
#include "cleng/cleng.h"
#include "teng.h"
//#include "newton.h"
#include "output.h"
#include "segment.h"
#include "state.h"
#include "reaction.h"
#include "system.h"
#include "tools.h"
#include "variate.h"
#include "sfnewton.h"
#include "solve_scf.h"

string version = "2.1.1.1.1.1.1";
// meaning:
// newton version number =2
// system version number =1
// lattice version number =1
// molecule version number =1
// segment version number =1
// alias version number =1
// output version number =1
Real check = 0.4534345;
Real e = 1.60217e-19;
Real T = 298.15;
Real k_B = 1.38065e-23;
Real k_BT = k_B * T;
Real eps0 = 8.85418e-12;
Real PIE = 3.14159265;
int DEBUG_BREAK = 1;
//Used for command line switches
bool debug = false;
bool suppress = false;

//Output when the user malforms input. Update when adding new command line switches.
void improperInput() {
  cerr << "Improper usage: namics [filename] [-options]." << endl << "Options available:" << endl;
  cerr << "-d Enables debugging mode." << endl;
  cerr << "-s Suppresses newton's extra output." << endl;
}

int main(int argc, char* argv[]) {
	vector<string> args(argv, argv+argc);
//Output error if no filename has been specified.
 	if (argc == 1) {
    		improperInput();
    		return 1;
  	}
//Output error if user starts with a commandline switch. (Also catches combination with forgotten filename)
  	if ( (args[1])[0] == '-' ) {
    		improperInput();
    		return 1;
  	}

// If the specified filename has no extension: add the extension specified below.
	string extension = "in";
	ostringstream filename;
	filename << argv[1];
	bool hasNoExtension = (filename.str().substr(filename.str().find_last_of(".")+1) != extension);
 	if (hasNoExtension) filename << "." << extension;

//If the switch -d is given, enable debug. Add new switches by copying and replacing -d and debug = true.
  	if ( find(args.begin(), args.end(), "-d") != args.end() ) {
    		debug = true;
 	 }

  	if ( find(args.begin(), args.end(), "-s") != args.end() ) {
    		suppress = true;
  	}

  	bool cuda;
  	int start = 0;
  	int n_starts = 0;

  	string initial_guess;
  	string final_guess;
  	string METHOD = "";
  	Real* X = NULL;
  	int MX = 0, MY = 0, MZ = 0;
  	int fjc_old = 0;
  	bool CHARGED = false;
  	vector<string> MONLIST;
	vector<string> STATELIST;

#ifdef CUDA
  GPU_present();
  cuda = true;
#else
  cuda = false;
#endif

// All class instances are stored in the following vectors
  	vector<Input*> In;     // Inputs read from file
  	vector<Output*> Out;   // Outputs written to file
  	vector<Lattice*> Lat;  // Properties of the lattice
  	vector<Molecule*> Mol; // Properties of entire molecule
  	vector<Segment*> Seg;  // Properties of molecule segments
	vector<State*> Sta;
	vector<Reaction*> Rea;
  	vector<Solve_scf*> New;   // Solvers and iteration schemes
  	vector<System*> Sys;
  	vector<Variate*> Var;
  	vector<Mesodyn*> Mes;
  	vector<Cleng*> Cle; //enginge for clampled molecules
  	vector<Teng*> Ten; //enginge for clampled molecules

// Create input class instance and handle errors(reference above)
  	In.push_back(new Input(filename.str()) );
  	if (In[0]->Input_error) {
    		return 0;
  	}
  	n_starts = In[0]->GetNumStarts();
  	if (n_starts == 0)
    		n_starts++; // Default to 1 start..

/******** This while loop basically contains the rest of main, initializes all classes and performs  ********/
/******** calculations for a given number (start) of cycles ********/

	while (start < n_starts) {

    		start++;
   		In[0]->MakeLists(start);
    		cout << "Problem nr " << start << " out of " << n_starts << endl;

/******** Class creation starts here ********/

// Create lattice class instance and check inputs (reference above)
    		Lat.push_back(new Lattice(In, In[0]->LatList[0]));
    		if (!Lat[0]->CheckInput(start)) {
      			return 0;
    		}

// Create segment class instance and check inputs (reference above)
    		int n_seg = In[0]->MonList.size();
    		for (int i = 0; i < n_seg; i++)
      			Seg.push_back(new Segment(In, Lat, In[0]->MonList[i], i, n_seg));

//Create state class instance and check inputs
		int n_stat=In[0]->StateList.size();
		for (int i=0; i<n_stat; i++)
			Sta.push_back(new State(In,Seg,In[0]->StateList[i]));


    		for (int i = 0; i < n_seg; i++) {
      			for (int k = 0; k < n_seg; k++) {
        			Seg[i]->PutChiKEY(Seg[k]->name);
     			}
			for (int k = 0; k < n_stat; k++) {
        			Seg[i]->PutChiKEY(Sta[k]->name);
     			}
      			if (!Seg[i]->CheckInput(start)) return 0;
    		}
		for (int i=0; i<n_stat; i++) {
      			for (int k = 0; k < n_seg; k++) {
        			Sta[i]->PutChiKEY(Seg[k]->name);
     			}
			for (int k = 0; k < n_stat; k++) {
        			Sta[i]->PutChiKEY(Sta[k]->name);
     			}
			if (!Sta[i]->CheckInput(start)) return 0;
		}

//Create reaction class instance and check inputs
		int n_rea=In[0]->ReactionList.size();
		for (int i=0; i<n_rea; i++) {
			Rea.push_back(new Reaction(In,Seg,Sta,In[0]->ReactionList[i]));
			if (!Rea[i]->CheckInput(start)) return 0;
		}


// Create segment class instance and check inputs (reference above)
    		int n_mol = In[0]->MolList.size();
    		for (int i = 0; i < n_mol; i++) {
      			Mol.push_back(new Molecule(In, Lat, Seg, In[0]->MolList[i]));
      			if (!Mol[i]->CheckInput(start))
        			return 0;
   		}

// Create system class instance and check inputs (reference above)
    		Sys.push_back(new System(In, Lat, Seg, Sta, Rea, Mol, In[0]->SysList[0]));
    		Sys[0]->cuda = cuda;
    		if (!Sys[0]->CheckInput(start)) {
      			return 0;
    		}
    		if (!Sys[0]->CheckChi_values(n_seg))
      			return 0;

// Prepare variables used in variate class creation
// TODO: What do these variables do?
    		int n_var = In[0]->VarList.size();
    		int n_search = 0;
    		int n_scan = 0;
    		int n_ets = 0;
    		int n_etm = 0;
    		int n_target = 0;
    		int search_nr = -1, scan_nr = -1, target_nr = -1, ets_nr = -1, etm_nr = -1;

// Create variate class instance and check inputs (reference above)
    		for (int k = 0; k < n_var; k++) {
      			Var.push_back(new Variate(In, Lat, Seg, Sta, Rea, Mol, Sys, In[0]->VarList[k]));
      			if (!Var[k]->CheckInput(start)) {
        		return 0;
      			}

      			if (Var[k]->scanning > -1) {
        			scan_nr = k;
        			n_scan++;
     			}
      			if (Var[k]->searching > -1) {
        			search_nr = k;
        			n_search++;
      			}
      			if (Var[k]->targeting > -1) {
        			target_nr = k;
        			n_target++;
      			}
      			if (Var[k]->eq_to_solvating > -1) {
        			ets_nr = k;
        			n_ets++;
      			}
      			if (Var[k]->eq_to_mu > -1) {
        			etm_nr = k;
        			n_etm++;
      			}
		}


// Error code for faulty variate class creation
		if (n_etm > 1) {
      			cout << "too many equate_to_mu's in var statements. The limit is 1 " << endl;
      			return 0;
    		}
    		if (n_ets > 1) {
      			cout << "too many equate_to_solvent's in var statements. The limit is 1 " << endl;
      			return 0;
    		}
    		if (n_search > 1) {
      			cout << "too many 'search'es in var statements. The limit is 1 " << endl;
      			return 0;
    		}
    		if (n_scan > 1) {
      			cout << "too many 'scan's in var statements. The limit is 1 " << endl;
      			return 0;
    		}
    		if (n_target > 1) {
      			cout << "too many 'target's in var statements. The limit is 1 " << endl;
      			return 0;
    		}
    		if (n_search > 0 && n_target == 0) {
      			cout << "lonely search. Please specify in 'var' a target function, e.g. 'var : sys-NN : grand_potential : 0'" << endl;
      			return 0;
    		}
    		if (n_target > 0 && n_search == 0) {
      			cout << "lonely target. Please specify in 'var' a search function, e.g. 'var : mol-lipid : search : theta'" << endl;
      			return 0;
    		}

// Create newton class instance and check inputs (reference above)
    		New.push_back(new Solve_scf(In, Lat, Seg, Sta, Rea, Mol, Sys, Var, In[0]->NewtonList[0]));
    		if (!New[0]->CheckInput(start)) {
      			return 0;
    		}

//Guesses geometry
    		if (Sys[0]->initial_guess == "file") {
      			MONLIST.clear();
			      STATELIST.clear();
      			if (!Lat[0]->ReadGuess(Sys[0]->guess_inputfile, X, METHOD, MONLIST,STATELIST, CHARGED, MX, MY, MZ, fjc_old, 0)) {
// last argument 0 is to first checkout sizes of system.
        			return 0;
      			}
			int nummon = MONLIST.size();
			int numstate=STATELIST.size();
      			int m;
      			if (MY == 0) {
      				m = MX + 2;
      			} else {
        			if (MZ == 0) {
         				m = (MX + 2) * (MY + 2);
        			} else {
          				m = (MX + 2) * (MY + 2) * (MZ + 2);
        			}
      			}
      			int IV = (nummon+numstate) * m;

      			if (CHARGED)
        			IV += m;
      			if (start > 1){
        			free(X);
              			X = (Real*)malloc(IV * sizeof(Real));
            		}
      			MONLIST.clear();
			STATELIST.clear();
      			Lat[0]->ReadGuess(Sys[0]->guess_inputfile, X, METHOD, MONLIST, STATELIST,CHARGED, MX, MY, MZ, fjc_old, 1);
// last argument 1 is to read guess in X.
		}

    		int substart = 0;
    		int subloop = 0;
    		if (scan_nr < 0)
      			substart = 0;
   		else
      			substart = Var[scan_nr]->num_of_cals;
    		if (substart < 1)
      			substart = 1; // Default to 1 substart

		EngineType TheEngine;
		TheEngine=SCF;
		if (In[0]->MesodynList.size()>0) {TheEngine=MESODYN;}
		if (In[0]->ClengList.size()>0) {TheEngine=CLENG;}
		if (In[0]->TengList.size()>0) {TheEngine=TENG;}

		int ii;
		int n_out=0;
		switch(TheEngine) {
			case SCF:
				// Prepare, catch errors for output class creation
    				n_out = In[0]->OutputList.size();
    				if (n_out == 0)
      					cout << "Warning: no output defined!" << endl;

				// Create output class instance and check inputs (reference above)
    				for (int ii = 0; ii < n_out; ii++) {
      					Out.push_back(new Output(In, Lat, Seg,Sta,Rea, Mol, Sys, New, In[0]->OutputList[ii], ii, n_out));
      					if (!Out[ii]->CheckInput(start)) {
        					cout << "input_error in output " << endl;
        					return 0;
      					}
    				}

				while (subloop < substart) {
      					if (scan_nr > -1)
        					Var[scan_nr]->PutVarScan(subloop);
      					New[0]->AllocateMemory();
      					New[0]->Guess(X, METHOD, MONLIST,STATELIST,CHARGED, MX, MY, MZ,fjc_old);
					if (search_nr < 0 && ets_nr < 0 && etm_nr < 0) {
						bool print=true;
       	  				New[0]->Solve(print);
          				} else {
            					if (!debug) cout << "Solve towards superiteration " << endl;
            					New[0]->SuperIterate(search_nr, target_nr, ets_nr, etm_nr);
         				}
      					New[0]->PushOutput();

      					for (ii = 0; ii < n_out; ii++) {
        					Out[ii]->WriteOutput(subloop);
      					}
      					subloop++;
    				}

				break;
			case MESODYN:
				New[0]->AllocateMemory();
      				New[0]->Guess(X, METHOD, MONLIST,STATELIST,CHARGED, MX, MY, MZ,fjc_old);
				if (debug) cout << "Creating mesodyn" << endl;
        			Mes.push_back(new Mesodyn(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->MesodynList[0]));
        			if (!Mes[0]->CheckInput(start)) {
          				return 0;
        			}
              try {
        			  Mes[0]->mesodyn();
              } catch (error RC) {
                cout << "Exiting Mesodyn with error: " << RC << ". Please check the documentation for details." << endl;
                exit(RC);
              }
              delete Mes[0];
              Mes.clear();
				break;
			case CLENG:
				New[0]->AllocateMemory();
      				New[0]->Guess(X, METHOD, MONLIST,STATELIST, CHARGED, MX, MY, MZ,fjc_old);
				if (!debug) cout << "Creating Cleng module" << endl;
				Cle.push_back(new Cleng(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->ClengList[0]));
				if (!Cle[0]->CheckInput(start)) {return 0;}
				break;
			case TENG:
				New[0]->AllocateMemory();
      				New[0]->Guess(X, METHOD, MONLIST,STATELIST,CHARGED, MX, MY, MZ,fjc_old);
				cout << "Solving Teng problem" << endl;
				Ten.push_back(new Teng(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->TengList[0]));
				if (!Ten[0]->CheckInput(start)) {return 0;}
				break;
			default:
				cout <<"TheEngine is unknown. Programming error " << endl; return 0;
				break;
		}

        for (auto all_segments : Seg)
          all_segments->prepared = false;

    		if (scan_nr > -1)
      			Var[scan_nr]->ResetScanValue();
    		if (Sys[0]->initial_guess == "previous_result") {
      			METHOD = New[0]->SCF_method; //check this..
      			MX = Lat[0]->MX;
      			MY = Lat[0]->MY;
      			MZ = Lat[0]->MZ;
      			CHARGED = Sys[0]->charged;
            		int IV_new=New[0]->iv; //check this
      			if (start > 1 || (start == 1 && Sys[0]->initial_guess == "file"))
                		free(X);
            		X = (Real *)malloc(IV_new * sizeof(Real));
                #ifdef CUDA
                TransferDataToHost(X, New[0]->xx, IV_new);
                #else
            		for (int i = 0; i < IV_new; i++) X[i] = New[0]->xx[i];
                #endif
      			fjc_old=Lat[0]->fjc;
      			int mon_length = Sys[0]->ItMonList.size();
			int state_length=Sys[0]->ItStateList.size();
      			MONLIST.clear();
			STATELIST.clear();
      			for (int i = 0; i < mon_length; i++) {
       				MONLIST.push_back(Seg[Sys[0]->ItMonList[i]]->name);
      			}
      			for (int i = 0; i < state_length; i++) {
       				STATELIST.push_back(Sta[Sys[0]->ItStateList[i]]->name);
      			}
    		}
    		if (Sys[0]->final_guess == "file") {
      			MONLIST.clear();
			STATELIST.clear();
      			int mon_length = Sys[0]->ItMonList.size();
			int state_length=Sys[0]->ItStateList.size();
      			for (int i = 0; i < mon_length; i++) {
        			MONLIST.push_back(Seg[Sys[0]->ItMonList[i]]->name);
      			}
      			for (int i = 0; i < state_length; i++) {
        			STATELIST.push_back(Sta[Sys[0]->ItStateList[i]]->name);
      			}
      			Lat[0]->StoreGuess(Sys[0]->guess_outputfile, New[0]->xx, New[0]->SCF_method, MONLIST, STATELIST, Sys[0]->charged, start);
   		}

/******** Clear all class instances ********/

    		for (int i = 0; i < n_out; i++) delete Out[i];
		    Out.clear();
   	 	  for (int i = 0; i < n_var; i++) delete Var[i];
    		Var.clear();
   	 	  delete New[0];
    		New.clear();
   	 	  delete Sys[0];
    		Sys.clear();
   	 	  for (int i = 0; i < n_mol; i++) delete Mol[i];
    		Mol.clear();
    		for (int i = 0; i < n_seg; i++) delete Seg[i];
   	 	  Seg.clear();
		for (int i=0; i<n_stat; i++) delete Sta[i];
		  Sta.clear();
		for (int i=0; i<n_rea; i++) delete Rea[i];
	          Rea.clear();
    		delete Lat[0];
    		Lat.clear();
	} //loop over starts.
  free(X);
  delete In[0];
  In.clear();
return 0;
}
