#define MAINxH
#include "alias.h"
#include "engine.h"
#include "input.h"
#include "lattice.h"
#include "molecule.h"
#include "namics.h"
#include "newton.h"
#include "output.h"
#include "segment.h"
#include "system.h"
#include "tools.h"
#include "variate.h"

string version = "1.1.1.1.1.1.1.1";
// meaning:
// engine version number =1
// newton version number =1
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

bool debug = false;

int main(int argc, char *argv[]) {

  string fname;
  string filename;
  if (argc == 2)
    fname = argv[1];
  else {
    printf("Use: namics filename -without extension- \n");
    return 1;
  }
  // argc counts no. of input arguments ./namics fname(2 arguments).... and
  // makes fname from argv[1] else shows error and stops.
  filename = fname + ".in";

  // TODO: I think these two are no longer used, but I'm too scared to delete
  // them now - Daniel
  string line;
  ofstream out_file;

  bool cuda;
  int start = 0;
  int n_starts = 0;

  string initial_guess;
  string final_guess;
  string METHOD = "";
  Real *X = NULL;
  int MX = 0, MY = 0, MZ = 0;
  int fjc_old; 
  bool CHARGED = false;
  vector<string> MONLIST;

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
  vector<Newton*> New;   // Solvers and iteration schemes
  vector<System*> Sys;
  vector<Engine*> Eng;
  vector<Variate*> Var;

  // Create input class instance and handle errors(reference above)
  In.push_back(new Input(filename));
  if (In[0]->Input_error) return 0;
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
    for (int i = 0; i < n_seg; i++) {
      for (int k = 0; k < n_seg; k++) {
        Seg[i]->PutChiKEY(Seg[k]->name);
      }
      if (!Seg[i]->CheckInput(start))
        return 0;
    }

    // Create segment class instance and check inputs (reference above)
    int n_mol = In[0]->MolList.size();
    for (int i = 0; i < n_mol; i++) {
      Mol.push_back(new Molecule(In, Lat, Seg, In[0]->MolList[i]));
      if (!Mol[i]->CheckInput(start))
        return 0;
    }

    // Create system class instance and check inputs (reference above)
    Sys.push_back(new System(In, Lat, Seg, Mol, In[0]->SysList[0]));
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
      Var.push_back(new Variate(In, Lat, Seg, Mol, Sys, In[0]->VarList[k]));
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
      cout
          << "too many equate_to_solvent's in var statements. The limit is 1 " << endl;
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
    New.push_back( new Newton(In, Lat, Seg, Mol, Sys, Var, In[0]->NewtonList[0]) );
    if (!New[0]->CheckInput(start)) {
      return 0;
    }

    // Create engine class instance and check inputs (reference above)
    Eng.push_back( new Engine(In, Lat, Seg, Mol, Sys, New, In[0]->EngineList[0]) );
    if (!Eng[0]->CheckInput(start)) {
      return 0;
    }

    // Prepare, catch errors for output class creation
    int n_out = In[0]->OutputList.size();
    if (n_out == 0)
      cout << "Warning: no output defined" << endl;

    // Create output class instance and check inputs (reference above)
    for (int i = 0; i < n_out; i++) {
      Out.push_back( new Output(In, Lat, Seg, Mol, Sys, New, Eng, In[0]->OutputList[i], i, n_out) );
      if (!Out[i]->CheckInput(start)) {
        cout << "input_error in output " << endl;
        return 0;
      }
    }

    // TODO: 1. sys.initial_guess is set to "previous_result" by default at
    // system.cpp:301
    //	2. What does this block do, exactly? It seems to handle some code for
    // charged molecules.
    //	3. Can we move this to the system class bit?

    if (Sys[0]->initial_guess == "file") {
      MONLIST.clear();
      if (!Lat[0]->ReadGuess(Sys[0]->guess_inputfile, X, METHOD, MONLIST, CHARGED, MX, MY, MZ, fjc_old, 0)) {
        // last argument 0 is to first checkout sizes of system.
        return 0;
      }
      int nummon = MONLIST.size();
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
      int IV = nummon * m;

      if (CHARGED)
        IV += m;
      if (start > 1)
        free(X);
      X = (Real *)malloc(IV * sizeof(Real));
      MONLIST.clear();
      Lat[0]->ReadGuess(Sys[0]->guess_inputfile, X, METHOD, MONLIST, CHARGED, MX, MY, MZ, fjc_old, 1);
      // last argument 1 is to read guess in X.
    }

/********** All classes have been created and input gathered, time to start some calculations *********/
		int substart=0; int subloop=0;

		if (scan_nr<0) substart=0;
		else substart = Var[scan_nr]->num_of_cals;
		if (substart<1) substart=1; // Default to 1 substart

		while(subloop < substart){
			if (scan_nr >-1) Var[scan_nr]->PutVarScan(subloop);
			New[0]->AllocateMemory();
			New[0]->Guess(X,METHOD,MONLIST,CHARGED,MX,MY,MZ,fjc_old);

			// This is the starting point of all calculations. Solve provides the flow through computation schemes.
			if (search_nr<0 && ets_nr < 0 && etm_nr <0) {
				New[0]->Solve(true);
			} else {
				New[0]->SuperIterate(search_nr,target_nr,ets_nr,etm_nr);
			}

/********** Output information about all classes to file *********/
      Lat[0]->PushOutput();
      New[0]->PushOutput();
      Eng[0]->PushOutput();
      int length = In[0]->MonList.size();
      for (int i = 0; i < length; i++)
        Seg[i]->PushOutput();
      length = In[0]->MolList.size();
      for (int i = 0; i < length; i++) {
        int length_al = Mol[i]->MolAlList.size();
        for (int k = 0; k < length_al; k++) {
          Mol[i]->Al[k]->PushOutput();
        }
        Mol[i]->PushOutput();
      }
      // length = In[0]->AliasList.size();
      // for (int i=0; i<length; i++) Al[i]->PushOutput();
      Sys[0]->PushOutput(); // needs to be after pushing output for seg.

      for (int i = 0; i < n_out; i++) {
        Out[i]->WriteOutput(subloop);
      }
      subloop++;
    }

/******** Clear all class instances ********/
    // TODO: Again some code for ..charged molecules?
    if (scan_nr > -1)
      Var[scan_nr]->ResetScanValue();
    if (Sys[0]->initial_guess == "previous_result") {
      int IV_new;
      METHOD = New[0]->GetNewtonInfo(IV_new);
      MX = Lat[0]->MX;
      MY = Lat[0]->MY;
      MZ = Lat[0]->MZ;
      CHARGED = Sys[0]->charged;
      if (start > 1 || (start == 1 && Sys[0]->initial_guess == "file"))
        free(X);
      X = (Real *)malloc(IV_new * sizeof(Real));
      for (int i = 0; i < IV_new; i++) X[i] = New[0]->xx[i];
      fjc_old=Lat[0]->fjc; 
      int length = Sys[0]->SysMonList.size();
      MONLIST.clear();
      for (int i = 0; i < length; i++) {
        MONLIST.push_back(Seg[Sys[0]->SysMonList[i]]->name);
      }
    }
    if (Sys[0]->final_guess == "file") {
      MONLIST.clear();
      int length = Sys[0]->SysMonList.size();
      for (int i = 0; i < length; i++) {
        MONLIST.push_back(Seg[Sys[0]->SysMonList[i]]->name);
      }
      Lat[0]->StoreGuess(Sys[0]->guess_outputfile, New[0]->xx, New[0]->method, MONLIST, Sys[0]->charged, start);
    }
    for (int i = 0; i < n_out; i++)
      delete Out[i];
    Out.clear();
    delete Eng[0];
    Eng.clear();
    for (int i = 0; i < n_var; i++)
      delete Var[i];
    Var.clear();
    delete New[0];
    New.clear();
    delete Sys[0];
    Sys.clear();
    for (int i = 0; i < n_mol; i++) {
      // Mol[i]->DeleteAl();
      // int length_al = Mol[i]->MolAlList.size();
      // for (int k=0; k<length_al; k++) {
      // delete Mol[i]->Al[k]; Mol[i]->Al.clear();
      //}
      delete Mol[i];
    }
    Mol.clear();
    for (int i = 0; i < n_seg; i++)
      delete Seg[i];
    Seg.clear();
    delete Lat[0];
    Lat.clear();
  }
  return 0;
}
