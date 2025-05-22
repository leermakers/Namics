#include "microemulsion.h"
#include "tools.h"

Microemulsion::Microemulsion(vector<Input *> In_, vector<Output *> Out_, vector<Lattice *> Lat_, vector<Segment *> Seg_, vector<State *> Sta_, vector<Reaction *> Rea_, vector<Molecule *> Mol_, vector<System *> Sys_,vector<Solve_scf*> New_,vector<Variate*> Var_,string name_):
	name{name_}, In{In_}, Out{Out_}, Lat{Lat_}, Mol{Mol_}, Seg{Seg_}, Sta{Sta_}, Rea{Rea_}, Sys{Sys_}, New{New_}, Var{Var_}
{
	if (debug)cout << "Constructor for microemulsion " << endl;
  	KEYS.push_back("surfactant");
	KEYS.push_back("co_solvent");
}
Microemulsion::~Microemulsion(){}

bool Microemulsion::CheckInput(int start_)
{
	if (debug)
		cout << "CheckInput for system " << endl;
	start=start_;
	bool success = true;
	int surfactant=-1;
	int co_solvent=-1;
	string molname;
	int length = In[0]->MolList.size();
	success = In[0]->CheckParameters("micro", name, start, KEYS, PARAMETERS, VALUES);
	if (success) {
		if (GetValue("surfactant").size()>0) {
			molname=GetValue("surfactant");
			for (int i=0; i<length; i++) {
				if (In[0]->MolList[i]==molname) surfactant=i;
			}
			if (surfactant<0) {
				success=false; cout <<"In microemulsions you need to specify the surfactant: use : micro : 'name' : surfactant : 'molname' " << endl;
			}
		}
		if (GetValue("co_solvent").size()>0) {
			molname=GetValue("co_solvent");
			for (int i=0; i<length; i++) {
				if (In[0]->MolList[i]==molname) co_solvent=i;
			}
			if (co_solvent<0) {
				success=false; cout <<"In microemulsions you need to specify the co_solvent: use : micro : 'name' : co_solvent : 'molname' " << endl;
			}
		}
	}
	return success;
}

string Microemulsion::GetValue(string parameter)
{
	if (debug)
		cout << "GetValue " + parameter + " for microemulsion " << endl;
	int length = PARAMETERS.size();
	for (int i = 0; i < length; ++i)
	{
		if (parameter == PARAMETERS[i])
		{
			return VALUES[i];
		}
	}
	return "";
}

bool Microemulsion:: FixedPoint() {
	if (debug) cout << "Fixed point in microemulsions " << endl;
	Sys[0]->MakeItsLists();
	New[0]->AllocateMemory();
	New[0]->Guess(X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);
	if (search_nr < 0 && ets_nr < 0 && etm_nr < 0) {
		New[0]->Solve(true);
	} else {
		if (debug) cout << "Solve towards superiteration " << endl;
		New[0]->SuperIterate(search_nr, target_nr, ets_nr, etm_nr, bm_nr);
	}
	return true;
}

bool Microemulsion:: WriteResults() {
	bool kal_append=false;
	int n_out = In[0]->OutputList.size();

	for (int ii = 0; ii < n_out; ii++) {
		Out.push_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[ii], ii, n_out));
		if (!Out[ii]->CheckInput(start)){
			cout << "input_error in output " << endl;
			return 0;
		} else {
			if (Out[ii]->name=="kal") { //this is to make sure that append is set to true when 'kal'-file is not initiated for the first time.
				if (kal_append) Out[ii]->append=true;
				else kal_append=true;
			}
		}
	}
	New[0]->PushOutput();

	for (int ii = 0; ii < n_out; ii++){
		Out[ii]->WriteOutput(subloop);
	}
	if (Sys[0]->final_guess == "file"){ //if iv have changed: see namics variant.
		Lat[0]->StoreGuess(Sys[0]->guess_outputfile, New[0]->xx , METHOD, MONLIST, STATELIST, CHARGED, start);
	}
	subloop++;
	return true;
}

bool Microemulsion::Doit(Real* X_,string METHOD_,vector<string> MONLIST_,vector<string> STATELIST_,bool CHARGED_,int MX_,int MY_,int MZ_,int fjc_old_,int search_nr_,int ets_nr_,int etm_nr_,int target_nr_,int bm_nr_,int subloop_) {
	X=X_; METHOD=METHOD_; MONLIST=MONLIST_; STATELIST=STATELIST_; CHARGED=CHARGED_; MX=MX_; MY=MY_; MZ=MZ_; fjc_old=fjc_old_; search_nr=search_nr_; ets_nr=ets_nr_;  etm_nr=etm_nr_; target_nr=target_nr_; bm_nr=bm_nr_; subloop=subloop_;
	if (debug) cout <<"I'll do it " << endl;

	FixedPoint();
	WriteResults();

	return true;
}

//now borrow algorithm from jupyter nootbook.... and go


