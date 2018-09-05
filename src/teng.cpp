#include "teng.h"
#include "output.h"

Teng::Teng(vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<State*> Sta_,vector<Reaction*> Rea_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Solve_scf*> New_,  string name_)
    : name{name_},
      In{In_},
      Lat{Lat_},
      Mol{Mol_},
      Seg{Seg_},
      Sta{Sta_},
      Rea{Rea_},
      Sys{Sys_},
      New{New_}

{
	if (debug) cout << "Teng initialized" << endl;
 	KEYS.push_back("MCS");
  	KEYS.push_back("save_interval");
  	KEYS.push_back("save_filename");
  	KEYS.push_back("seed");
}

Teng::~Teng() {
}

// MonteCarlo Engine to drive the tagged molecules.
bool Teng::MonteCarlo() {
  if (debug) cout << "Monte Carlo in Teng" << endl;
	bool success;
 	success=CP(to_teng);
	New[0]->Solve(true);
	Real F_bm = Sys[0]->FreeEnergy;
	Real F_am;
	WriteOutput(0);
		for (int time=1; time<=MCS; time++){
			//Copy x,y,z to x_bm, y_bm, z_bm
			cout << "Storing a copy of molecular positions"<< endl;
			success=CP(to_bm);
			//Make a montecarlo move
			cout << "Changing the mode of the particles"<< endl;
			ChangeMode();
			//copy moved positions to segment
			cout << "Copying new mode molecular positions into segment class"<< endl;
			success=CP(to_segment);
			//Solve SCF
			New[0]->Solve(true);
			F_am = Sys[0]->FreeEnergy;
			//Decide to accept or reject
			//Accept
			cout << F_bm << "	:" << F_am << "	:" << GetRandom(1.0) << endl; 
			if (F_am-F_bm <= 0 || GetRandom(1.0) < exp(F_am-F_bm)){
				//Change F_bm to sys->Freeenergy
				cout << "Monte Carlo move is accepted" << endl;
				F_bm=Sys[0]->FreeEnergy;
			//Reject
			} else{
				//Copy pos_bm to segment (basically reset configuration)
				cout << "Monte Carlo move is rejected" << endl;
				success = CP(reset);
				//Sys->Freeenergy should be set to F_bm
				Sys[0]->FreeEnergy = F_bm;
			}
			if (time%save_interval ==0)WriteOutput(time);
		}
	return success;
}

int Teng::GetRandom(int maxvalue) {
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> distance(1,maxvalue);
        int randomnumber = distance(gen);
	return randomnumber;
}

Real Teng::GetRandom(Real maxvalue){
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> distance(0,maxvalue);
        Real randomnumber = distance(gen);
	return randomnumber;
}

bool Teng::ChangeMode(){
	bool success=false;
	while(!success){
	Real Amplitude;
	Real Wavenumber;
	Real pi=4.0*atan(1.0);
	for(int i=0;i<n_particles;i++){
		Amplitude=GetRandom(3.0);
		Wavenumber=GetRandom(2.0);
		X[i]=X[i]+round(0.5-round(GetRandom(1.0)));
		Y[i]=Y[i]+round(0.5-round(GetRandom(1.0)));
		Z[i]=Lat[0]->MZ/2 + round(Amplitude*(sin(Wavenumber*pi*X[i]/Lat[0]->MX)*sin(Wavenumber*pi*Y[i]/Lat[0]->MY))); 		//Move should be dependent on the system size rather than arbitrary 10.
		cout << "Position of particle " << i << ": (" <<X[i] << ","<< Y[i] <<"," << Z[i] <<")" << endl; //Not necessary. Should be in debug.
		}
	success=CP(to_segment);
	success=IsLegal();
	}
	return success;
}

// Checks for illegal montecarlo moves. i.e., particle collisions alone. Boundary checky is still not implemented. 
// However, boundary checks might be redundant in mode based montecarlo as the modes can never be moving particles out of boundary.
bool Teng::IsLegal(){
	bool success=true;
	int i,j;
	for(i=0; i<n_particles; i++){
		for (j=0; j<i; j++){
		if(i!=j){
			if(Seg[tag_seg]->H_P[i]==Seg[tag_seg]->H_P[j]) {success=false; cout << "PARTICLES COLLIDED" << endl;}
			}
		}
	}
	return success;
}

// Transfer the particle locations from segment to teng, and vice versa.,
bool Teng::CP(transfer tofrom) {
	int JX=Lat[0]->JX;
	int JY=Lat[0]->JY;
	int M = Lat[0]->M;

	bool success=true;
	int i;
	switch(tofrom) {
		case to_teng:
			for (i=0; i<n_particles; i++) {
				PX=Seg[tag_seg]->H_P[i]/JX;
				PY=(Seg[tag_seg]->H_P[i]-PX*JX)/JY;
				PZ=(Seg[tag_seg]->H_P[i]-PX*JX-PY*JY);
				X.push_back(PX); Y.push_back(PY); Z.push_back(PZ);
			}
		break;
		case to_segment:
			Zero(Seg[tag_seg]->H_MASK,M);  //TODO: (Ram) Not sure If I have to remove the mask, all calculations use H_P after this.
			for (i=0; i<n_particles; i++) {
				Seg[tag_seg]->H_P[i]=X[i]*JX+Y[i]*JY+Z[i];
				Seg[tag_seg]->H_MASK[X[i]*JX+Y[i]*JY+Z[i]]=1;
			}
		break;
		case to_bm:
			for (i=0; i<n_particles; i++){
				PX=X[i];PY=Y[i];PZ=Z[i];
				X_bm.push_back(PX); Y_bm.push_back(PY); Z_bm.push_back(PZ);
			}
		break;
		case reset:
			for (i=0; i<n_particles; i++){
				PX=X_bm[i];PY=Y_bm[i];PZ=Z_bm[i];
				X.push_back(PX); Y.push_back(PY); Z.push_back(PZ);
			}
		break;
		default:
			success=false;
			cout <<"error in tranfer" << endl;
		break;
	}
	return success;
}

// Push outputs from this and other classes to Output class
void Teng::WriteOutput(int subloop){
      	PushOutput();
       New[0]->PushOutput();
      	for (int i = 0; i < n_out; i++) {
        	Out[i]->WriteOutput(subloop);
	}
}

// Additional information to pass outputs from Teng class.
void Teng::PushOutput() {
	int* point;
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
			Out[i]->PointerVectorReal.push_back(Mol[tag_mol]->gn); //this is the pointer to the start of the 'vector' that is reported to output.
			Out[i]->SizeVectorReal.push_back(sizeof(Mol[tag_mol]->gn)); //this is the size of the 'vector' that is reported to output
		}
		if (Out[i]->name=="ana" || Out[i]->name=="pos") { //example for putting 3 array's of Integers of arbitrary length to output
			string s="array;0";
			Out[i]->push("X",s);
			point=X.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(X.size());
			s="array;1";
			Out[i]->push("Y",s);
			point = Y.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(Y.size());
			s="array;2";
			Out[i]->push("Z",s);
			point=Z.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(Z.size());
		}

	}
}


// Standard procedures.
void Teng::PutParameter(string new_param) {
	KEYS.push_back(new_param);
}

string Teng::GetValue(string parameter) {
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

// Procedure that acquires the given inputs and checks
// Primary call to all engines (montecarlo or langevin dynamics) are also located here.
bool Teng::CheckInput(int start) {
  	if (debug) cout << "CheckInput in Teng" << endl;
  	bool success = true;

  	success = In[0]->CheckParameters("teng", name, start, KEYS, PARAMETERS, VALUES);
  	if (success) {
    		vector<string> options;
    		if (GetValue("MCS").size() > 0) {
      			success = In[0]->Get_int(GetValue("MCS"), MCS, 1, 10000, "The number of timesteps should be between 1 and 10000");
    		}
    		if (debug)
      			cout << "MCS is " << MCS << endl;

    		if (GetValue("save_interval").size() > 0) {
      			success = In[0]->Get_int(GetValue("save_interval"), save_interval,1,MCS,"The save interval nr should be between 1 and 100");
    		}
    		if (debug) cout << "Save_interval " << save_interval << endl;
    		if (Sys[0]->SysTagList.size() <1) {cout <<"Teng needs to have tagged molecules in the system" << endl; success=false;}
		else {tag_seg=Sys[0]->SysTagList[0]; if (Sys[0]->SysTagList.size()>1) {success=false; cout <<"Currently the Tagging is limited to one molecule per system. " << endl; }}
    		if (success) {
			n_particles = Seg[tag_seg]->n_pos;
    		}
    		tag_mol=-1;
    		int length = In[0]->MolList.size();
    		for (int i=0; i<length; i++) if (Mol[i]->freedom =="tagged") tag_mol=i;
  	}
  	if (success) {
   		n_out = In[0]->OutputList.size();
    		if (n_out == 0) cout << "Warning: no output defined!" << endl;
    		for (int i =  0; i < n_out; i++) {
      			Out.push_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[i], i, n_out));
       			if (!Out[i]->CheckInput(start)) {
        		cout << "input_error in output " << endl;
        		success=false;
      			}
     		}
     		//TODO: There should be a for loop over the required number of timesteps.
     		MonteCarlo();
   	}
  return success;
}
