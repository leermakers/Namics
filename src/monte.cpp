#include "monte.h"

Monte::Monte(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,vector<System*> Sys_,vector<Newton*> New_, string name_) {
	cout << "constructor in Monte" << endl;
	In=In_; name=name_;   Lat=Lat_; Mol=Mol_; Seg = Seg_; Sys=Sys_; New=New_; 
	KEYS.push_back("timesteps"); KEYS.push_back("walks"); 
	KEYS.push_back("timebetweensaves");
}
Monte::~Monte() {
}

void Monte::AllocateMemory() {
	if (debug) cout <<"nothing to allocate in Monte" << endl; 
}


void Monte::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

bool Monte::CheckInput(int start) {
if (debug) cout <<"Check Monte" << endl;
	bool success=true;
	success=In[0]->CheckParameters("monte",name,start,KEYS,PARAMETERS,VALUES);
	if (success) {
		vector<string> options;
		if (GetValue("timesteps").size()>0) {
                   success=In[0]->Get_int(GetValue("timesteps"),timesteps, 1, 10000, "The number of timesteps should be between 1 and 10000");
		} 
	cout <<"timesteps is " << timesteps << endl;

	}
	return success; 
}

int Monte::Positions1Dto3D(int* H_P, int n_pos, int JX, int JY){
// Usually position is flattened into 1D. To make monte carlo moves we need particle co-ordinates in 3D
// This procedure converts flattend array into expanded 3D array. 
/*	int pos_x[n_pos];
	int pos_y[n_pos];
	int pos_z[n_pos];
	for(int i=0; i<n_pos; i++){
	pos_x[i]=H_P[i]/JX;
	pos_y[i]=(H_P[i]-pos_x[i]*JX)/JY;
	pos_z[i]=(H_P[i]-pos_x[i]*JX)%JY;
	}
*/
	return 0;
}

bool Monte::Simulate(){
	bool success=true;
	for (int i=1; i<3; i++){
	cout << "Running Monte-Carlo simulations" << endl;
	New[0]->Solve(true);
	cout << "Change in Free energy, computed" << endl;}
	return success;
}








string Monte::GetValue(string parameter){
	int i=0;
	int length = PARAMETERS.size();
	while (i<length) {
		if (parameter==PARAMETERS[i]) {
			return VALUES[i]; 
		}
		i++;
	}
	return "" ; 
}

void Monte::push(string s, Real X) {
	Reals.push_back(s);
	Reals_value.push_back(X); 
}
void Monte::push(string s, int X) {
	ints.push_back(s);
	ints_value.push_back(X); 
}
void Monte::push(string s, bool X) {
	bools.push_back(s);
	bools_value.push_back(X); 
}
void Monte::push(string s, string X) {
	strings.push_back(s);
	strings_value.push_back(X); 	
}
void Monte::PushOutput() {
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();  
}
Real* Monte::GetPointer(string s) {
	//vector<string> sub;
	//nothing yet
	return NULL;
}

int Monte::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
	int i=0;
	int length = ints.size();
	while (i<length) {
		if (prop==ints[i]) { 
			int_result=ints_value[i];
			return 1;
		}
		i++;
	}
	i=0;
	length = Reals.size();
	while (i<length) {
		if (prop==Reals[i]) { 
			Real_result=Reals_value[i];
			return 2;
		}
		i++;
	}
	i=0;
	length = bools.size();
	while (i<length) {
		if (prop==bools[i]) { 
			if (bools_value[i]) string_result="true"; else string_result="false"; 
			return 3;
		}
		i++;
	}
	i=0;
	length = strings.size();
	while (i<length) {
		if (prop==strings[i]) { 
			string_result=strings_value[i]; 
			return 3;
		}
		i++;
	}
	return 0; 
}
