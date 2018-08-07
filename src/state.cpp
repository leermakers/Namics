#include "state.h"

State::State(vector<Input*> In_,vector<Segment*> Seg_, string name_) {
	In=In_; name=name_;  Seg=Seg_;
	KEYS.push_back("alphabulk");
	KEYS.push_back("valence");
	KEYS.push_back("mon");
}
State::~State() {
	DeAllocateMemory();
}
void State::DeAllocateMemory(){
if (debug) cout <<"Destructor for State " + name << endl; 
#ifdef CUDA
#else
#endif	
}

void State::AllocateMemory(int Clamp_nr, int n_box) {
if (debug) cout <<"AllocateMemory in State " + name << endl;

#ifdef CUDA
#else 
#endif	
}

void State::PrepareForCalculations() {
if (debug) cout <<"PrepareForCalculations in State " + name << endl; 

}


void State::PutParameter(string new_param) {
if (debug) cout <<"PutParameter in State " + name << endl;
	KEYS.push_back(new_param); 
}

bool State::CheckInput(int start) {
if (debug) cout <<"CheckInput in State " + name << endl;
	bool success=true;
	fixed=false;
	in_reaction =false;
	alphabulk =-1;
	valence=0;
	success= In[0]->CheckParameters("state",name,start,KEYS,PARAMETERS,VALUES);
	int length=In[0]->MonList.size();
	for (int k=0; k<length; k++) {
		if (Seg[k]->name == name) { cout << "name of state can not be the same as the name of any mon in the system" << endl; 
		success=false;
		}
	}
	if (success) {
		if (GetValue("mon").size()==0) {
			cout << "Please specify the 'mon' for state " << name << endl; 
			success=false;
		} else {
			string mon_name;
			mon_name=GetValue("mon"); 
			bool mon_found=false;
			for (int k=0; k<length; k++) {
				if (Seg[k]->name==mon_name) {
					mon_found=true; 
					mon_nr=k;
					
				} 
			}
			if (!mon_found) {
				cout <<"in state " << name << " mon name " << mon_name << " not found in system " << endl; 
				success=false;
			}
			if (GetValue("alphabulk").size()>0) {
				fixed=true;
				alphabulk=In[0]->Get_Real(GetValue("alphabulk"),alphabulk);
				if (alphabulk <0 || alphabulk > 1) {
					cout << "for state " << name << " value for alphabulk is out of range 0 ... 1 " << endl; 
					success=false;
				}
			}
			if (GetValue("valence").size()>0) {
				valence=In[0]->Get_Real(GetValue("valence"),valence);
				if (valence <-10 || valence > 10) {
					cout << "for state " << name << " value for valence " << GetValue("valence") << " is out of range -10 ... 10 " << endl; 
					success=false;
				}
			}
		} 
		state_nr=Seg[mon_nr]->AddState(name,alphabulk,valence);
	}
	return success;
}
 
string State::GetValue(string parameter){
if (debug) cout <<"GetValue in State " + name << endl;
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
 
void State::push(string s, Real X) {
if (debug) cout <<"push (Real) in State " + name << endl;
	Reals.push_back(s);
	Reals_value.push_back(X); 
}
void State::push(string s, int X) {
if (debug) cout <<"push (int) in State " + name << endl;
	ints.push_back(s);
	ints_value.push_back(X); 
}
void State::push(string s, bool X) {
if (debug) cout <<"push (boool) in State " + name << endl;
	bools.push_back(s);
	bools_value.push_back(X); 
}
void State::push(string s, string X) {
if (debug) cout <<"push (string) in State " + name << endl;
	strings.push_back(s);
	strings_value.push_back(X); 	
}
void State::PushOutput() {
if (debug) cout <<"PushOutput in State " + name << endl;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();  
#ifdef CUDA
#endif
}

Real* State::GetPointer(string s) {
if (debug) cout <<"GetPointer in State " + name << endl;
	return NULL;
}


int State::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
if (debug) cout <<"GetValue (long)  in State " + name << endl;
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
