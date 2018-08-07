#include "reaction.h"

Reaction::Reaction(vector<Input*> In_,vector<State*> Sta_, string name_) {
	In=In_; name=name_;   Sta=Sta_; 
	KEYS.push_back("K"); 
	KEYS.push_back("pK");
	KEYS.push_back("equation");  
}
Reaction::~Reaction() {
	DeAllocateMemory();
}
void Reaction::DeAllocateMemory(){
if (debug) cout <<"Destructor for Reaction " + name << endl; 

#ifdef CUDA

#else
#endif	
}

void Reaction::AllocateMemory(int Clamp_nr, int n_box) {
if (debug) cout <<"AllocateMemory in Reaction " + name << endl;
#ifdef CUDA
#else 
#endif
	
}

void Reaction::PrepareForCalculations() {
if (debug) cout <<"PrepareForCalculations in Reaction " + name << endl; 

}


void Reaction::PutParameter(string new_param) {
if (debug) cout <<"PutParameter in Reaction " + name << endl;
	KEYS.push_back(new_param); 
}

bool Reaction::CheckInput(int start) {
if (debug) cout <<"CheckInput in Reaction " + name << endl;
	bool success=true;
	K=-1;
	pK=-100;
	Sto.clear();
	State_nr.clear();
	success= In[0]->CheckParameters("reaction",name,start,KEYS,PARAMETERS,VALUES);
	if (success) {
		if (GetValue("K").size()==0 && GetValue("pK").size()==0)  {
			cout <<" reaction " << name << " has no K nor pK value" << endl; success=false;
		} else {
			if (GetValue("K").size() ==0) {
				pK=In[0]->Get_Real(GetValue("pK"),pK);
				if (pK==-100) {cout <<" reaction " << name << " no valid pK value found " << endl; success=false; }
				K=pow(10,-pK); 
			} else {
				K=In[0]->Get_Real(GetValue("K"),pK);
				if (K<0) {
					cout <<" reaction " << name << " has not a positive value for 'K' " << endl; success=false;
				} else {
					pK=-log10(K); 
				}
			}
		}
		equation=GetValue("equation"); 
		if (equation.size()==0) {
			success=false; cout <<" reaction " << name << " has no equation specified" <<endl; 
		} else {
			string s=equation; 
			vector<string>sub_equal;
			In[0]->split(s,'=',sub_equal);
			if (sub_equal.size()!=2) {
				cout <<" reaction : " << name << " equation : " << equation << "should have one '=' sign" << endl; 
				success=false; 
			} else {
				for (int k=0; k<2; k++) {
					vector<string>sub_plus;
					In[0]->split(sub_equal[k],'+',sub_plus);
					int sub_l=sub_plus.size();
					for (int l=0; l<sub_l; l++) { 
						vector<int>open;
						vector<int>close;
						In[0]->EvenBrackets(sub_plus[l], open, close);
						int length=open.size();
						if (length !=1) {
							cout <<" reaction : " << name << " equation " << equation << " has too many mon types in between '+' signs " << endl; 
							success=false; 
						} else  {
							string state_name=sub_plus[l].substr(open[0]+1,close[0]-open[0]-1);
							int num_states=In[0]->StateList.size(); 
							bool found=false;
							for (int i=0; i<num_states; i++) {
								string s_name=In[0]->StateList[i];
								if (state_name == s_name) {
									found = true; 
									State_nr.push_back(i); 
									Sta[i]->in_reaction=true;
									Seg_nr.push_back(Sta[i]->mon_nr);
								}
							}
							if (!found) {cout << " reaction : " << name << " equation " << equation << " state " << state_name << " not found " << endl; success=false;  }
							
							int sto=In[0]->Get_int(sub_plus[l].substr(0,open[0]),0);
							if (sto<1) {
								if (sto==0) cout << " reaction : " << name << " equation : " << equation << " has a zero as stocheometry number " << endl;
								else
								cout << " reaction : " << name << " equation : " << equation << " has negative stocheometry number " << endl; 
								success=false; 
							}
							if (k>0) sto*=-1;
							Sto.push_back(sto);
						}
					}
				}
				
			}

		} 
	}
	int length=In[0]->MonList.size();
	int LENGTH=Sto.size();  
	for (int i=0; i<length; i++) {
		int sum=0;
		
		for (int j=0; j<LENGTH; j++) {
			if (Seg_nr[j]==i) sum+=Sto[j]; 
		}
		if (sum !=0 ) {
			success=false;
			cout <<" reaction : " << name << " equation : " << equation << " not balanced for internal states of mon type: " << In[0]->MonList[i] << endl; 
		}
	}
	Real charge=0;
	for (int j=0; j<LENGTH; j++) {
		charge += Sto[j]*Sta[State_nr[j]]->valence;
	}
	if (abs(charge)>1e-10) {
		success=false;
		cout <<" reaction : " << name << " equation : " << equation << " not electroneutral. Absolute value appears : " << abs(charge) << endl; 
	}
	
	return success;
}
 
string Reaction::GetValue(string parameter){
if (debug) cout <<"GetValue in Reaction " + name << endl;
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
 
void Reaction::push(string s, Real X) {
if (debug) cout <<"push (Real) in Reaction " + name << endl;
	Reals.push_back(s);
	Reals_value.push_back(X); 
}
void Reaction::push(string s, int X) {
if (debug) cout <<"push (int) in Reaction " + name << endl;
	ints.push_back(s);
	ints_value.push_back(X); 
}
void Reaction::push(string s, bool X) {
if (debug) cout <<"push (boool) in Reaction " + name << endl;
	bools.push_back(s);
	bools_value.push_back(X); 
}
void Reaction::push(string s, string X) {
if (debug) cout <<"push (string) in Reaction " + name << endl;
	strings.push_back(s);
	strings_value.push_back(X); 	
}
void Reaction::PushOutput() {
if (debug) cout <<"PushOutput in Reaction " + name << endl;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();  
	push("equation",equation);
	push("pK",pK);
#ifdef CUDA

#endif
}

Real* Reaction::GetPointer(string s) {
if (debug) cout <<"GetPointer in Reaction " + name << endl;
	return NULL;
}


int Reaction::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
if (debug) cout <<"GetValue (long)  in Reaction " + name << endl;
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

Real Reaction::Residual_value() {
	Real res_value=0;
	Real alphab=0;
	int length=Sto.size();
	for (int i=0; i<length; i++) {
		alphab=Seg[Sta_nr[i]]Sta[State_nr[i]]->alphabulk;
		if (alphab<0 || alphab>1) cout <<"alphabulk out of range: programming error; returning undetermined value to residuals: do not continue" << endl; 
		else res_value+=Sto[i]*log10(alphab); 
	}
	return res_value+pK;
}
