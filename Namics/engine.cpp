#include "engine.h"

Engine::Engine(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,vector<System*> Sys_,vector<Newton*> New_, string name_) {
	In=In_; name=name_;   Lat=Lat_; Mol=Mol_; Seg = Seg_; Sys=Sys_; New=New_; 
	KEYS.push_back("MC"); KEYS.push_back("SN_MD"); KEYS.push_back("brand"); 
}
Engine::~Engine() {
}

void Engine::AllocateMemory() {
	cout <<"nothing to allocate in Engine" << endl; 
}


void Engine::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

bool Engine::CheckInput(int start) {
	bool success=true;
//right now all the checks for engine are done in input.cpp. could be move back here.
	success=In[0]->CheckParameters("engine",name,start,KEYS,PARAMETERS,VALUES);
	if (success) {
		vector<string> options;
		options.push_back("sfbox"); options.push_back("var"); options.push_back("search");
		if (GetValue("brand").size()>0) {
                        if (!In[0]->Get_string(GetValue("brand"),brand,options,"In engine " + name + " value of brand " + brand + " is not recognised"));
		} else brand="sfbox";
		
		if(brand=="search"){
			success=In[0]->LoadItems(brand,VAR_key,VAR_param,VAR_val);	
			int nmol=In[0]->MolList.size(); bool molfound=false;
			for(int j=0; j<nmol; j++) {
				if(VAR_key[0]=="mol" && VAR_param[0]==In[0]->MolList[j]) {
					cout << "For brand 'search', search will be performed over '"+VAR_param[0]+"'. "+VAR_val[0]+" will be modified."<< endl;	
					molfound=true; 
				}
			}
			if(VAR_key[0]=="mol" && !molfound) {cout << "Error. Search undefined. Execution stopped." << endl; success=false; return 0;}
			if(VAR_key[1]=="mol" && VAR_param[1]!="phibulk") {cout << "Target for 'mol' can only be 'phibulk'."<< endl; success=false;}
			if(VAR_key[1]=="sys" && VAR_param[1]!="GrandPotential") {cout << "Target for 'sys' can only be 'GrandPotential as of now." << endl; success=false;}
		}
	}
	return success; 
}
 
string Engine::GetValue(string parameter){
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

void Engine::push(string s, Real X) {
	Reals.push_back(s);
	Reals_value.push_back(X); 
}
void Engine::push(string s, int X) {
	ints.push_back(s);
	ints_value.push_back(X); 
}
void Engine::push(string s, bool X) {
	bools.push_back(s);
	bools_value.push_back(X); 
}
void Engine::push(string s, string X) {
	strings.push_back(s);
	strings_value.push_back(X); 	
}
void Engine::PushOutput() {
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();  
}
Real* Engine::GetPointer(string s) {
	//vector<string> sub;
	//nothing yet
	return NULL;
}

int Engine::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
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

int Engine::SubProblemNum(){
// Apart from the determination of number of sub problems this routine also
// checks the input values if the engine brand var is chosen.
// This could be placed in checkinput of engine, but right now I opt to place it here. - Ram.

	bool success = true;
	bool varkeyfound = false;
	int len = 0;

	if(brand=="var"){
		success = In[0]->LoadItems(brand, VAR_key, VAR_param, VAR_val);
		vector<string> varkeys;
		varkeys.push_back("mol"); varkeys.push_back("sys"); varkeys.push_back("mon");
		int keylength = VAR_key.size();
		for(int i=0; i<keylength-1; i++) {
			if(VAR_key[i+1] != VAR_key[i]){
				cout << "Only one parameter can be varied. Please do not vary mulitple parameters at the same time." << endl; success= false; return 0;
			}
		}
		if (success){
			int varkeylength = varkeys.size();
			for (int i=0; i<varkeylength; i++){
				if(varkeys[i] == VAR_key[0]) {varkeyfound = true; break;
				} else { varkeyfound = false;
				}
			}
			if (!varkeyfound){
				cout << "Can't vary anything in '" + VAR_key[0] + "'. Choose from:" << endl;
				for (int i=0; i<varkeylength; i++){
				cout << i << ": " << varkeys[i] << endl;
				} success =false;
			}
		}

// Checking the parameters provided against the parameters that can be varied. 
// Also prints error if the paramter provided is not a standard parameter.
		if (success) {
			if (VAR_key[0] == "mol") {
			vector<string> var_param;
			var_param.push_back("theta"); var_param.push_back("n"); var_param.push_back("phibulk"); var_param.push_back("chainlength");
			bool parameterfound=false; int paramlength = var_param.size();
				for(int i=0; i<paramlength; i++){
					if(VAR_val[0] == var_param[i]) { parameterfound=true; break;
					} else { parameterfound=false; 
					}
				}
				if(!parameterfound){
				cout << "Parmeter '" + VAR_val[0] + "' provided can't be varied. Choose from:" << endl;
					for(int i=0; i<paramlength; i++){
						cout << i << ": " << var_param[i] << endl;
						success=false;
					}
				}
			int nmol = In[0]->MolList.size(); bool molfound=false;
			for (int i=0; i<nmol; i++) if(VAR_param[0]==In[0]->MolList[i]) molfound=true;
			if(!molfound) {cout << "Molecule entered in 'var' statement doesn't exist. May be a typo! Can you change it? " << endl; success=false;}
			vector<string> loopers; loopers.push_back("start"); loopers.push_back("end"); loopers.push_back("step");
			bool loopformat=false;
			if(VAR_param[1]=="start" && VAR_param[2]=="end" && VAR_param[3]=="step") loopformat=true;
			if(!loopformat) {cout << "Looping in var has a format. Refer authors or documentation. I see a violation of format." << endl; success=false;}
		

// Checks should be implemented to notify unreasonable values given as starting and ending values.
// Omitted right now. And the code doesn't check them.
			}
                        if (VAR_key[0] == "sys") {
                        vector<string> var_param;
                        var_param.push_back("null"); 
                        bool parameterfound=false; int paramlength = var_param.size();
                                for(int i=0; i<paramlength; i++){
                                        if(VAR_val[0] == var_param[i]) { parameterfound=true; break;
                                        } else { parameterfound=false;
                                        }
                                }
                                if(!parameterfound){
                                cout << "Parmeter '" + VAR_val[0] + "' provided can't be varied. Choose from:" << endl;
                                        for(int i=0; i<paramlength; i++){
                                                cout << i << ": " << var_param[i] << endl;
                                                success=false;
                                        }
                                }
                        }
                        if (VAR_key[0] == "mon") {
                        vector<string> var_param;
			int n_seg = In[0]->MonList.size();
			
			for (int i=0; i<n_seg; i++) {
                        	var_param.push_back("chi-"+Seg[i]->name);
			}
 
                        bool parameterfound=false; int paramlength = var_param.size();
                                for(int i=0; i<paramlength; i++){
                                        if(VAR_val[0] == var_param[i]) { parameterfound=true; break;
                                        } else { parameterfound=false;
                                        }
                                }
                                if(!parameterfound){
                                cout << "Parmeter '" + VAR_val[0] + "' provided can't be varied. Choose from:" << endl;
                                        for(int i=0; i<paramlength; i++){
                                                cout << i << ": " << var_param[i] << endl;
                                                success=false;
                                        }
                                }
                        }
			
		}
		if (success) {
                        
                        Real start = In[0]->Get_Real(VAR_val[1],10*Lat[0]->volume);
                        Real end = In[0]->Get_Real(VAR_val[2],10*Lat[0]->volume);
                        Real step = In[0]->Get_Real(VAR_val[3],10*Lat[0]->volume);
                        len = (end-start)/step +1;
			 
		}
	}
	else len=1;
	return len;
}

bool Engine::Doit(int sub){
	bool success=true;
	if (brand=="sfbox") {
		New[0]->AllocateMemory();
		New[0]->Guess();
		New[0]->Solve();
	} else if (brand == "var"){
		if (success) {
                        int i=sub;
			New[0]->AllocateMemory();
                        if(VAR_key[0]=="mol") success=VarMol(i); 
			if(VAR_key[0]=="mon") success=VarMon(i);
			if(!success) {cout << "failure" <<endl; return 0;}
	                New[0]->Solve(); i++; 
                } else {
                cout << "There is a problem in loading the 'var' parameters" << endl;  return 0;
                }

	} else if (brand == "search") {
		//search just searches, doesn't print output anything during loops. 
		//search procedure is called n times until success, if n exceeds 2000 problem is stopped.
		int searchlimit=2000; int i=0; 
		success = In[0]->LoadItems(brand,VAR_key,VAR_param,VAR_val);
		int nmol=In[0]->MolList.size();int id=0;
		for(int i=0; i<nmol; i++){if(VAR_param[0]==In[0]->MolList[i]) {id=i;}}
		New[0]->AllocateMemory();
		while (i<searchlimit) {
			if(VAR_param[1]=="GrandPotential"){
				success=Search(id);
				if(success) {cout << "search is successfull" << endl; break;}
			} else if (VAR_param[1]=="phibulk"){
				//To be implemented when necessary.
			}
			i++;
		}

	} else {
		cout <<"Sorry to date only 'sfbox' & 'var' only are implemented.." << endl; return 0; 
	}
	return success;
}

bool Engine::VarMol(int sub){
	bool success = true;
	int nmol=In[0]->MolList.size(); int i=sub; int j=0;
	while (j < nmol) {
		if(In[0]->MolList[j]==VAR_param[0]){
		Real start = In[0]->Get_Real(VAR_val[1],10*Lat[0]->volume);
		Real end = In[0]->Get_Real(VAR_val[2],10*Lat[0]->volume);
		Real step = In[0]->Get_Real(VAR_val[3],10*Lat[0]->volume);
		int length = (end-start)/step;
			if (VAR_val[0] == "theta") {
				Mol[j]->theta = start +(i)*step; Mol[j]->n=Mol[j]->theta/Mol[j]->chainlength; Lat[0]->GuessVar(New[0]->xx,Mol[j]->theta,Sys[0]->GuessType,Seg[Sys[0]->MonA]->guess_u,Seg[Sys[0]->MonB]->guess_u);
				cout << "Subproblem number: " << i << " out of " << length <<"; Paramter varied is 'theta' of Molecule: " <<In[0]->MolList[j] << "; current value is: "<< Mol[j]->theta << endl;
				success=true; break;
			} else if (VAR_val[0]=="n") {
                	       	Mol[j]->n = start +(i)*step; Mol[j]->theta=Mol[j]->n*Mol[j]->chainlength;New[0]->Guess();
				cout << "Subproblem number: " << i << " out of " << length <<"; Paramter varied is 'n' of Molecule: " <<In[0]->MolList[j] << "; current value is: "<< Mol[j]->n << endl; 
				success=true; break;
			} else if (VAR_val[0]=="chainlength") { 
        	               	Mol[j]->chainlength = start +(i)*step; Mol[j]->n=Mol[j]->theta/Mol[j]->chainlength;New[0]->Guess();
				cout << "Subproblem number: " << i << " out of " << length <<"; Paramter varied is 'chainlenth' of Molecule: " <<In[0]->MolList[j] << "; current value is: "<< Mol[j]->chainlength << endl;
				success=true; break;
			} else if (VAR_val[0] == "phibulk"){
        	               	Mol[j]->phibulk = start +(i)*step;New[0]->Guess();
				cout << "Subproblem number: " << i << " out of " << length <<"; Paramter varied is 'phibulk' of Molecule: " <<In[0]->MolList[j] << "; current value is: "<< Mol[j]->phibulk << endl;
				success=true;
			}
		} else {
		success=false;
		}j++;
	}
	return success;
} 

bool Engine::VarMon(int sub){
        bool success = true;
        int nseg=In[0]->MonList.size(); int i=sub; int j=0;
	vector<string> chimon;
	for (int k=0; k<nseg; k++) chimon.push_back("chi-"+Seg[k]->name);
        while (j < nseg) {
                if(chimon[j]==VAR_val[0]){
                Real start = In[0]->Get_Real(VAR_val[1],10*Lat[0]->volume);
		Real end = In[0]->Get_Real(VAR_val[2],10*Lat[0]->volume);
                Real step = In[0]->Get_Real(VAR_val[3],10*Lat[0]->volume);
		int length = (end-start)/step;
			for(int r=0; r<nseg; r++) if(r!=j) {Sys[0]->CHI[r*nseg+j]=start+i*step; Sys[0]->CHI[j*nseg+r]=Sys[0]->CHI[r*nseg+j];}
			New[0]->Guess();
			cout << "Subproblem number: " << i << " out of " << length <<"; Paramter varied is '" + VAR_val[0] + "' of Monomer: " <<VAR_param[0] << endl; 
			success=true;
			break;
                } else {
                success=false;
                }j++;
        }
	for(int r=0; r<nseg; r++) for(int s=0; s<nseg; s++) cout<< r*nseg+s <<" : " <<Sys[0]->CHI[r*nseg+s] << endl;
        return success;
}

bool Engine::Search(int id){
	bool success=true; int m=id;
	success = In[0]->LoadItems(brand,VAR_key,VAR_param,VAR_val);
	Real value = In[0]->Get_Real(VAR_val[1],10*Lat[0]->volume);
	Real tolerance = pow(10.0 ,-3.0);
	Real ll=value-tolerance; Real ul=value+tolerance;
	New[0]->Guess(); New[0]->Solve();
	Real lower; Real upper; Real mean;
	if(Sys[0]->GrandPotential < ll || Sys[0]->GrandPotential > ul) {
		Mol[m]->theta=Mol[m]->theta+0.1*Mol[m]->theta;	
		Mol[m]->n=Mol[m]->theta/Mol[m]->chainlength;
		cout << Sys[0]->GrandPotential << " : " << Mol[m]->theta << endl;
		success=false;
	} else {
		lower=Mol[m]->theta; upper=Mol[m]->theta+0.1*Mol[m]->theta;
		for(int i=0; i<100; i++){
			mean=(upper+lower)/2.0;Mol[m]->theta=mean;Mol[m]->n=Mol[m]->theta/Mol[m]->chainlength;
			New[0]->Guess(); New[0]->Solve();
			if(Sys[0]->GrandPotential > value) lower=mean; else upper=mean;
		}
		success=true;
	}
	
return success;
}

/*
void Engine::subdomain(int *fBx,int *fBy,int *fBz, int *fPx,int *fPy,int *fPz, int zloc){
    int loop = 0;
		fPx[loop] = 26 ; fPy[loop] = 26 ; fPz[loop] = 26 ;
		fBx[loop] = 26-10; fBy[loop] = 26-10 ; fBz[loop] = 26-10 ;
    for (int xcord=1; xcord<=10; xcord++){
		for(int ycord=1; ycord<=10; ycord++){	
		fPx[loop] = (xcord)*5-2 ; fPy[loop] = (ycord)*5-2 ; fPz[loop] = zloc ;
		fBx[loop] = (xcord*5-2)-10; fBy[loop] = (ycord*5-2)-10 ; fBz[loop] = zloc-10 ;	
		loop = loop+1;
		}
	}
	for (int diag=1; diag<=12; diag++){
		fPx[loop] = diag*4 ; fPy[loop] = diag*4 ; fPz[loop] = 26 ;
		fBx[loop] = (diag*4)-10; fBy[loop] = (diag*4)-10 ; fBz[loop] = 16 ;
		loop = loop+1;
	}
} 
*/
