#include "engine.h"

Engine::Engine(vector<Input*> In_,vector<System*> Sys_,vector<Newton*> New_, string name_) {
	In=In_; name=name_;   Sys=Sys_; New=New_; 
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

bool Engine::CheckInput() {
	bool success=true;
	success=In[0]->CheckParameters("engine",name,KEYS,PARAMETERS,VALUES);
	if (success) {
		vector<string> options;
		options.push_back("sfbox"); //options.push_back("zero-tension");
		if (GetValue("brand").size()>0) {
			if (!In[0]->Get_string(GetValue("sfbox"),brand,options,"In engine "+ name + "value of brand not recognised; default 'sfbox' used. ")) brand="sfbox"; 
		} else brand="sfbox"; 
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

void Engine::push(string s, double X) {
	doubles.push_back(s);
	doubles_value.push_back(X); 
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
	doubles.clear();
	doubles_value.clear();
	ints.clear();
	ints_value.clear();  
}
double* Engine::GetPointer(string s) {
	//vector<string> sub;
	//nothing yet
	return NULL;
}

int Engine::GetValue(string prop,int &int_result,double &double_result,string &string_result){
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
	length = doubles.size();
	while (i<length) {
		if (prop==doubles[i]) { 
			double_result=doubles_value[i];
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


bool Engine::Doit(){
	bool success=true;
	if (brand=="sfbox") {
		New[0]->AllocateMemory(); 
		New[0]->Solve();
	} else {
		cout <<"Sorry to date only 'sfbox' is implemented" << endl; 
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
