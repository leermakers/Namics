#include "engine.h"

Engine::Engine(vector<Input*> In_,vector<System*> Sys_, string name_) {
	In=In_; name=name_;   Sys=Sys_; 
	KEYS.push_back("MC"); KEYS.push_back("SN_MD");
}
Engine::~Engine() {
}
void Engine::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

bool Engine::CheckInput() {
	return In[0]->CheckParameters("engine",name,KEYS,PARAMETERS,VALUES);
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
 

/*
void subdomain(int *fBx,int *fBy,int *fBz, int *fPx,int *fPy,int *fPz, int zloc){
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