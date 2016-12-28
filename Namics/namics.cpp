#define MAINxH
#include "namics.h"
#include "tools.h" 
#include "input.h"
#include "lattice.h"
#include "segment.h"
#include "alias.h"
#include "molecule.h"
#include "system.h"
#include "newton.h"
#include "engine.h"
#include "output.h"

string version="0.0.0.0.0.0.0.0";
//meaning:
// engine version number =0
// newton version number =0
// system version number =0
// lattice version number =0
// molecule version number =0
// segment version number =0
// alias version number =0
// output version number =0

double b_length=5e-10;
double e=1.60217e-19;
double T=298.15;
double k_B=1.38065e-23;
double k_BT=k_B*T;
double eps0=8.85418e-12;
double eps=80;
double factor=e*e/(eps*eps0*b_length*k_BT); 

int main(int argc, char *argv[]) {
	string fname;
	string filename;
	string line;
	ofstream out_file;
	bool cuda; 

#ifdef CUDA
	cudaDeviceReset();
	stat = cublasCreate(&handle); if (stat !=CUBLAS_STATUS_SUCCESS) {printf("CUBLAS failed \n");}
	cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
	BlasResult= (double*)AllOnDev(1);
	cuda=true;
#else 
	cuda=false; 
#endif
	if (argc == 2) fname = argv[1]; else {printf("Use: namics filename -without extension- \n"); return 1;}
	filename = fname + ".in";
	vector<Input*> In; In.push_back(new Input(filename)); if (In[0]->Input_error) {return 0;}
	vector<Lattice*> Lat; Lat.push_back(new Lattice(In,In[0]->LatList[0])); if (!Lat[0]->CheckInput()) {return 0;}
	vector<Segment*> Seg; int n_seg=In[0]->MonList.size();	
	for (int i=0; i<n_seg; i++) Seg.push_back(new Segment(In,Lat,In[0]->MonList[i],i,n_seg));
	for (int i=0; i<n_seg; i++) { for (int k=0; k<n_seg; k++) Seg[i]->PutChiKEY(Seg[k]->name); if (!Seg[i]->CheckInput()) return 0;}
	vector<Alias*> Al; int n_al = In[0]->AliasList.size();
	for (int i=0; i<n_al; i++) {Al.push_back(new Alias(In,In[0]->AliasList[i])); if (!Al[i]->CheckInput()) return 0;} 
	int n_mol = In[0]->MolList.size();  
	vector<Molecule*> Mol; for (int i=0; i<n_mol; i++) {Mol.push_back(new Molecule(In,Lat,Seg,Al,In[0]->MolList[i])); if (!Mol[i]->CheckInput()) return 0;}
	vector<System*> Sys; Sys.push_back(new System(In,Lat,Seg,Mol,In[0]->SysList[0])); Sys[0]->cuda=cuda;  
	if (!Sys[0]->CheckInput()) {return 0;} if (!Sys[0]->CheckChi_values(n_seg)) return 0; 
	vector<Newton*> New; New.push_back(new Newton(In,Lat,Seg,Mol,Sys,In[0]->NewtonList[0])); if (!New[0]->CheckInput()) {return 0;}
	vector<Engine*> Eng; Eng.push_back(new Engine(In,Sys,New,In[0]->EngineList[0])); if (!Eng[0]->CheckInput()) {return 0;} 
	int n_out = In[0]->OutputList.size(); 
	if (n_out ==0) cout <<"Warning: no output defined" << endl; 
	vector<Output*> Out; for (int i=0; i<n_out; i++) { Out.push_back(new Output(In,Lat,Seg,Mol,Sys,New,Al,Eng,In[0]->OutputList[i],i,n_out)); if (Out[i]->input_error) return 0;}

	//New[0]->AllocateMemory(); 
	//New[0]->PrepareForCalculations(); 
	//New[0]->Solve(); 
	
 	Eng[0]->Doit(); 
	for (int i=0; i<n_out; i++) Out[i]->WriteOutput(); 	

//	Out[0]->density();
//	Out[0]->printlist();
	
	return 0;
}

