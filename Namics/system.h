#ifndef SYSTEMxH
#define SYSTEMxH
#include "namics.h"
#include "segment.h"
#include "molecule.h"
#include "lattice.h"
class System {
public:
	System(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,string);

~System();

	string name;
	vector<Input*> In; 
	double* CHI;
	vector<Segment*> Seg; 
	vector<Molecule*> Mol;
	vector<Lattice*> Lat; 
	vector<int> SysMonList; 
	vector<int> FrozenList;
	vector<int> SysTagList; 
	double FreeEnergy;
	double GrandPotential;
	double* phitot;  
	int* KSAM;
	double* H_GrandPotentialDensity;
	double* H_FreeEnergyDensity;
	double* H_alpha;
	double* GrandPotentialDensity;
	double* FreeEnergyDensity;
	double* alpha;
	double* TEMP;
	bool GPU; 
	int n_mol; 
	int solvent; 
	int tag_segment; 
	bool input_error; 
	bool cuda; 

	int MonA,MonB; 

	vector<string> ints;
	vector<string> doubles;
	vector<string> bools;
	vector<string> strings;
	vector<double> doubles_value;
	vector<int> ints_value;
	vector<bool> bools_value;
	vector<string> strings_value;
	void push(string,double);
	void push(string,int);
	void push(string,bool);
	void push(string,string);
	void PushOutput();
	double* GetPointer(string);
	int GetValue(string,int&,double&,string&); 
	string CalculationType; // {micro_emulsion,micro_phasesegregation};
	string GuessType; // {lamellae,Im3m,FCC,BCC,HEX,gyroid,double_gyroid,double_diamond,perforated_lamellae};
		
	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string); 
	string GetValue(string); 	
	string GetMonName(int );
	bool CheckChi_values(int);
	void AllocateMemory();
	bool PrepareForCalculations();
	bool ComputePhis();
	bool CheckResults();
	double GetFreeEnergy();
	double GetGrandPotential();
	bool CreateMu();	
};
#endif

 
