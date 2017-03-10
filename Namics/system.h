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
	Real* CHI;
	vector<Segment*> Seg; 
	vector<Molecule*> Mol;
	vector<Lattice*> Lat; 
	vector<int> SysMonList; 
	vector<int> FrozenList;
	vector<int> SysTagList; 
	Real FreeEnergy;
	Real GrandPotential;
	Real* phitot;  
	int* KSAM;
	Real* H_GrandPotentialDensity;
	Real* H_FreeEnergyDensity;
	Real* H_alpha;
	Real* GrandPotentialDensity;
	Real* FreeEnergyDensity;
	Real* alpha;
	Real* TEMP;
	bool GPU; 
	int n_mol; 
	int solvent; 
	int tag_segment; 
	bool input_error; 
	bool cuda; 

	int MonA,MonB; 

	vector<string> ints;
	vector<string> Reals;
	vector<string> bools;
	vector<string> strings;
	vector<Real> Reals_value;
	vector<int> ints_value;
	vector<bool> bools_value;
	vector<string> strings_value;
	void push(string,Real);
	void push(string,int);
	void push(string,bool);
	void push(string,string);
	void PushOutput();
	Real* GetPointer(string);
	int GetValue(string,int&,Real&,string&); 
	string CalculationType; // {micro_emulsion,micro_phasesegregation};
	string GuessType; // {lamellae,Im3m,FCC,BCC,HEX,gyroid,Real_gyroid,Real_diamond,perforated_lamellae};
		
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
	Real GetFreeEnergy();
	Real GetGrandPotential();
	bool CreateMu();	
};
#endif

 
