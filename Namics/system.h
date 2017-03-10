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
	rene* CHI;
	vector<Segment*> Seg; 
	vector<Molecule*> Mol;
	vector<Lattice*> Lat; 
	vector<int> SysMonList; 
	vector<int> FrozenList;
	vector<int> SysTagList; 
	rene FreeEnergy;
	rene GrandPotential;
	rene* phitot;  
	int* KSAM;
	rene* H_GrandPotentialDensity;
	rene* H_FreeEnergyDensity;
	rene* H_alpha;
	rene* GrandPotentialDensity;
	rene* FreeEnergyDensity;
	rene* alpha;
	rene* TEMP;
	bool GPU; 
	int n_mol; 
	int solvent; 
	int tag_segment; 
	bool input_error; 
	bool cuda; 

	int MonA,MonB; 

	vector<string> ints;
	vector<string> renes;
	vector<string> bools;
	vector<string> strings;
	vector<rene> renes_value;
	vector<int> ints_value;
	vector<bool> bools_value;
	vector<string> strings_value;
	void push(string,rene);
	void push(string,int);
	void push(string,bool);
	void push(string,string);
	void PushOutput();
	rene* GetPointer(string);
	int GetValue(string,int&,rene&,string&); 
	string CalculationType; // {micro_emulsion,micro_phasesegregation};
	string GuessType; // {lamellae,Im3m,FCC,BCC,HEX,gyroid,rene_gyroid,rene_diamond,perforated_lamellae};
		
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
	rene GetFreeEnergy();
	rene GetGrandPotential();
	bool CreateMu();	
};
#endif

 
