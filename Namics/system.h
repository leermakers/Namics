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
	int M; 
	int MX,MY,MZ;
	vector<int> SysMonList; 
	vector<int> FrozenList;
	vector<int> SysTagList; 
	double FreeEnergy;
	double GrandPotential;
	double* phitot; 
	double* phib; 
	double* KSAM;
	double* H_KSAM;
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
	string calculation_type; 
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput();
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

 
