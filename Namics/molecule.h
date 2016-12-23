#ifndef MOLECULExH
#define MOLECULExH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "alias.h"
#include "lattice.h"
#include "tools.h"
class Molecule {
public:
	Molecule(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Alias*>,string); 
~Molecule();
  
	string name;  
	vector<Input*> In;
	vector<Segment*> Seg; 
	vector<Lattice*> Lat;
	vector<Alias*> Al;
	vector<int> MolMonList;
	int n_mol; 
	int M; 
	int MX,MY,MZ;
	int BX1,BY1,BZ1,BXM,BYM,BZM;
	int JX;
	int JY; 	
	int mol_nr; 
	double Mu; 
	double theta;
	double phibulk;
	double fraction(int); 
	string freedom;
	MoleculeType MolType; 
	double n; 
	double GN,GN1,GN2; 
	int chainlength;
	bool save_memory; 
	string composition;
	vector<int> mon_nr;
	vector<int> n_mon; 
	vector<int> molmon_nr; 
	double *phi;
	double *H_phi;
	double *phitot;
	double *Gg_f;
	double *Gg_b; 
	int tag_segment; 
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(void);
	void PutParameter(string);
	bool Decomposition(string); 
	int GetChainlength(void); 
	double Theta(void);
	string GetValue(string); 
	int GetMonNr(string);
	bool MakeMonList(void);  
	bool IsPinned(void);
	bool IsTagged(void); 
	bool IsCharged(void); 
	void AllocateMemory(void);
	bool PrepareForCalculations(); 
	bool ComputePhi(); 
	bool ComputePhiMon();
	bool ComputePhiLin();

	
};

#endif
