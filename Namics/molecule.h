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
	Molecule(vector<Input*>,vector<Lattice*>,vector<Segment*>,string); 
~Molecule();
  
	string name;  
	vector<Input*> In;
	vector<Segment*> Seg; 
	vector<Lattice*> Lat;
	vector<Alias*> Al;
	vector<int> MolMonList;
	vector<int> MolAlList;
	int start;
	int n_mol; 
	int mol_nr; 
	Real Mu; 
	Real theta;
	Real phibulk;
	Real fraction(int); 
	string freedom;
	MoleculeType MolType; 
	Real n; 
	Real GN,GN1,GN2; 
	Real norm;
	int chainlength;
	bool save_memory; 
	bool compute_phi_alias; 
	string composition;
	vector<int> Gnr;
	vector<int> first_s;
	vector<int> last_s;
	vector<int> first_b;
	vector<int> last_b;
	vector<int> mon_nr;
	vector<int> n_mon; 
	vector<int> molmon_nr; 
	vector<int> memory;
	vector<int> last_stored;
	Real *phi;
	Real *H_phi;
	Real *phitot;
	Real *H_phitot; 
	Real *Gg_f;
	Real *Gg_b; 
	Real *Gs; 
	Real *UNITY;
	int tag_segment; 

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
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string);
	bool ExpandAlias(vector<string>,string&);
	bool ExpandBrackets(string&);
	bool Interpret(string,int);
	bool GenerateTree(string,int,int&,vector<int>,vector<int>);
	bool Decomposition(string); 
	int GetChainlength(void); 
	Real Theta(void);
	string GetValue(string); 
	int GetMonNr(string);
	int GetAlNr(string);
	bool MakeMonList(void);  
	bool IsPinned(void);
	bool IsTagged(void); 
	bool IsCharged(void); 
	void AllocateMemory(void);
	bool PrepareForCalculations(int*); 
	bool ComputePhi(); 
	void propagate_forward(Real*, Real*, int&, int,int);
	Real* propagate_forward(int&,int,int);
	void propagate_backward(Real*, Real*, Real*,int&,int,int);
	void propagate_backward(int&,int,int);
	Real* Forward(int, int&);
	void Backward(Real*, int, int&);
	bool ComputePhiMon();
	bool ComputePhiLin();
	bool ComputePhiBra();
};

#endif
