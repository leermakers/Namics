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
	rene Mu; 
	rene theta;
	rene phibulk;
	rene fraction(int); 
	string freedom;
	MoleculeType MolType; 
	rene n; 
	rene GN,GN1,GN2; 
	rene norm;
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
	rene *phi;
	rene *H_phi;
	rene *phitot;
	rene *H_phitot; 
	rene *Gg_f;
	rene *Gg_b; 
	rene *Gs; 
	rene *UNITY;
	int tag_segment; 

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
	rene Theta(void);
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
	void propagate_forward(rene*, rene*, int&, int,int);
	rene* propagate_forward(int&,int,int);
	void propagate_backward(rene*, rene*, rene*,int&,int,int);
	void propagate_backward(int&,int,int);
	rene* Forward(int, int&);
	void Backward(rene*, int, int&);
	bool ComputePhiMon();
	bool ComputePhiLin();
	bool ComputePhiBra();
};

#endif
