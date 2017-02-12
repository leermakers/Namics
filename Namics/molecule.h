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
	double Mu; 
	double theta;
	double phibulk;
	double fraction(int); 
	string freedom;
	MoleculeType MolType; 
	double n; 
	double GN,GN1,GN2; 
	double norm;
	int chainlength;
	bool save_memory; 
	bool compute_phi_alias; 
	string composition;
	vector<int> Gnr;
	vector<int> mon_nr;
	vector<int> n_mon; 
	vector<int> molmon_nr; 
	vector<int> memory;
	double *phi;
	double *H_phi;
	double *phitot;
	double *H_phitot; 
	double *Gg_f;
	double *Gg_b; 
	double *Gs; 
	int tag_segment; 

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
	double Theta(void);
	string GetValue(string); 
	int GetMonNr(string);
	int GetAlNr(string);
	bool MakeMonList(void);  
	bool IsPinned(void);
	bool IsTagged(void); 
	bool IsCharged(void); 
	void AllocateMemory(void);
	bool PrepareForCalculations(); 
	bool ComputePhi(); 
	void propagate_forward(double*, double*, int&, int,int);
	void propagate_backward(double*, double*, double*,int&,int,int);
	bool ComputePhiMon();
	bool ComputePhiLin();

	
};

#endif
