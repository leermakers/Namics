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
virtual ~Molecule();

	string name;
	bool all_molecule;
	vector<Input*> In;
	vector<Segment*> Seg;
	vector<Lattice*> Lat;
	vector<Alias*> Al;
	vector<int> MolMonList;
	vector<int> FillRangesList;
	vector<int> MolAlList;
	int start;
	int n_mol;
	int mol_nr;
	int n_box;
	int var_al_nr;
	Real KStiff;
	Real Kw;
	Real Mu;
	Real theta;
	Real theta_range,n_range;
	Real thata_in_range;
	Real phibulk;
	Real phib1;
	string freedom;
	MoleculeType MolType;
	Real n;
	Real GN,GN1,GN2;
	Real norm;
	Real phi1,phiM,width,Dphi,pos_interface,phi_av;
	int chainlength,N;
	bool save_memory;
	bool compute_phi_alias;
	bool sym_dend;
	string composition;
	vector<int> Gnr; //generation-number
	vector<int> first_s;
	vector<int> last_s;
	vector<int> first_b;
	vector<int> last_b;
	vector<int> first_a;
	vector<int> last_a;
	vector<int> n_arm;
	vector<int> mon_nr;
	vector<int> n_mon;
	vector<int> d_mon; //degeneracy of mon : for dendrimer case.
	vector<int> molmon_nr;
	vector<int> memory;
	vector<int> last_stored;
	vector<Real> mu_state;
	vector<Real> block;
	int *Bx;
	int *By;
	int *Bz;
	int *Px1;
	int *Py1;
	int *Pz1;
	int *Px2;
	int *Py2;
	int *Pz2;
	int *H_Bx;
	int *H_By;
	int *H_Bz;
	int *H_Px1;
	int *H_Py1;
	int *H_Pz1;
	int *H_Px2;
	int *H_Py2;
	int *H_Pz2;
	Real k_stiff;
	Real *phi;
	//Real *G1;
	//Real *u;
	Real *rho;
	Real *H_phi;
	Real *H_mask1;
	Real *H_mask2;
	//Real* H_u;
	Real *H_gn;
	Real *mask1;
	Real *mask2;
	int *R_mask;
	Real *gn;
	Real *g1;
	Real *phitot;
	Real *H_phitot;
	Real *Gg_f;
	Real *Gg_b;
	Real *Gs;
	Real *UNITY;
	Real *P;
	int size;
    int Markov;
	int tag_segment;
	int Var_steps;
	Real Var_step;
	Real Var_end_value;
	Real Var_start_value;
	Real Var_start_search_value;
	int num_of_steps;
	int Var_target;
	int Var_scan_value;
	int Var_search_value;
	string Var_type;
	string scale;
	Real Var_target_value;
	int n_generations;
	bool ring;
	Real phi_LB_X;
	Real phi_UB_X;
	Real J; 
	Real Delta_MU;

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
	bool DeleteAl();
	Real* GetPointer(string,int&);
	int* GetPointerInt(string,int&);
	int GetValue(string,int&,Real&,string&);
	bool PutVarInfo(string,string,Real);
	int PutVarScan(Real,Real,int,string);

	bool ResetInitValue();
	bool UpdateVarInfo(int);
	Real GetError();
	Real GetValue();
	void PutValue(Real);
	bool ComputeWidth();
	void SetThetaBlocks(int);
	void NormPerBlock(int);
	bool Filling;

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int,bool);
	void PutParameter(string);
	bool ExpandAlias(vector<string>,string&);
	bool ExpandBrackets(string&);
	bool Interpret(string,int);
	bool GenerateTree(string,int,int&,vector<int>,vector<int>);
	//bool GenerateDend(string,int);
	bool Decomposition(string);
	int GetChainlength(void);
	Real Theta(void);
	string GetValue(string);
	int GetMonNr(string);
	int GetAlNr(string);
	bool MakeMonList(void);
	bool IsClamped(void);
	bool IsPinned(void);
	int GetPinnedSeg(void);
	bool IsTagged(void);
	bool IsCharged(void);
	Real Charge(void);
	void DeAllocateMemory(void);
	void AllocateMemory(void);
	bool PrepareForCalculations(int*);
	bool ComputePhi(Real*,int);
	virtual bool ComputePhi();
	virtual Real GetPhib1();
	virtual void AddToGP(Real*);
	virtual void AddToF(Real*);


	virtual Real fraction(int);

	Real* propagate_forward(Real*,int&,int,int,int);
	void propagate_backward(Real*,int&,int,int,int);
	Real* propagate_forward(Real*,int&,int,Real*,int,int);
	void propagate_backward(Real*,int&,int,Real*,int&,int);

};

#endif
