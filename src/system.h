#ifndef SYSTEMxH
#define SYSTEMxH
#include "namics.h"
#include "segment.h"
#include "state.h"
#include "reaction.h"
#include "molecule.h"
#include "lattice.h"

class System {
public:
	System(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<State*>,vector<Reaction*>,vector<Molecule*>,string);

~System();

	string name;
	vector<Input*> In;
	Real* CHI;
	vector<Segment*> Seg;
	vector<State*> Sta;
	vector<Reaction*> Rea;
	vector<Molecule*> Mol;
	vector<Lattice*> Lat;
	vector<int> SysMonList;
	vector<int> ItMonList;
	vector<int> ItStateList;
	vector<int> StatelessMonList;
	vector<int> SysMolMonList;
	vector<int> FrozenList;
	vector<int> SysTagList;
	vector<int> SysClampList;
	vector<int> XmolList;
	vector<int> XstateList_1;
	vector<int> XstateList_2;
	vector<int> DeltaMolList;
	vector<int> ExternsMolList;
	vector<int> px;
	vector<int> py;
	vector<int> pz;
	vector<int> pointx;
	vector<int> liney;
	vector<int> surfacez;
	vector<int> Xn_1;
	Real FreeEnergy; 
	Real GrandPotential;
	Real KJ0;
	Real Kbar;
	Real* phitot;
	int* KSAM;
	Real* eps;
	Real* H_psi;
	Real* H_q;
	Real* psi;
	Real* EE;
	Real* E;
	int* psiMask;
	bool fixedPsi0;
	bool grad_epsilon;
	bool constraintfields;
	bool externsfields;	
	Real* q;
	Real* H_GrandPotentialDensity;
	Real* H_FreeEnergyDensity;
	Real* H_alpha;
	Real* GrandPotentialDensity;
	Real* FreeEnergyDensity;
	Real* alpha;
	Real* TEMP;
	Real* H_BETA;
	Real* BETA;
	Real* H_mu_ex;
	Real* mu_ex;
	int* H_beta;
	int* beta; 
	bool GPU;
	int n_mol;
	int solvent;
	int neutralizer;
	int tag_segment;
	int volume;
	int boundaryless_volume;
	bool input_error;
	bool cuda;
	bool charged;
	bool internal_states;
	string initial_guess;
	string guess_inputfile;
	string final_guess;
	string guess_outputfile;
	string ConstraintType;
	string externstype;
	string delta_inputfile; 
	int Var_target;
	Real Var_target_value;

	bool prepared;

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
	Real* GetPointer(string,int&);
	int* GetPointerInt(string,int&);
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
	bool IsCharged();
	bool IsUnique(int,int);
	void AllocateMemory();
	bool PrepareForCalculations();
	bool generate_mask();
	bool ComputePhis();
	void DoElectrostatics(Real*,Real*);
	bool CheckResults(bool);
	Real GetFreeEnergy();
	Real GetGrandPotential();
	Real GetSpontaneousCurvature();
	bool CreateMu();
	bool PutVarInfo(string,string,Real);
	Real GetError();
};
#endif
