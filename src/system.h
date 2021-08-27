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
	bool all_system;
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
	vector<int> px;
	vector<int> py;
	vector<int> pz;
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
	Real phi_ratio;
	int* H_beta;
	int* beta;
	bool GPU;
	bool first_pass;
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
	bool local_solution;
	bool do_blocks;
	int split;
	string initial_guess;
	string guess_inputfile;
	string final_guess;
	string guess_outputfile;
	string ConstraintType;
	string delta_inputfile;
	int Var_target;
	Real Var_target_value;
	int start;
	int extra_constraints;
	Real old_residual;
	int progress;

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
	void DeAllocateMemory();
	bool CheckInput(int);
	void PutParameter(string);
	string GetValue(string);
	string GetMonName(int );
	int GetMonNr(string);
	bool CheckChi_values(int);
	bool MakeItsLists();
	bool IsCharged();
	bool IsUnique(int,int);
	void AllocateMemory();
	bool PrepareForCalculations(bool);
	bool generate_mask();
	bool ComputePhis(Real);
	void ComputePhis(Real*,bool,Real);
	bool Put_U(Real*);
	bool PutU(Real*);
	void Classical_residual(Real* ,Real*,Real,int,int);

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
