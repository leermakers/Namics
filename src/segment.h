#ifndef SEGMENTxH
#define SEGMENTxH
#include "namics.h"
#include "input.h"
#include "lattice.h"
#include "tools.h"
class Segment {
public:
	Segment(vector<Input*>,vector<Lattice*>,string,int,int);

~Segment();

	string name;
	vector<Input*> In;
	vector<Lattice*> Lat;
	vector<int> px1;
	vector<int> px2;
	vector<int> py1;
	vector<int> py2;
	vector<int> pz1;
	vector<int> pz2;
	vector<int> bx;
	vector<int> by;
	vector<int> bz;
	vector<int> px;
	vector<int> py;
	vector<int> pz;
	int n;
	int R;
	vector<int> constraint_z;
	vector<Real> constraint_phi;
	vector<Real> constraint_beta;
	bool constraints;

	vector<string> chi_name;
	vector<Real> chi;
	int clamp_nr;
	int n_seg;
	int seg_nr;
	bool unique;
	int seg_nr_of_copy;
	int state_nr_of_copy;
	bool prepared;

	Real theta_exc;
	Real M1,M2,Fl;
	Real epsilon;
	Real valence;
	Real PSI0;
	bool fixedPsi0;
	Real phibulk;
	string freedom;
	Real guess_u;
	vector<int>state_change;
	vector<Real>state_valence;
	vector<int>state_id;
	vector<string>state_name;
	vector<int>state_nr;
	vector<Real>state_alphabulk;
	vector<Real>state_phibulk;
	vector<Real>state_theta;

	string filename;
	string s_freedom;
	bool block;
	bool all_segment;
	int ns;

	int n_pos;
	int n_box;
	int* r;
	int mx,my,mz,m;
	string scale;
	int Var_steps;
	Real Var_step;
	Real Var_end_value;
	Real Var_start_value;
	int num_of_steps;
	int Var_target;
	string Var_type;
	int chi_var_seg;
	int chi_var_state;
	Real chi_value;
	Real Amplitude;
	int labda;
	int seed;
	int start;
	int var_pos;

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
	Real Volume_particles(void);
	bool Overlap(int,int);
	void PushOutput();
	Real* GetPointer(string,int&);
	int* GetPointerInt(string,int&);
	int GetValue(string,int&,Real&,string&);
	bool SetExternalPotentials();
	Real Get_g(int ) ;
	void Put_beta(int, Real );


	int* H_P;
	int* H_MASK;
	Real* H_u;
	Real* H_phi;
	Real* H_u_ext;

	Real* H_phi_state;
	Real* H_alpha;

	int* P;
	int* MASK;
	Real* G1;
	Real* phi;
	Real* phi_state;
	Real* phi_side;
	Real* u;
	Real* u_ext;

	Real* alpha;


	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
 	 std::vector<string> VALUES;
	bool CheckInput(int);
	void PutChiKEY(string);
	string GetValue(string);
	string GetFreedom();
	bool GetClamp(string);
	bool IsFree();
	bool IsPinned();
	bool IsFrozen();
	bool IsTagged();
	bool IsClamp();
	int* GetMASK();
	Real* GetPhi();
	void DeAllocateMemory();
	void AllocateMemory();
	bool PrepareForCalculations(int*,bool);
	bool PutAdsorptionGuess(Real,int*);
	bool PutTorusPotential(int);
	bool PutMembranePotential(int);
	void UpdateValence(Real*,Real*,Real*,Real*,bool);
	bool PutVarInfo(string,string,Real);
	int PutVarScan(Real,Real,int,string);
	bool ResetInitValue() ;
	bool UpdateVarInfo(int);
	void PutValue(Real);
	Real GetValue();
	int AddState(int,Real,Real,bool);
	void SetPhiSide();
	void PutAlpha(Real*,int&);
	bool CanBeReached(int, int, int, int);
};
#endif
