#ifndef SOLVE_SCFxH
#define SOLVE_SCFxH
#include "namics.h"
#include "input.h"
#include "system.h"
#include "segment.h"
#include "state.h"
#include "reaction.h"
#include "lattice.h"
#include "molecule.h"
#include "tools.h"
#include "variate.h"
#include "sfnewton.h"
#include <functional>

class Solve_scf : public SFNewton {
public:
	Solve_scf() {};

	Solve_scf(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<State*>,vector<Reaction*>,vector<Molecule*>,vector<System*>,vector<Variate*>,string);

	~Solve_scf();

	enum rescue {
		NONE,
		ZERO,
		M,
		DELTA_MAX,
	};

	bool attempt_DIIS_rescue();
	rescue rescue_status;

	string name;
	vector<Input*> In;
	vector<System*> Sys;
	vector<Segment*> Seg;
	vector<Lattice*> Lat;
	vector<Molecule*> Mol;
	vector<Variate*> Var;
	vector<State*> Sta;
	vector<Reaction*> Rea;

	int start;
	Real Value_tolerance;
	string SCF_method;
	string gradients;
	string StoreFileGuess;
	string ReadFileGuess;
	string stop_criterion;
	int iv;
	int m;

	Real super_tolerance,tolerance;
	Real super_deltamax,deltamax;
	Real super_deltamin,deltamin;

	int super_iterationlimit,iterationlimit;
	bool super_e_info, value_e_info;
	bool super_s_info, value_s_info;
	int super_i_info, value_i_info;

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
	int GetValue(string,int&,Real&,string&);
	enum iteration_method {HESSIAN,PSEUDOHESSIAN,PICARD,diis,conjugate_gradient};
	enum inner_iteration_method {super,proceed};
	enum gradient_method {classical, MESODYN, Picard, custum, WEAK};
	iteration_method solver;
	gradient_method gradient;
	inner_iteration_method control;


#ifdef CUDA
	Real *x0;
	Real *g;
	Real *xR;
	Real *x_x0;
#endif
	Real *xx;
	int *SIGN;
	Real* alpha;
	bool mesodyn;
	Real* RHO;
	int value_search;
	int value_target;
	int value_ets,old_value_ets;
	int value_etm,old_value_etm;

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string);
	string GetValue(string);
	void Copy(Real*,Real*,int,int,int,int);
	bool Guess(Real*,string,vector<string>,vector<string>,bool,int,int,int,int);

	bool Solve(bool);
	bool SolveMesodyn(function< void(Real*, size_t) >, function< Real*() >); //first argument should contain rho
	function< Real*() > mesodyn_flux;
	function< void(Real*, size_t) > mesodyn_load_alpha;
	bool SuperIterate(int,int,int,int);
	void DeAllocateMemory();
	void AllocateMemory();
	bool PrepareForCalculations(void);
	void ComputePhis();
	bool PutU();
	void residuals(Real*,Real*);
	void gradient_log(Real*, int, int, int, int);
	void gradient_quotient(Real*, int, int, int, int);
	void gradient_minus(Real*, int, int, int, int);
	function<void(Real*, int, int, int, int)> target_function;

	void inneriteration(Real*,Real*,float*,Real,Real&,Real,int);

};
#endif
