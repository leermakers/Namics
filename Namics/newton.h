#ifndef NEWTONxH
#define NEWTONxH
#include "namics.h"
#include "input.h"
#include "system.h"
#include "segment.h"
#include "lattice.h"
#include "molecule.h"
#include "tools.h"
#include "variate.h"

class Newton {
public:
	Newton(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,vector<System*>,vector<Variate*>,string);

~Newton();

	string name;
//*************scheutjens **************

	int lineiterations,linetolerance,resetiteration;
	int iterations;
	int print_hessian_at_it;
	int linesearchlimit;
	Real smallAlpha;
	int maxNumSmallAlpha;
	int numIterationsSinceHessian;
	int smallAlphaCount;
	int reverseDirectionRange;
	int numReverseDirection;
	Real maxFrReverseDirection;
	Real minAccuracyForHessian;
	Real minAccuracySoFar;
	Real resetHessianCriterion;
	Real trustregion;
	Real trustfactor;
	bool newtondirection;
	bool reset_pseudohessian;
	bool pseudohessian;
	bool hessian;
	bool samehessian;
	Real accuracy;
	Real max_accuracy_for_hessian_scaling;
	int n_iterations_for_hessian;
	Real ALPHA;
	bool ignore_newton_direction;
	Real delta_min;
	int nbits;
	int trouble;
	int hessianwidth;
	Real alphaMax,alphaMin,alphabound,minimum,normg;
//***************************
	int *reverseDirection;
	Real delta_max;
	bool e_info;
	bool s_info;
	bool super_e_info;
	bool super_s_info;
	vector<Input*> In;
	vector<System*> Sys;
	vector<Segment*> Seg;
	vector<Lattice*> Lat;
	vector<Molecule*> Mol;
	vector<Variate*> Var;
	int iterationlimit,m,i_info;
	int super_iterationlimit;
	int k_diis,it;
	int n_tr,n_reset,n_ignore;
	int start;
	Real tolerance;
	Real super_tolerance;
	Real residual;
	Real epsilon;
	string method;
	string StoreFileGuess;
	string ReadFileGuess;
	string stop_criterion;
	int iv;

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


	Real* xx;
	Real* x0;
	Real* g;
	Real* xR;
	Real* x_x0;
	Real* alpha;
	Real* Aij;
	Real* Ci;
	Real* Apij;
	Real* p;
	Real* p0;
	Real* g0;
	float* h;
	int* mask;
	Real* r;
	Real* d;
	Real* q;
	Real* s;

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string);
	string GetValue(string);
	void Copy(Real*,Real*,int,int,int);
	bool Guess(Real*,string,vector<string>,bool,int,int,int);
	string GetNewtonInfo(int&);
	bool Solve(bool);
	bool SolveMesodyn(Real*); //first argument should contain rho
	bool SuperIterate(int,int,int,int);
	void DeAllocateMemory();
	void AllocateMemory();
	bool PrepareForCalculations(void);
	void Ax(Real* , Real* , int );
	void DIIS(Real* , Real* , Real* , Real*, Real* ,Real* , int ,int, int , int );

	void COMPUTEG(Real*,Real*,int);
	void ResetX(Real*,int);
	void ComputePhis();
	void ComputeG(Real*);
	void ComputeG_ext();
	void ComputeG_mesodyn(Real*);
	void MinvTimesQ(Real*, Real*, int);
	bool SolvePsi(Real*, Real*, Real*);
	bool Iterate_Picard();
	bool Iterate_DIIS(Real*);
	bool Iterate_DIIS();
	void Message(bool,bool,int, int,Real, Real,string);
	bool PutU();

//**********Scheutjens****************

	void iterate(Real*,int); //there is only one iterate;
	void newdirection(float*, Real*, Real*, Real*, Real*, Real*, int, Real); //there is only one of this.
	void inneriteration(float*, Real*,Real*,Real, int);
	void direction(float*, Real*, Real*, Real*, Real*, int, Real);
	void newhessian(float*,Real*,Real*,Real*,Real*,int);
	void resethessian(float*, Real*, Real*, int);
	void startderivatives(float*,Real*,Real*,int);

	void newtrustregion(Real*,Real*,Real*,Real*,int); //there is only one.
	Real linesearch(Real*,Real*,Real*,Real*,Real*,int, Real);  //there is only one.
	Real zero(Real*,Real*,Real*,Real*,Real*,int,Real);
	Real stepchange(Real*,Real*,Real*,Real*,Real*,Real*,int,Real&);
	Real linecriterion(Real*, Real*, Real*, Real*,int);
	void numhessian(float*, Real*, Real*, int);
	void findhessian(float* ,Real*,Real*,int);
	void decomposition(float*,int,int&);
	Real norm2(Real*,int);
	void decompos(float*, int, int&);
	int signdeterminant(float*, int);
	void multiply(Real*, Real, float*, Real*, int);
	void updateneg(float* ,Real* , int, Real);
	void updatpos(float*, Real*, Real*, int, Real);
	void gausa(float*, Real*, Real*, int);
	void gausb(float*, Real*, int);
	Real newfunction(Real*, Real*, int);
	Real residue(Real*, Real*, Real*, int, Real);
//************************************
};
#endif
