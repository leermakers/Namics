#ifndef NEWTONxH
#define NEWTONxH
#include "namics.h"
#include "input.h"
#include "system.h"
#include "segment.h"
#include "lattice.h"
#include "molecule.h"
#include "tools.h"

class Newton {
public:
	Newton(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,vector<System*>,string);

~Newton();

	string name; 
//*************scheutjens **************

	int lineiterations,linetolerance,resetiteration;
	int iterations;
	int print_hessian_at_it;
	int linesearchlimit;
	rene smallAlpha;
	int maxNumSmallAlpha;
	int numIterationsSinceHessian;
	int smallAlphaCount; 
	int reverseDirectionRange;
	int numReverseDirection;
	rene maxFrReverseDirection;
	rene minAccuracyForHessian;
	rene minAccuracySoFar; 
	rene resetHessianCriterion;
	rene trustregion;
	rene trustfactor;
	bool newtondirection;
	bool reset_pseudohessian;
	bool pseudohessian;
	bool hessian; 
	bool samehessian;
	rene accuracy;
	rene max_accuracy_for_hessian_scaling;
	int n_iterations_for_hessian; 
	rene ALPHA;
	bool ignore_newton_direction;
	rene delta_min;
	int nbits;
	int trouble;
	int hessianwidth;  
	rene alphaMax,alphaMin,alphabound,minimum,normg;
//***************************
	int *reverseDirection; 
	rene delta_max;
	bool e_info;
	bool s_info;
	vector<Input*> In; 
	vector<System*> Sys;
	vector<Segment*> Seg;
	vector<Lattice*> Lat;  
	vector<Molecule*> Mol;
	int iterationlimit,m,i_info;
	int k_diis,it; 
	int n_tr,n_reset,n_ignore;
	int start;
	rene tolerance;
	rene residual;
	rene epsilon;
	string method;
	string StoreFileGuess;
	string ReadFileGuess;  
	string stop_criterion;
	int iv; 

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
	int GetValue(string,int&,rene&,string&);

	rene* xx;
	rene* x0;
	rene* g;
	rene* xR;
	rene* x_x0;
	rene* alpha;
	rene* Aij;
	rene* Ci;
	rene* Apij; 
	rene* p;
	rene* p0;
	rene* g0;
	float* h;
	int* mask;
		
	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string); 
	string GetValue(string); 
	bool Guess();
	bool Solve();
	void AllocateMemory(); 
	bool PrepareForCalculations(void);
	void Ax(rene* , rene* , int );
	void DIIS(rene* , rene* , rene* , rene*, rene* ,rene* , int , int , int );
	void ComputeG(rene*); 
	void COMPUTEG(rene*,rene*,int); 
	void ComputePhis();
	void ComputeG_ext();
	bool Iterate_Picard();
	bool Iterate_DIIS();
	void Message(int, int,rene, rene); 
	bool PutU();
	
//**********Scheutjens****************
	
	void iterate(rene*,int); //there is only one iterate;
	void newdirection(float*, rene*, rene*, rene*, rene*, rene*, int, rene); //there is only one of this.
	void inneriteration(float*, rene*,rene*,rene, int);  
	void direction(float*, rene*, rene*, rene*, rene*, int, rene);
	void newhessian(float*,rene*,rene*,rene*,rene*,int);
	void resethessian(float*, rene*, rene*, int);
	void startderivatives(float*,rene*,rene*,int);

	void newtrustregion(rene*,rene*,rene*,rene*,int); //there is only one. 
	rene linesearch(rene*,rene*,rene*,rene*,rene*,int, rene);  //there is only one. 
	rene zero(rene*,rene*,rene*,rene*,rene*,int,rene);
	rene stepchange(rene*,rene*,rene*,rene*,rene*,rene*,int,rene&);
	rene linecriterion(rene*, rene*, rene*, rene*,int); 
	void numhessian(float*, rene*, rene*, int);
	void findhessian(float* ,rene*,rene*,int);
	void decomposition(float*,int,int&); 
	rene norm2(rene*,int);
	void decompos(float*, int, int&);
	int signdeterminant(float*, int);
	void multiply(rene*, rene, float*, rene*, int);
	void updateneg(float* ,rene* , int, rene);
	void updatpos(float*, rene*, rene*, int, rene);
	void gausa(float*, rene*, rene*, int);
	void gausb(float*, rene*, int);	
	rene newfunction(rene*, rene*, int);	
	rene residue(rene*, rene*, rene*, int, rene);	
//************************************
};
#endif
