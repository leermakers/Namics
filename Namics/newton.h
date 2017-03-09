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
	double smallAlpha;
	int maxNumSmallAlpha;
	int numIterationsSinceHessian;
	int smallAlphaCount; 
	int reverseDirectionRange;
	int numReverseDirection;
	double maxFrReverseDirection;
	double minAccuracyForHessian;
	double minAccuracySoFar; 
	double resetHessianCriterion;
	double trustregion;
	double trustfactor;
	bool newtondirection;
	bool reset_pseudohessian;
	bool pseudohessian;
	bool hessian; 
	bool samehessian;
	double accuracy;
	double max_accuracy_for_hessian_scaling;
	int n_iterations_for_hessian; 
	double ALPHA;
	bool ignore_newton_direction;
	double delta_min;
	int nbits;
	int trouble;
	int hessianwidth;  
	double alphaMax,alphaMin,alphabound,minimum,normg;
//***************************
	int *reverseDirection; 
	double delta_max;
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
	double tolerance;
	double residual;
	double epsilon;
	string method;
	string StoreFileGuess;
	string ReadFileGuess;  
	string stop_criterion;
	int iv; 

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
	int GetValue(string,int&,double&,string&);

	double* xx;
	double* x0;
	double* g;
	double* xR;
	double* x_x0;
	double* alpha;
	double* Aij;
	double* Ci;
	double* Apij; 
	double* p;
	double* p0;
	double* g0;
	float* h;
	int* mask;
		
	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string); 
	string GetValue(string); 
	bool Solve();
	void AllocateMemory(); 
	bool PrepareForCalculations(void);
	void Ax(double* , double* , int );
	void DIIS(double* , double* , double* , double*, double* ,double* , int , int , int );
	void ComputeG(double*); 
	void COMPUTEG(double*,double*,int); 
	void ComputePhis();
	void ComputeG_ext();
	bool Iterate_Picard();
	bool Iterate_DIIS();
	void Message(int, int,double, double); 
	bool PutU();
	
//**********Scheutjens****************
	
	void iterate(double*,int); //there is only one iterate;
	void newdirection(float*, double*, double*, double*, double*, double*, int, double); //there is only one of this.
	void inneriteration(float*, double*,double*,double, int);  
	void direction(float*, double*, double*, double*, double*, int, double);
	void newhessian(float*,double*,double*,double*,double*,int);
	void resethessian(float*, double*, double*, int);
	void startderivatives(float*,double*,double*,int);

	void newtrustregion(double*,double*,double*,double*,int); //there is only one. 
	double linesearch(double*,double*,double*,double*,double*,int, double);  //there is only one. 
	double zero(double*,double*,double*,double*,double*,int,double);
	double stepchange(double*,double*,double*,double*,double*,double*,int,double&);
	double linecriterion(double*, double*, double*, double*,int); 
	void numhessian(float*, double*, double*, int);
	void findhessian(float* ,double*,double*,int);
	void decomposition(float*,int,int&); 
	double norm2(double*,int);
	void decompos(float*, int, int&);
	int signdeterminant(float*, int);
	void multiply(double*, double, float*, double*, int);
	void updateneg(float* ,double* , int, double);
	void updatpos(float*, double*, double*, int, double);
	void gausa(float*, double*, double*, int);
	void gausb(float*, double*, int);	
	double newfunction(double*, double*, int);	
	double residue(double*, double*, double*, int, double);	
//************************************
};
#endif
