#ifndef SFNEWTONxH
#define SFNEWTONxH
#include <vector>
#include <stdlib.h>  
#include <math.h>
#include <float.h>
#include "namics.h"

class SFNewton {
public: 
	SFNewton();

~SFNewton();
	int nbits;	
	int it, iterations;
	int lineiterations,linetolerance,numIterationsSinceHessian,resetiteration;
	int print_hessian_at_it;
	int linesearchlimit;
	Real smallAlpha;
	int maxNumSmallAlpha;
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

	int trouble;
	Real alphaMax,alphaMin,alphabound,minimum,normg;
	Real delta_max;
	bool e_info;
	bool d_info;
	bool g_info;
	bool h_info;
	bool x_info; 
	bool s_info;
	bool t_info;
	bool filter;

	int iterationlimit,i_info,iv;
	Real tolerance;
	Real residual;
	Real epsilon;
	int* reverseDirection;
	int* mask;
	int IV;

	virtual void residuals(Real*,Real*); //x,g
	virtual void inneriteration(Real*,Real*,float*, Real, int); //x g accuracy nvar 

	string GetNewtonInfo(int&);
	void COMPUTEG(Real*,Real*,int);
	void ResetX(Real*,int);
	bool Message(bool,bool,int, int,Real, Real,string);

	 
	void newdirection(float*, Real*,Real*, Real*,Real*, Real*, int, Real); //there is only one of this.
	//void inneriteration(Real*,Real*,float*,Real,int);
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
	bool iterate(Real*,int);
	bool iterate_Picard(Real*,int);
	bool iterate_DIIS(Real*,int,int);
	void Ax(Real*, Real*, int);
	void DIIS(Real* , Real*, Real*, Real* , Real* ,Real*, int, int , int, int );
 
};
#endif
