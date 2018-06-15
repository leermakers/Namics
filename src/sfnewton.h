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
	int numIterationsForHessian;

	bool ignore_newton_direction;
	Real delta_min; 

	int trouble;
	Real normg;
	bool e_info;
	bool d_info;
	bool g_info;
	bool h_info;
	bool x_info; 
	bool s_info;
	bool t_info;

	int i_info,iv;
	Real residual;
	Real epsilon;
	int* reverseDirection;
	int* mask;
	int IV;
	int iterations; 
	Real minimum;

	virtual void residuals(Real*,Real*); //x,g
	virtual void inneriteration(Real*,Real*,float*, Real, Real, Real, int); //x g accuracy nvar 
	bool getnewtondirection();
	int getiterations();
	bool ispseudohessian();

	string GetNewtonInfo(int&);
	void COMPUTEG(Real*,Real*,int,bool);
	void ResetX(Real*,int,bool);
	bool Message(bool,bool,int, int,Real, Real,string);
	 
	Real newdirection(float*, Real*,Real*, Real*,Real*, Real*, int, Real,bool); //there is only one of this.
	//void inneriteration(Real*,Real*,float*,Real,int);
	void direction(float*, Real*, Real*, Real*, Real*, int, Real,Real,bool);
	void newhessian(float*,Real*,Real*,Real*,Real*,int,Real,Real,bool);
	void resethessian(float*, Real*, Real*, int);
	void startderivatives(float*,Real*,Real*,int);
	//void newtrustregion(Real*,Real*,Real*,Real*,int); //there is only one.
	void newtrustregion(Real*,Real,Real&,Real&,Real,Real,int); //there is only one.
	Real linesearch(Real*,Real*,Real*,Real*,Real*,int, Real,bool);  //there is only one.
	Real zero(Real*,Real*,Real*,Real*,Real*,int,Real,bool);
	Real stepchange(Real*,Real*,Real*,Real*,Real*,Real*,int,Real&,bool);
	Real linecriterion(Real*, Real*, Real*, Real*,int);
	void numhessian(float*, Real*, Real*, int,bool);
	void findhessian(float* ,Real*,Real*,int,bool);
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
	bool iterate(Real*,int,int,Real,Real,Real,bool);
	bool iterate_Picard(Real*,int,int,Real,Real);
	bool iterate_DIIS(Real*,int,int, int, Real, Real);
	void Ax(Real*, Real*, int);
	void DIIS(Real* , Real*, Real*, Real* , Real* ,Real*, int, int , int, int );
 
};
#endif
