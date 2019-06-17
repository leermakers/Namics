#ifndef SFNEWTONxH
#define SFNEWTONxH
#include <vector>
#include <cstdlib>
#include <cfloat>
#include "namics.h"

class SFNewton {
public:
	SFNewton();

 	virtual ~SFNewton();
	bool max_g;
	int nbits;
	int lineiterations,numIterationsSinceHessian,resetiteration;
	Real linetolerance;
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

	virtual void residuals(Real*,Real*) = 0; //x,g
	virtual void inneriteration(Real*,Real*,Real*, Real, Real&, Real, int) = 0; //x g accuracy nvar
	bool getnewtondirection();
	int getiterations();
	bool ispseudohessian();

	string GetNewtonInfo(int&);
	void COMPUTEG(Real*,Real*,int,bool);
	void ResetX(Real*,int,bool);
	bool Message(bool,bool,int, int,Real, Real,string);

	Real newdirection(Real*, Real*,Real*, Real*,Real*, Real*, int, Real,bool); //there is only one of this.
	//void inneriteration(Real*,Real*,Real*,Real,int);
	void direction(Real*, Real*, Real*, Real*, Real*, int, Real,Real,bool);
	void newhessian(Real*,Real*,Real*,Real*,Real*,int,Real,Real,bool);
	void resethessian(Real*, Real*, Real*, int);
	void startderivatives(Real*,Real*,Real*,int);
	//void newtrustregion(Real*,Real*,Real*,Real*,int); //there is only one.
	void newtrustregion(Real*,Real,Real&,Real&,Real,Real,int); //there is only one.
	Real linesearch(Real*,Real*,Real*,Real*,Real*,int, Real,bool);  //there is only one.
	Real zero(Real*,Real*,Real*,Real*,Real*,int,Real,bool);
	Real stepchange(Real*,Real*,Real*,Real*,Real*,Real*,int,Real&,bool);
	Real linecriterion(Real*, Real*, Real*, Real*,int);
	void numhessian(Real*, Real*, Real*, int,bool);
	void findhessian(Real* ,Real*,Real*,int,bool);
	void decomposition(Real*,int,int&);
	Real norm2(Real*,int);
	void decompos(Real*, int, int&);
	int signdeterminant(Real*, int);
	void multiply(Real*, Real, Real*, Real*, int);
	void updateneg(Real* ,Real* , int, Real);
	void updatpos(Real*, Real*, Real*, int, Real);
	void gausa(Real*, Real*, Real*, int);
	void gausb(Real*, Real*, int);
	Real newfunction(Real*, Real*, int);
	Real residue(Real*, Real*, Real*, int, Real);
//************************************
	bool iterate(Real*,int,int,Real,Real,Real,bool);
	bool iterate_Picard(Real*,int,int,Real,Real);
	bool iterate_DIIS(Real*,int,int, int, Real, Real);
	bool iterate_RF(Real*,int,int,Real,Real,string);
	void Ax(Real*, Real*, int);
	void DIIS(Real* , Real*, Real*, Real* , Real* ,Real*, int, int , int, int);
	void conjugate_gradient(Real *, int, int ,Real ) ;
	void Hd(Real *, Real *, Real *, Real *, Real *, Real*, Real);
private:
	Real computeresidual(Real*, int);

};
#endif
