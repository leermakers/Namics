#ifndef SFNEWTONxH
#define SFNEWTONxH

class Newton {
public:
	Newton();

~Newton();
	bool pseudohessian;
	bool samehessian;
	bool e_info;
	bool s_info;
	bool t_info;
	bool ignore_newton_direction;

	int lineiterations,linetolerance,resetiteration;
	int iterations;
	int iterationlimit,i_info;

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
	
	bool hessian;

	Real accuracy;
	Real max_accuracy_for_hessian_scaling;
	int n_iterations_for_hessian;
	Real ALPHA;
	int nbits;
	int trouble;
	Real alphaMax,alphaMin,alphabound,minimum,normg;
	int *reverseDirection;
	Real delta_max;
	Real delta_min;
	int n_tr,n_reset,n_ignore;
	int iv;
	int start;
	Real tolerance;
	Real residual;
	Real epsilon;
	string method

	Real* x;
	Real* x0;
	Real* g;
	Real* g0;

	Real* p;
	Real* p0;
	Real* alpha;

	float* h;
	virtual void inneriteration(float*, Real*,Real*,Real, int);
	virtual void iterate(Real*,int); 
	virtual void residuals(Real*,Real*);

	Real residue(Real*, Real*, Real*, int, Real);
	void newdirection(float*, Real*, Real*, Real*, Real*, Real*, int, Real); //there is only one of this.
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
	void ResetX(Real*,int);
	void COMPUTEG(Real, Real, int) ;
	void Message(bool, bool, int, int,Real, Real, string);
	string GetNewtonInfo(int);

};
#endif
