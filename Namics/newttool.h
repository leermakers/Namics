#ifndef NEWTTOOLxH
#define NEWTTOOLxH

#include <math.h>
#include <float.h>

void decompos(float * const h, const int nvar, int ntr);
void updatpos(float * const l, double * const w, double * const v,const int nvar, double alpha);
void gausa(const float *const l, double *const dup, const double *const g,int nvar);
void gausb(const float *const du,double *const p, const int nvar);

void multiply(double *const v,double alpha,const float *const h,
		const double *const w, int n);
double innerproduct(const double *const, const double *const,int);
double norm2(const double *const,int);
void updateneg(float * const,double * const,int,double);
int signdeterminant(const float *const,int);

#endif
