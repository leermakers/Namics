#ifndef NEWTTOOLxH
#define NEWTTOOLxH

#include <math.h>
#include <float.h>
#include "namics.h"

Real norm2(Real*,int);
void decompos(float*, int, int, int&);
void decompos(float*, int, int&);
int signdeterminant(float*, int, int);
int signdeterminant(float*, int);
void multiply(Real*, Real, float*, Real*, int, int);
void multiply(Real*, Real, float*, Real*, int);
void updateneg(float*,Real*, int, int, Real);
void updateneg(float* ,Real* , int, Real);
void updatpos(float*, Real*, Real*, int, int, Real);
void updatpos(float*, Real*, Real*, int, Real);
void gausa(float*, Real*, Real*, int, int);
void gausa(float*, Real*, Real*, int);
void gausb(float*, Real*, int, int);
void gausb(float*, Real*, int);
		

Real linecriterion(Real*, Real*, Real*, Real*, int);
Real newfunction(Real*, Real*, int);	
Real residue(Real*, Real*, Real*, int, Real);

#endif
