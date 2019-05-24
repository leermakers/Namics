#ifndef NEWTTOOLxH
#define NEWTTOOLxH

#include <math.h>
#include <float.h>
#include "namics.h"

Real norm2(Real*,int);
void decompos(Real*, int, int, int&);
void decompos(Real*, int, int&);
int signdeterminant(Real*, int, int);
int signdeterminant(Real*, int);
void multiply(Real*, Real, Real*, Real*, int, int);
void multiply(Real*, Real, Real*, Real*, int);
void updateneg(Real*,Real*, int, int, Real);
void updateneg(Real* ,Real* , int, Real);
void updatpos(Real*, Real*, Real*, int, int, Real);
void updatpos(Real*, Real*, Real*, int, Real);
void gausa(Real*, Real*, Real*, int, int);
void gausa(Real*, Real*, Real*, int);
void gausb(Real*, Real*, int, int);
void gausb(Real*, Real*, int);
		

Real linecriterion(Real*, Real*, Real*, Real*, int);
Real newfunction(Real*, Real*, int);	
Real residue(Real*, Real*, Real*, int, Real);

#endif
