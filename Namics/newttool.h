#ifndef NEWTTOOLxH
#define NEWTTOOLxH

#include <math.h>
#include <float.h>

double norm2(double*,int);
void decompos(float*, int, int, int&);
void decompos(float*, int, int&);
int signdeterminant(float*, int, int);
int signdeterminant(float*, int);
void multiply(double*, double, float*, double*, int, int);
void multiply(double*, double, float*, double*, int);
void updateneg(float*,double*, int, int, double);
void updateneg(float* ,double* , int, double);
void updatpos(float*, double*, double*, int, int, double);
void updatpos(float*, double*, double*, int, double);
void gausa(float*, double*, double*, int, int);
void gausa(float*, double*, double*, int);
void gausb(float*, double*, int, int);
void gausb(float*, double*, int);
		

double linecriterion(double*, double*, double*, double*, int);
double newfunction(double*, double*, int);	
double residue(double*, double*, double*, int, double);

#endif
