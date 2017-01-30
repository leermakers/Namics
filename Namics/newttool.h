#ifndef NEWTTOOLxH
#define NEWTTOOLxH

#include <math.h>
#include <float.h>

double norm2(double*,int);
void decompos(double*, int, int, int&);
void decompos(double*, int, int&);
int signdeterminant(double*, int, int);
int signdeterminant(double*, int);
void multiply(double*, double, double*, double*, int, int);
void multiply(double*, double, double*, double*, int);
void updateneg(double*,double*, int, int, double);
void updateneg(double* ,double* , int, double);
void updatpos(double*, double*, double*, int, int, double);
void updatpos(double*, double*, double*, int, double);
void gausa(double*, double*, double*, int, int);
void gausa(double*, double*, double*, int);
void gausb(double*, double*, int, int);
void gausb(double*, double*, int);
		

double linecriterion(double*, double*, double*, double*, int);
double newfunction(double*, double*, int);	
double residue(double*, double*, double*, int, double);

#endif
