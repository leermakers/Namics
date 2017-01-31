#ifndef NAMICSxH
#define NAMICSxH

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream> 
#include <cstring> 
#include <string>
#include <fstream>
#include <sstream>
#define SSTR( x ) static_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x )).str()
#include <vector>
#include <iterator>
#include <algorithm>
#include <f2c.h>
#include <clapack.h>
#include <iomanip> 
using namespace std;
#ifdef CUDA
//#include <cuda.h>
//#include <cublas_v2.h>
//#include <cuda_runtime.h>
#endif

//nvcc Open.cu -lm -lcuda -lcudart -llapack -lblas -lf2c -lcublas -arch=sm_20 -o box
//I.V. Ionova, E.A. Carter, "Error vector choice in direct inversion in the iterative subspace method, J. Compt. Chem. 17, 1836-1847, 1996.

#ifndef MAINxH //here define global variables. 
extern double* BlasResult;
extern string version;
extern double e;
extern double T;
extern double k_B;  
extern double k_BT;
extern double PIE; 
extern double eps0;
extern double eps;
extern bool debug;
//extern double factor;
#endif

enum MoleculeType {monomer, linear};


#endif
