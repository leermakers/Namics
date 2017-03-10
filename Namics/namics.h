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
#ifdef _WIN32
    //define lapack for Windows (32-bit and 64-bit, this part is common)
   #ifdef _WIN64
      //define lapack for Windows (64-bit only)
   #else
      //define lapack for Windows (32-bit only)
   #endif
#elif __APPLE__
    //define lapack for apple
    #if TARGET_IPHONE_SIMULATOR
         // iOS Simulator
    #elif TARGET_OS_IPHONE
        // iOS device
    #elif TARGET_OS_MAC
        // Other kinds of Mac OS
    #else
    #   error "Unknown Apple platform"
    #endif
#elif __linux
    #include <lapacke.h>
#elif __unix__ // all unices not caught above
    #include <clapack.h>

#elif defined(_POSIX_VERSION)
    // POSIX
#else
#   error "Unknown compiler"
#endif

#include <iomanip> 
using namespace std;
typedef double rene;

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

enum MoleculeType {monomer, linear, branched};


#endif
