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
#include <cuda.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#endif


//#include <f2c.h>
//#include <clapack.h>

//nvcc Open.cu -lm -lcuda -lcudart -llapack -lblas -lf2c -lcublas -arch=sm_20 -o box
//I.V. Ionova, E.A. Carter, "Error vector choice in direct inversion in the iterative subspace method, J. Compt. Chem. 17, 1836-1847, 1996.

float b_length=5e-10,e=1.60217e-19,k_BT=1.38065e-23*298.15,eps0=8.85418e-12,eps=80,factor=e*e/(eps*eps0*b_length*k_BT);
//conversion from volume fraction salt to molar salt concentration is in this setting approximately c_s=12*phib_s.
int N,N1,N2,n_seg,n_sol,n_mol,n_cosol,n_box,iterations=1000;
int N_A,N_B;
float CHI,alpha_seg,phib_s;
float Theta_A, GNA, GNB;
int i,j,k,k_diis,m=20,s,it,iv,Mx,My,Mz,M,jx,jy,MX,MY,MZ,JX,JY,MM,bx1,by1,bz1,bxm,bym,bzm,BX1,BY1,BZ1,BXM,BYM,BZM;
float error = 1, tolerance = 1e-3, eta = 0.1;//because of float, tolerance can not set to very low values...
float *Aij,*Ci,*Apij,*phi,*rho,*g1,*phitot,*G1,*alpha,*Gg_f,*Gg_b,*phi_side,*x,*x0,*g,*xR,*x_x0,*mask,*GN_A,*GN_B,*MASK,*KSAM;
float *GG_F;
int *Px,*Py,*Pz,*Bx, *By, *Bz;
float *u;
float *phib;
bool MEmulsion,Membrane;
bool charges = true;
float *H_mask, *H_psi, *H_phi, *H_u, *H_MASK, *H_KSAM, *H_g1, *H_G1, *H_rho, *H_GN_A, *H_GN_B, *H_PHI;
float *q, *psi, *psi_0, *psi_side, *phi_na, *phi_cl, *PHI;
int *H_Px,*H_Py,*H_Pz, *H_Bx, *H_By, *H_Bz;
std::vector<string> cal_types;
string cal_type; 
enum MoleculeType {monomer, linear};

