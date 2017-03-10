#ifndef TOOLSxH
#define TOOLSxH
#include "namics.h"
#ifdef CUDA
#include <cuda.h>
//#include <cublas_v2.h>
#include <cuda_runtime.h>

//extern cublasStatus_t stat;
//extern cublasHandle_t handle;
extern const int block_size;
 
__global__ void distributeg1(Real, Real, int*,  int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int );
__global__ void collectphi(Real*, Real*,Real*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
__global__ void sum(Real*, Real*, int);
__global__ void dot(Real*, Real*, Real*, int);
__global__ void composition(Real*, Real*, Real*,Real*,Real,int);
__global__ void times(Real*, Real*, Real*, int);
__global__ void times(Real*, Real*, int*, int);
__global__ void addtimes(Real*, Real*, Real*, int);
__global__ void norm(Real*, Real , int );
__global__ void zero(Real*, int);
__global__ void zero(int*, int);
__global__ void unity(Real*, int);
__global__ void cp (Real*, Real*, int);
__global__ void cp (Real*, int*, int);
__global__ void yisaplusctimesb(Real*, Real*,Real*,Real,int);
__global__ void yisaminb(Real*, Real*,Real*, int);
__global__ void yisaplusc(Real*, Real*, Real, int);
__global__ void yisaplusb(Real*, Real*,Real*, int);
__global__ void yplusisctimesx(Real*, Real*, Real, int);
__global__ void updatealpha(Real*, Real*, Real, int);
__global__ void picard(Real*, Real*, Real, int);
__global__ void add(Real*, Real*, int);
__global__ void add(int*, int*, int);
__global__ void dubble(Real*, Real*, Real, int);
__global__ void boltzmann(Real*, Real*, int);
__global__ void invert(int*, int*, int);
__global__ void invert(Real*, Real*, int);
__global__ void putalpha(Real*,Real*,Real*,Real,Real,int);
__global__ void div(Real*,Real*,int);
__global__ void oneminusphitot(Real*, Real*, int);
__global__ void addg(Real*, Real*, Real*, int);
__global__ void computegn(Real*, Real*, int, int);
void TransferDataToHost(Real*, Real*, int);
void TransferDataToDevice(Real*, Real*, int);
void TransferIntDataToHost(int*, int*, int);
void TransferIntDataToDevice(int*, int*, int);
__global__ void bx(Real*, int, int, int, int, int, int, int);
__global__ void b_x(Real*, int, int, int, int, int, int, int);
__global__ void by(Real*, int, int, int, int, int, int, int);
__global__ void b_y(Real*, int, int, int, int, int, int, int);
__global__ void bz(Real*, int, int, int, int, int, int, int);
__global__ void b_z(Real*, int, int, int, int, int, int, int);
__global__ void bx(int*, int, int, int, int, int, int, int);
__global__ void b_x(int*, int, int, int, int, int, int, int);
__global__ void by(int*, int, int, int, int, int, int, int);
__global__ void b_y(int*, int, int, int, int, int, int, int);
__global__ void bz(int*, int, int, int, int, int, int, int);
__global__ void b_z(int*, int, int, int, int, int, int, int);
void Dot(Real&, Real*,Real*,int);
void Sum(Real&,Real*,int);
bool GPU_present();
Real *AllOnDev(int);
int *AllIntOnDev(int);
void AddTimes(Real*, Real*, Real*, int);
void Times(Real*, Real*, Real*, int);
void Times(Real*, Real*, int*, int);
void Composition(Real*, Real*, Real*,Real*,Real,int);
void Norm(Real*, Real, int);
void Zero(Real*, int);
void Zero(int*, int);
void Unity(Real*,int);
void Cp(Real *, Real *, int);
void Cp(Real *, int *, int);
void YisAplusCtimesB(Real*, Real*, Real*, Real,int);
void YisAminB(Real*, Real*, Real*, int);
void YisAplusC(Real*, Real*, Real, int);
void YisAplusB(Real*, Real*, Real*, int);
void YplusisCtimesX(Real*, Real*, Real, int);
void UpdateAlpha(Real*, Real*, Real, int);
void Picard(Real*, Real*, Real, int);
void Add(Real*, Real*, int);
void Add(int*, int*, int);
void Dubble(Real*, Real*, Real,int);
void Boltzmann(Real*, Real*, int);
void Invert(int*, int*, int);
void Invert(Real*, Real*, int);
void PutAlpha(Real*, Real*, Real*, Real, Real, int);
void Div(Real*, Real*, int);
void AddG(Real*, Real*, Real*, int);
void OneMinusPhitot(Real*, Real*, int);
void ComputeGN(Real*, Real*, int, int);
void SetBoundaries(Real*,    int, int, int, int, int, int, int, int, int, int, int);
void RemoveBoundaries(Real*, int, int, int, int, int, int, int, int, int, int, int);
void SetBoundaries(int*,       int, int, int, int, int, int, int, int, int, int, int);
void RemoveBoundaries(int*,    int, int, int, int, int, int, int, int, int, int, int);
#else
void bx(Real *, int, int, int, int, int, int, int);
void b_x(Real*, int, int, int, int, int, int, int);
void by(Real *, int, int, int, int, int, int, int);
void b_y(Real*, int, int, int, int, int, int, int);
void bz(Real*, int, int, int, int, int, int, int);
void b_z(Real*, int, int, int, int, int, int, int);
void bx(int *, int, int, int, int, int, int, int);
void b_x(int*, int, int, int, int, int, int, int);
void by(int *, int, int, int, int, int, int, int);
void b_y(int*, int, int, int, int, int, int, int);
void bz(int*, int, int, int, int, int, int, int);
void b_z(int*, int, int, int, int, int, int, int);
void Dot(Real&,Real*, Real* ,int);
void Sum(Real&,Real*,int);
bool GPU_present();
Real *AllOnDev(int);
int *AllIntOnDev(int);
void AddTimes(Real*, Real*, Real*, int);
void Times(Real*, Real*, Real*, int);
void Times(Real*, Real*, int*, int);
void Composition(Real*, Real*, Real*,Real*,Real,int);
void Norm(Real*, Real, int);
void Zero(Real*, int);
void Zero(int*, int);
void Unity(Real*,int);
void Cp(Real *, Real *, int);
void Cp(Real *, int *, int);
void YisAplusCtimesB(Real*, Real*, Real*, Real,int);
void YisAminB(Real*, Real*, Real*, int);
void YisAplusC(Real*, Real*, Real, int);
void YisAplusB(Real*, Real*, Real*, int);
void YplusisCtimesX(Real*, Real*, Real, int);
void UpdateAlpha(Real*, Real*, Real, int);
void Picard(Real*, Real*, Real, int);
void Add(Real*, Real*, int);
void Add(int*, int*, int);
void Dubble(Real*, Real*, Real,int);
void Boltzmann(Real*, Real*, int);
void Invert(int*, int*, int);
void Invert(Real*, Real*, int);
void PutAlpha(Real*, Real*, Real*, Real, Real, int);
void Div(Real*, Real*, int);
void AddG(Real*, Real*, Real*, int);
void OneMinusPhitot(Real*, Real*, int);
void ComputeGN(Real*, Real*, int, int);
void SetBoundaries(Real*,    int, int, int, int, int, int, int, int, int, int, int);
void RemoveBoundaries(Real*, int, int, int, int, int, int, int, int, int, int, int);
void SetBoundaries(int*,       int, int, int, int, int, int, int, int, int, int, int);
void RemoveBoundaries(int*,    int, int, int, int, int, int, int, int, int, int, int);
#endif

void H_Zero(Real*, int);
void H_Zero(float*, int);
void H_Zero(int*, int);
Real H_Sum(Real*, int);
Real H_Dot(Real*, Real*, int);
void H_Invert(int*, int*, int);
void H_Invert(Real*, Real*, int);

#endif

Real pythag(Real, Real);
void svdcmp(Real **, int , int , Real *, Real **);





