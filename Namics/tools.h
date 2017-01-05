#ifndef TOOLSxH
#define TOOLSxH
#ifdef CUDA

#include <cuda.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>


extern cublasStatus_t stat;
extern cublasHandle_t handle;

__global__ void distributeg1(double, double, int*,  int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int );
__global__ void collectphi(double*, double*,double*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
__global__ void times(double*, double*, double*, int);
__global__ void addtimes(double*, double*, double*, int);
__global__ void norm(double*, double , int );
__global__ void zero(double*, int);
__global__ void cp (double*, double*, int);
__global__ void yisaminb(double*, double*,double*, int);
__global__ void yisaplusc(double*, double*, double, int);
__global__ void yisaplusb(double*, double*,double*, int);
__global__ void yplusisctimesx(double*, double*, double, int);
__global__ void updatealpha(double*, double*, double, int);
__global__ void picard(double*, double*, double, int);
__global__ void add(double*, double*, int);
__global__ void dubble(double*, double*, double, int);
__global__ void boltzmann(double*, double*, int);
__global__ void invert(double*, double*, int);
__global__ void putalpha(double*,double*,double*,double,double,int);
__global__ void div(double*,double*,int);
__global__ void oneminusphitot(double*, double*, int);
__global__ void addg(double*, double*, double*, int);
__global__ void computegn(double*, double*, int, int);
extern void TransferDataToHost(double*, double*, int);
extern void TransferDataToDevice(double*, double*, int);
extern void TransferIntDataToHost(int*, int*, int);
extern void TransferIntDataToDevice(int*, int*, int);
__global__ void bx(double*, int, int, int, int, int, int, int);
__global__ void b_x(double*, int, int, int, int, int, int, int);
__global__ void by(double*, int, int, int, int, int, int, int);
__global__ void b_y(double*, int, int, int, int, int, int, int);
__global__ void bz(double*, int, int, int, int, int, int, int);
__global__ void b_z(double*, int, int, int, int, int, int, int);
extern double Dot(double*,double*,int);
extern double Sum(double*,int);
extern bool GPU_present();
extern double *AllOnDev(int);
extern int *AllIntOnDev(int);
extern void AddTimes(double*, double*, double*, int);
extern void Times(double*, double*, double*, int);
extern void Norm(double*, double, int);
extern void Zero(double*, int);
extern void Cp(double *, double *, int);
extern void YisAminB(double*, double*, double*, int);
extern void YisAplusC(double*, double*, double, int);
extern void YisAplusB(double*, double*, double*, int);
extern void YplusisCtimesX(double*, double*, double, int);
extern void UpdateAlpha(double*, double*, double, int);
extern void Picard(double*, double*, double, int);
extern void Add(double*, double*, int);
extern void Dubble(double*, double*, double,int);
extern void Boltzmann(double*, double*, int);
extern void Invert(double*, double*, int);
extern void PutAlpha(double*, double*, double*, double, double, int);
extern void Div(double*, double*, int);
extern void AddG(double*, double*, double*, int);
extern void OneMinusPhitot(double*, double*, int);
extern void ComputeGN(double*, double*, int, int);
extern void SetBoundaries(double*, int, int, int, int, int, int, int, int, int, int, int);
extern void RemoveBoundaries(double*, int, int, int, int, int, int, int, int, int, int, int);

#else
void bx(double *, int, int, int, int, int, int, int);
void b_x(double*, int, int, int, int, int, int, int);
void by(double *, int, int, int, int, int, int, int);
void b_y(double*, int, int, int, int, int, int, int);
void bz(double*, int, int, int, int, int, int, int);
void b_z(double*, int, int, int, int, int, int, int);
double Dot(double*, double* ,int);
double Sum(double*,int);
bool GPU_present();
double *AllOnDev(int);
int *AllIntOnDev(int);
void AddTimes(double*, double*, double*, int);
void Times(double*, double*, double*, int);
void Norm(double*, double, int);
void Zero(double*, int);
void Cp(double *, double *, int);
void YisAminB(double*, double*, double*, int);
void YisAplusC(double*, double*, double, int);
void YisAplusB(double*, double*, double*, int);
void YplusisCtimesX(double*, double*, double, int);
void UpdateAlpha(double*, double*, double, int);
void Picard(double*, double*, double, int);
void Add(double*, double*, int);
void Dubble(double*, double*, double,int);
void Boltzmann(double*, double*, int);
void Invert(double*, double*, int);
void PutAlpha(double*, double*, double*, double, double, int);
void Div(double*, double*, int);
void AddG(double*, double*, double*, int);
void OneMinusPhitot(double*, double*, int);
void ComputeGN(double*, double*, int, int);
void SetBoundaries(double*, int, int, int, int, int, int, int, int, int, int, int);
void RemoveBoundaries(double*, int, int, int, int, int, int, int, int, int, int, int);
 
#endif

void H_Zero(double*, int);
double H_Sum(double*, int);
void H_Invert(double*, double*, int);

#endif







