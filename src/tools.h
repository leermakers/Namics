#ifndef TOOLSxH
#define TOOLSxH
#include "namics.h"
#include <numeric>
#include "lattice.h"

#ifdef CUDA
#include <cuda.h>
//#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <memory>

//extern cublasStatus_t stat;
//extern cublasHandle_t handle;
extern const int block_size;

__global__ void distributeg1(Real*, Real*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
__global__ void collectphi(Real*, Real*, Real*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
__global__ void sum(Real*, Real*, int);
__global__ void dot(Real*, Real*, Real*, int);
__global__ void composition(Real*, Real*, Real*, Real*, Real, int);
__global__ void times(Real*, Real*, Real*, int);
__global__ void times(Real*, Real*, int*, int);
__global__ void addtimes(Real*, Real*, Real*, int);
__global__ void norm(Real*, Real, int);
__global__ void zero(Real*, int);
__global__ void zero(int*, int);
__global__ void unity(Real*, int);
__global__ void cp(Real*, Real*, int);
__global__ void cp(Real*, int*, int);
__global__ void yisaplusctimesb(Real*, Real*, Real*, Real, int);
__global__ void yisaminb(Real*, Real*, Real*, int);
__global__ void yisaplusc(Real*, Real*, Real, int);
__global__ void yisaplusb(Real*, Real*, Real*, int);
__global__ void yplusisctimesx(Real*, Real*, Real, int);
__global__ void updatealpha(Real*, Real*, Real, int);
__global__ void picard(Real*, Real*, Real, int);
__global__ void add(Real*, Real*, int);
__global__ void add(int*, int*, int);
__global__ void dubble(Real*, Real*, Real, int);
__global__ void minlog(Real*, Real*, int);
__global__ void boltzmann(Real*, Real*, int);
__global__ void invert(int*, int*, int);
__global__ void invert(Real*, Real*, int);
__global__ void addgradsquare(Real*, Real*, Real*, Real*, int);
__global__ void putalpha(Real*, Real*, Real*, Real, Real, int);
__global__ void putalpha(Real*, Real*, Real, Real, int);
__global__ void div(Real*, Real*, int);
__global__ void oneminusphitot(Real*, Real*, int);
__global__ void addg(Real*, Real*, Real*, int);
__global__ void computegn(Real*, Real*, int, int);
__global__ void overwritec(Real*, int*, Real, int);
__global__ void overwritea(Real*, int*, Real*, int);
__global__ void upq(Real*, Real*, Real*, Real*, int, int, Real, int*, int);
__global__ void uppsi(Real*, Real*, Real*, Real*, int, int, Real, int*, int);
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
void Dot(Real&, Real*, Real*, int);
void Sum(Real&, Real*, int);
bool GPU_present();
Real* AllOnDev(int);
int* AllIntOnDev(int);
void AddTimes(Real*, Real*, Real*, int);
void Times(Real*, Real*, Real*, int);
void Times(Real*, Real*, int*, int);
void Composition(Real*, Real*, Real*, Real*, Real, int);
void Norm(Real*, Real, int);
void Zero(Real*, int);
void Zero(int*, int);
void Unity(Real*, int);
void Cp(Real*, Real*, int);
void Cp(Real*, int*, int);
void YisAplusCtimesB(Real*, Real*, Real*, Real, int);
void YisAminB(Real*, Real*, Real*, int);
void YisAplusC(Real*, Real*, Real, int);
void YisAplusB(Real*, Real*, Real*, int);
void YplusisCtimesX(Real*, Real*, Real, int);
void UpdateAlpha(Real*, Real*, Real, int);
void Picard(Real*, Real*, Real, int);
void Add(Real*, Real*, int);
void Add(int*, int*, int);
void Dubble(Real*, Real*, Real, int);
void MinLog(Real*, Real*, int);
void Boltzmann(Real*, Real*, int);
void Invert(int*, int*, int);
void Invert(Real*, Real*, int);
void AddGradSquare(Real*, Real*, Real*, Real*, int);
void PutAlpha(Real*, Real*, Real*, Real, Real, int);
void PutAlpha(Real*, Real*, Real, Real, int);
void Div(Real*, Real*, int);
void AddG(Real*, Real*, Real*, int);
void OneMinusPhitot(Real*, Real*, int);
void ComputeGN(Real*, Real*, int, int);
template <typename T>
inline void SetBoundaries(T*, int, int, int, int, int, int, int, int, int, int, int);
template <typename T>
inline void RemoveBoundaries(T*, int, int, int, int, int, int, int, int, int, int, int);
void DistributeG1(Real*, Real*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
void CollectPhi(Real*, Real*, Real*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
void OverwriteC(Real*, int*, Real, int);
void OverwriteA(Real*, int*, Real*, int);
void UpQ(Real*, Real*, Real*, Real*, int, int, Real, int*, int);
void UpPsi(Real*, Real*, Real*, Real*, int, int, Real, int*, int);

#else

template <typename T>
inline void Times(Real* P, Real* A, T* B, int M) {
  std::transform(A, A + M, B, P, std::multiplies<Real>());
}

template <typename T>
inline void Add(T* P, T* A, int M) {
  std::transform(P, P + M, A, P, std::plus<T>());
}

template <typename T>
inline void Sum(Real &result, T *x,int M)   {
  result = 0;
  for (int i=0; i<M; i++) result +=x[i];
}

template <typename T>
inline void Invert(T* KSAM, T* MASK, int M) {
  std::transform(MASK, MASK + M, KSAM, KSAM, [](Real A, Real B) { if (A==0) return 1.0; else return 0.0; });
}

template <typename T>
inline void Zero(T* P, int M) {
  std::fill(P, P + M, 0);
}

template <typename T>
inline void Cp(Real* P, T* A, int M) {
  std::copy(A, A + M, P);
}

void bx(Real*, int, int, int, int, int, int, int);
void b_x(Real*, int, int, int, int, int, int, int);
void by(Real*, int, int, int, int, int, int, int);
void b_y(Real*, int, int, int, int, int, int, int);
void bz(Real*, int, int, int, int, int, int, int);
void b_z(Real*, int, int, int, int, int, int, int);
void bx(int*, int, int, int, int, int, int, int);
void b_x(int*, int, int, int, int, int, int, int);
void by(int*, int, int, int, int, int, int, int);
void b_y(int*, int, int, int, int, int, int, int);
void bz(int*, int, int, int, int, int, int, int);
void b_z(int*, int, int, int, int, int, int, int);
void Dot(Real&, Real*, Real*, int);
bool GPU_present();
Real* AllOnDev(int);
int* AllIntOnDev(int);
void AddTimes(Real*, Real*, Real*, int);
void Composition(Real*, Real*, Real*, Real*, Real, int);
void Norm(Real*, Real, int);
void Unity(Real*, int);
void YisAplusCtimesB(Real*, Real*, Real*, Real, int);
void YisAminB(Real*, Real*, Real*, int);
void YisAplusC(Real*, Real*, Real, int);
void YisAplusB(Real*, Real*, Real*, int);
void YplusisCtimesX(Real*, Real*, Real, int);
void UpdateAlpha(Real*, Real*, Real, int);
void Picard(Real*, Real*, Real, int);
void Dubble(Real*, Real*, Real, int);
void MinLog(Real*, Real*, int);
void Boltzmann(Real*, Real*, int);
void AddGradSquare(Real*, Real*, Real*, Real*, int);
void PutAlpha(Real*, Real*, Real*, Real, Real, int);
void PutAlpha(Real*, Real*, Real, Real, int);
void Div(Real*, Real*, int);
void AddG(Real*, Real*, Real*, int);
void OneMinusPhitot(Real*, Real*, int);
void ComputeGN(Real*, Real*, int, int);
void DisG1(Real*, Real*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
void ColPhi(Real*, Real*, Real*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
void OverwriteC(Real*, int*, Real, int);
void OverwriteA(Real*, int*, Real*, int);
void UpQ(Real*, Real*, Real*, Real*, int, int, Real, int*, int);
void UpPsi(Real*, Real*, Real*, Real*, int, int, Real, int*, int);

template <typename T>
inline void SetBoundaries(T* P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
  bx(P, Mx + 1, My + 2, Mz + 2, bx1, bxm, jx, jy);
  by(P, Mx + 2, My + 1, Mz + 2, by1, bym, jx, jy);
  bz(P, Mx + 2, My + 2, Mz + 1, bz1, bzm, jx, jy);
}

template <typename T>
inline void RemoveBoundaries(T* P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
  b_x(P, Mx + 1, My + 2, Mz + 2, bx1, bxm, jx, jy);
  b_y(P, Mx + 2, My + 1, Mz + 2, by1, bym, jx, jy);
  b_z(P, Mx + 2, My + 2, Mz + 1, bz1, bzm, jx, jy);
}
#endif

template <typename T>
inline void H_Zero(T* H_P, int M) {
  std::fill(H_P, H_P + M, 0);
}

template <typename T>
inline Real H_Sum(T* H, int M) {
  Real Sum=0;
	for (int i=0; i<M; i++) Sum+=H[i];
  return Sum;
}

template <typename T>
inline Real H_Dot(T* A, T* B, int M) {
  Real Sum=0;
	for (int i=0; i<M; i++) Sum+=A[i]*B[i];
  return Sum;
}

template <typename T>
inline void H_Invert(T* KSAM, T* MASK, int M) {
  std::transform(MASK, MASK + M, KSAM, KSAM, [](Real A, Real B) { if (A==0) return 1.0; else return 0.0; });
}

Real pythag(Real, Real);
int svdcmp(Real**, int, int, Real*, Real**);

#endif
