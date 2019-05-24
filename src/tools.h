#ifndef TOOLSxH
#define TOOLSxH
#include "namics.h"
#include <numeric>
#include "lattice.h"

#ifdef PAR_MESODYN
	#include <thrust/extrema.h>
	#include <thrust/device_vector.h>
  #include <thrust/device_ptr.h>
  #include <thrust/pair.h>
#endif

#ifdef CUDA
#include <cuda.h>
//#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <memory>
#include <float.h>

//extern cublasStatus_t stat;
//extern cublasHandle_t handle;
extern const int block_size;

__global__ void distributeg1(Real*, Real*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
__global__ void collectphi(Real*, Real*, Real*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
__global__ void sum(Real*, Real*, int);
__global__ void sum(int*, int*, int);
__global__ void dot(Real*, Real*, Real*, int);
__global__ void max(Real *, Real *, int );
__global__ void composition(Real*, Real*, Real*, Real*, Real, int);
__global__ void times(Real*, Real*, Real*, int);
__global__ void times(Real*, Real*, int*, int);
__global__ void addtimes(Real*, Real*, Real*, int);
__global__ void norm(Real*, Real, int);
__global__ void zero(Real*, int);
__global__ void zero(int*, int);
__global__ void unity(Real*, int);
__global__ void flux_min(Real *, Real*, int, int);
__global__ void flux(Real*, Real*, Real*, int, int, int);
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
__global__ void propagate(Real *gs, Real *g_1, int JX, int JY, int JZ, int M);
__global__ void oneminusphitot(Real*, Real*, int);
__global__ void addg(Real*, Real*, Real*, int);
__global__ void computegn(Real*, Real*, int, int);
__global__ void overwritec(Real*, int*, Real, int);
__global__ void overwritea(Real*, int*, Real*, int);
__global__ void upq(Real*, Real*, Real*, Real*, int, int, Real, int*, int);
__global__ void uppsi(Real*, Real*, Real*, Real*, int, int, Real, int*, int);
template <typename T>
void TransferDataToHost(T*, T*, int);
template <typename T>
void TransferDataToDevice(T*, T*, int);
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
void Sum(int&, int*, int);
bool GPU_present(int);
int* AllIntOnDev(int);
Real* AllOnDev(int);
int* AllManagedIntOnDev(int);
Real* AllManagedOnDev(int);
Real* AllOnDev(int);
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
void Propagate(Real *gs, Real *g_1, int JX, int JY, int JZ, int M);
void AddG(Real*, Real*, Real*, int);
void OneMinusPhitot(Real*, Real*, int);
void ComputeGN(Real*, Real*, int, int);
template <typename T>
void SetBoundaries(T*, int, int, int, int, int, int, int, int, int, int, int);
template <typename T>
void RemoveBoundaries(T*, int, int, int, int, int, int, int, int, int, int, int);
namespace tools {
void DistributeG1(Real*, Real*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
void CollectPhi(Real*, Real*, Real*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int, int, int, int);
}
void OverwriteC(Real*, int*, Real, int);
void OverwriteA(Real*, int*, Real*, int);
void UpQ(Real*, Real*, Real*, Real*, int, int, Real, int*, int);
void UpPsi(Real*, Real*, Real*, Real*, int, int, Real, int*, int);

Real ComputeResidual(Real*, int);

struct saxpy_functor
{
    const double a;

    saxpy_functor(double _a) : a(_a) {}

    __host__ __device__
        double operator()(const double& x, const double& y) const { 
            return a * x + y;
        }
};

struct const_multiply_functor
{
    const double a;

    const_multiply_functor(double _a) : a(_a) {}

    __host__ __device__
        double operator()(const double& x, const double& y) const { 
            return a * x * y;
        }
};

struct order_param_functor
{

    order_param_functor() {}

    __host__ __device__
        double operator()(const double& x, const double& y) const { 
            return pow(x-y,2);
        }
};

struct is_negative_functor
{

  const double tolerance{0};

  is_negative_functor(double _tolerance=0) : tolerance(_tolerance) {}

  __host__ __device__
  bool operator()(const double &x) const
  {
    return x < 0-tolerance || x > 1+tolerance;
  }
};

struct is_not_unity_functor
{
  const double tolerance{0};

  is_not_unity_functor(double _tolerance=0) : tolerance(_tolerance) {}

  __host__ __device__
  bool operator()(const double &x) const
  {
    bool result{0};

    if (x > (1+tolerance) || x < (1-tolerance))
      result = 1;

    return result;
  }
};

#else

#include "tools_host.h"

#endif

template <typename T>
inline void H_Cp(Real* P, T* A, int M) {
  std::copy(A, A + M, P);
}

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
	for (int i = 0; i<M; i++)
    Sum += A[i]*B[i];
  return Sum;
}

template <typename T>
inline void H_Invert(T* KSAM, T* MASK, int M) {
  std::transform(MASK, MASK + M, KSAM, KSAM, [](Real A, Real B) {if (A==0) return 1.0; else return 0.0;});
}

template<typename T>
void H_PutAlpha(T *g, T *phitot, T *phi_side, T chi, T phibulk, int M)   {
	for (int i=0; i<M; i++) if (phitot[i]>0) g[i] = g[i] - chi*(phi_side[i]/phitot[i]-phibulk);
}

Real pythag(Real, Real);
int svdcmp(Real**, int, int, Real*, Real**);

#endif