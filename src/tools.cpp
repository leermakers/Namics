#include "tools.h"
#include "namics.h"
#define MAX_BLOCK_SZ 512


#ifdef CUDA
//cublasHandle_t handle;
//cublasStatus_t stat=cublasCreate(&handle);
const int block_size = 256;


__device__ void atomicAdd(Real* address, Real val)
{
    unsigned long long int* address_as_ull =
                             (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
	old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
    //return __longlong_as_double(old);
}

__global__ void distributeg1(Real *G1, Real *g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	int j = blockIdx.y*blockDim.y+threadIdx.y;
	int k = blockIdx.z*blockDim.z+threadIdx.z;
	if (k < Mz && j < My && i < Mx){
		int pM=jx+jy+1+k+jx*i+jy*j; 
		int pos_r=JX+JY+1;
		int ii;
		int MXm1 = MX-1;
		int MYm1 = MY-1;
		int MZm1 = MZ-1;
		for (int p = 0;p<n_box;++p){
			if (Bx[p]+i > MXm1) ii=(Bx[p]+i-MX)*JX; else ii=(Bx[p]+i)*JX;
			if (By[p]+j > MYm1) ii+=(By[p]+j-MY)*JY; else ii+=(By[p]+j)*JY;
			if (Bz[p]+k > MZm1) ii+=(Bz[p]+k-MZ); else ii+=(Bz[p]+k);
			g1[pM]=G1[pos_r+ii];
			pM+=M;
		}
	}	
}


__global__ void collectphi(Real* phi, Real* GN, Real* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	int j = blockIdx.y*blockDim.y+threadIdx.y;
	int k = blockIdx.z*blockDim.z+threadIdx.z;
	int pos_r=MM+(JX+JY+1);
	int pM=jx+jy+1+jx*i+jy*j+k;
	int ii,jj,kk;
	int MXm1 = MX-1;
	int MYm1 = MY-1;
	int MZm1 = MZ-1;
	if (k < Mz && j < My && i < Mx){
		for (int p = 0;p<n_box;++p){
			if (Bx[p]+i > MXm1)  ii=(Bx[p]+i-MX)*JX; else ii=(Bx[p]+i)*JX; 
			if (By[p]+j > MYm1)  jj=(By[p]+j-MY)*JY; else jj=(By[p]+j)*JY;
			if (Bz[p]+k > MZm1)  kk=(Bz[p]+k-MZ); else kk=(Bz[p]+k);
			//__syncthreads(); //will not work when two boxes are idential....
			//phi[pos_r+ii+jj+kk]+=rho[pM+jx*i+jy*j+k]*Inv_GNp;
			atomicAdd(&phi[pos_r+ii+jj+kk], rho[pM]/GN[p]);
			pM+=M;
		}		
	}
}

__global__ void dot(Real *X, Real *Y, Real *Z, int M)   {
   __shared__ Real tmp[MAX_BLOCK_SZ];
   int idx = blockDim.x * blockIdx.x + threadIdx.x;
   int l_idx = threadIdx.x;

   if (idx < M) tmp[l_idx] = X[idx]*Y[idx];
   __syncthreads();

   for (int s = blockDim.x/2; s > 0; s /= 2) {
      if (l_idx < s) tmp[l_idx] += tmp[l_idx + s];
      __syncthreads();
   }
   if (threadIdx.x == 0) Z[threadIdx.x] = tmp[0];
}

__global__ void sum(Real *X, Real *Z, int M)   {
   __shared__ Real tmp[block_size];
   int idx = blockDim.x * blockIdx.x + threadIdx.x;
   int l_idx = threadIdx.x;

   if (idx < M) tmp[l_idx] = X[idx];
   __syncthreads();

   for (int s = blockDim.x/2; s > 0; s /= 2) {
      if (l_idx < s) tmp[l_idx] += tmp[l_idx + s];
      __syncthreads();
   }
   if (threadIdx.x == 0) Z[threadIdx.x] = tmp[0];
}


__global__ void composition(Real *phi, Real *Gf, Real *Gb, Real* G1, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (G1[idx]>0) phi[idx] += C*Gf[idx]*Gb[idx]/G1[idx];
}
__global__ void times(Real *P, Real *A, Real *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void times(Real *P, Real *A, int *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void addtimes(Real *P, Real *A, Real *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx]*B[idx];
}
__global__ void norm(Real *P, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] *= C;
}
__global__ void zero(Real *P, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0.0;
}
__global__ void unity(Real *P, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 1.0;
}
__global__ void zero(int *P, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0;
}
__global__ void cp (Real *P, Real *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = A[idx];
}
__global__ void cp (Real *P, int *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 1.0*A[idx];
}

__global__ void yisaplusctimesb(Real *Y, Real *A,Real *B, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]+ C * B[idx];
}
__global__ void yisaminb(Real *Y, Real *A,Real *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]-B[idx];
}
__global__ void yisaplusc(Real *Y, Real *A, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]+C;
}
__global__ void yisaplusb(Real *Y, Real *A,Real *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]+B[idx];
}
__global__ void yplusisctimesx(Real *Y, Real *X, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*X[idx];
}
__global__ void updatealpha(Real *Y, Real *X, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*(X[idx]-1.0);
}
__global__ void picard(Real *Y, Real *X, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = C*Y[idx]+(1-C)*X[idx];
}
__global__ void add(Real *P, Real *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx];
}
__global__ void add(int *P, int *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx];
}
__global__ void dubble(Real *P, Real *A, Real norm, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]*=norm/A[idx];
}
__global__ void boltzmann(Real *P, Real *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=exp(-A[idx]);
}
__global__ void overwritec(Real* P, int* Mask, Real X,int M) {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (Mask[idx]==1) P[idx] = X ; else P[idx]=0;
}
__global__ void overwritea(Real* P, int* Mask, Real* A,int M) {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (Mask[idx]==1) P[idx] = A[idx] ; else P[idx]=0;
}

__global__ void upq(Real* g, Real* q, Real* psi, Real* eps, int jx, int jy, Real C, int* Mask, int M) {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	Real* Px=psi+jx;
	Real* P_x=psi-jx;
	Real* Py=psi+jy;
	Real* P_y=psi-jy;
	Real* Pz=psi+1;
	Real* P_z=psi-1;
	Real* ex=eps+jx;
	Real* e_x=eps-jx;
	Real* ey=eps+jy;
	Real* e_y=eps-jy;
	Real* ez=eps+1;
	Real* e_z=eps-1;
	Real* e=eps;
	if (idx<M && Mask[idx]==1)  {
		q[idx]=(((e_x[idx]+e[idx])*P_x[idx] + (ex[idx]+e[idx])*Px[idx] +
			(e_y[idx]+e[idx])*P_y[idx] + (ey[idx]+e[idx])*Py[idx] +
			(e_z[idx]+e[idx])*P_z[idx] + (ez[idx]+e[idx])*Pz[idx]) -
			(e_x[idx]+ex[idx]+e_y[idx]+ey[idx]+e_z[idx]+ez[idx]+6*e[idx])*psi[idx])/C;
		g[idx] -=q[idx];
	}
}
__global__ void uppsi(Real* q, Real* psi, Real* X, Real* eps, int jx, int jy, Real C, int* Mask, int M) {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	Real* Px=X+jx;
	Real* P_x=X-jx;
	Real* Py=X+jy;
	Real* P_y=X-jy;
	Real* Pz=X+1;
	Real* P_z=X-1;
	Real* ex=eps+jx;
	Real* e_x=eps-jx;
	Real* ey=eps+jy;
	Real* e_y=eps-jy;
	Real* ez=eps+1;
	Real* e_z=eps-1;
	Real* e=eps;
	if (idx<M && Mask[idx]==0)  {
		psi[idx]=((e_x[idx]+e[idx])*P_x[idx] + (ex[idx]+e[idx])*Px[idx] +
			(e_y[idx]+e[idx])*P_y[idx] + (ey[idx]+e[idx])*Py[idx] +
			(e_z[idx]+e[idx])*P_z[idx] + (ez[idx]+e[idx])*Pz[idx] +
			C*q[idx])/(e_x[idx]+ex[idx]+e_y[idx]+ey[idx]+e_z[idx]+ez[idx]+6*e[idx]);
		q[idx] -=psi[idx];
	}
}

__global__ void invert(int *SKAM, int *MASK, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) SKAM[idx]=(MASK[idx]-1)*(MASK[idx]-1);
}
__global__ void invert(Real *SKAM, Real *MASK, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) SKAM[idx]=(MASK[idx]-1)*(MASK[idx]-1);
}
__global__ void addgradsquare(Real *EE, Real* X,  Real* Y, Real* Z, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) EE[idx]+=pow(X[idx]-Y[idx],2)+pow(Y[idx]-Z[idx],2);
}
__global__ void putalpha(Real *g,Real *phitot,Real *phi_side,Real chi,Real phibulk,int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (phitot[idx]>0) g[idx] = g[idx] - chi*(phi_side[idx]/phitot[idx]-phibulk);
}

__global__ void putalpha(Real *g,Real *phi_side,Real chi,Real phibulk,int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) g[idx] = g[idx] - chi*(phi_side[idx]-phibulk);
}

__global__ void div(Real *P,Real *A,int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (A[idx]>0) P[idx] /=A[idx];
}
__global__ void oneminusphitot(Real *g, Real *phitot, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) g[idx]= 1/phitot[idx]-1;
}
__global__ void addg(Real *g, Real *phitot, Real *alpha, int M)    {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {
		if (phitot[idx]>0) g[idx]= g[idx] -alpha[idx] +1/phitot[idx]-1; else g[idx]=0;
	}
}

//__global__ void computegn(Real *GN, Real* G, int M, int n_box)   {
//	int p= blockIdx.x*blockDim.x+threadIdx.x;
//	if (p<n_box){
//		cublasDasum(handle,M,G+p*M,1,GN+p); //in case of float we need the cublasSasum etcetera.
//	}
//}

void TransferDataToHost(Real *H, Real *D, int M)    {
	cudaMemcpy(H, D, sizeof(Real)*M,cudaMemcpyDeviceToHost);
}
void TransferDataToDevice(Real *H, Real *D, int M)    {
	cudaMemcpy(D, H, sizeof(Real)*M,cudaMemcpyHostToDevice);
}
void TransferIntDataToHost(int *H, int *D, int M)   {
	cudaMemcpy(H, D, sizeof(int)*M,cudaMemcpyDeviceToHost);
}
void TransferIntDataToDevice(int *H, int *D, int M)     {
	cudaMemcpy(D, H, sizeof(int)*M,cudaMemcpyHostToDevice);
}
#endif

#ifdef CUDA
__global__ void bx(Real *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx, jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=P[bx1_jx+idx];
		P[jx_mmx+idx]=P[jx_bxm+idx];
	}
}
__global__ void b_x(Real *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=0;
		P[jx_mmx+idx]=0;
	}
}
__global__ void by(Real *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=P[jy_by1+idx];
		P[jy_mmy+idx]=P[jy_bym+idx];
	}
}
__global__ void b_y(Real *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=0;
		P[jy_mmy+idx]=0;
	}
}
__global__ void bz(Real *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=P[idx+bz1];
		P[idx+mmz]=P[idx+bzm];
	}
}
__global__ void b_z(Real *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=0;
		P[idx+mmz]=0;
	}
}
__global__ void bx(int *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx, jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=P[bx1_jx+idx];
		P[jx_mmx+idx]=P[jx_bxm+idx];
	}
}
__global__ void b_x(int *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=0;
		P[jx_mmx+idx]=0;
	}
}
__global__ void by(int *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=P[jy_by1+idx];
		P[jy_mmy+idx]=P[jy_bym+idx];
	}
}
__global__ void b_y(int *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=0;
		P[jy_mmy+idx]=0;
	}
}
__global__ void bz(int *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=P[idx+bz1];
		P[idx+mmz]=P[idx+bzm];
	}
}
__global__ void b_z(int *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=0;
		P[idx+mmz]=0;
	}
}
#else
void bx(Real *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i;
	int jx_mmx=jx*mmx;
	int jx_bxm=jx*bxm;
	int bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=P[bx1_jx+i];
		P[jx_mmx+i]=P[jx_bxm+i];
	}
}
void b_x(Real *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}
void by(Real *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
	}
}
void b_y(Real *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}
void bz(Real *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
	}
}
void b_z(Real *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=0;
		P[i+mmz]=0;
	}
}
void bx(int *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i;
	int jx_mmx=jx*mmx;
	int jx_bxm=jx*bxm;
	int bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=P[bx1_jx+i];
		P[jx_mmx+i]=P[jx_bxm+i];
	}
}
void b_x(int *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}
void by(int *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
	}
}
void b_y(int *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}
void bz(int *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
	}
}
void b_z(int *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=0;
		P[i+mmz]=0;
	}
}
#endif

#ifdef CUDA
bool GPU_present()    {
	cudaDeviceReset();
	int deviceCount =0; cuDeviceGetCount(&deviceCount);
   	if (deviceCount ==0) printf("There is no device supporting Cuda.\n"); else cudaDeviceReset();
	//if (deviceCount>0) {
	//	stat = cublasCreate(&handle);
	//	if (stat !=CUBLAS_STATUS_SUCCESS) {printf("CUBLAS handle creation failed \n");}
	//	cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
	//}
	return deviceCount > 0;
}
#else

bool GPU_present()    {
	return false;
}
#endif


#ifdef CUDA
Real *AllOnDev(int N)    {
	Real *X;
	if (cudaSuccess != cudaMalloc((void **) &X, sizeof(Real)*N))
	printf("Memory allocation on GPU failed.\n Please reduce size of lattice and/or chain length(s) or turn on CUDA\n");
	return X;
}
int *AllIntOnDev(int N)    {
	int *X;
	if (cudaSuccess != cudaMalloc((void **) &X, sizeof(int)*N))
	printf("Memory allocation on GPU failed.\n Please reduce size of lattice and/or chain length(s) or turn on CUDA\n");
	return X;
}
#else
Real *AllOnDev(int N)    {
	Real *X=NULL;
	printf("Turn on CUDA\n");
	return X;
}
int *AllIntOnDev(int N)    {
	int *X=NULL;
	printf("Turn on CUDA\n");
	return X;
}
#endif

#ifdef CUDA
void Dot(Real &result, Real *x,Real *y, int M)   {
/**/
	Real *H_XXX=(Real*) malloc(M*sizeof(Real));
	Real *H_YYY=(Real*) malloc(M*sizeof(Real));
	TransferDataToHost(H_XXX, x, M);
	TransferDataToHost(H_YYY, y, M);
	result=H_Dot(H_XXX,H_YYY,M);
	free(H_XXX);
	free(H_YYY);
/**/


/*	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	Real *D_XX=AllOnDev(n_blocks);
	Zero(D_XX,n_blocks);
	Real *H_XX=(Real*) malloc(n_blocks*sizeof(Real));
	dot<<<n_blocks,block_size>>>(x,y,D_XX,M);
	cudaThreadSynchronize();
	TransferDataToHost(H_XX, D_XX, n_blocks);
	result=0; for (int i=0; i<n_blocks; i++) result+=H_XX[i];
	free(H_XX);
	cudaFree(D_XX);
	if (cudaSuccess != cudaGetLastError()) cout <<"Problem at Dot" << endl;
*/
}
#else
void Dot(Real &result, Real *x,Real *y, int M)   {
	result=0.0;
 	for (int i=0; i<M; i++) result +=x[i]*y[i];
}
#endif

#ifdef CUDA
void Sum(Real &result, Real *x, int M)   {

/**/
	Real *H_XXX=(Real*) malloc(M*sizeof(Real));
	TransferDataToHost(H_XXX, x, M);
	result=H_Sum(H_XXX,M);
	if (debug) cout <<"Host sum =" << result << endl;
	for (int i=0; i<M; i++) if (isnan(H_XXX[i])) cout <<" At "  << i << " NaN" << endl;
 	free(H_XXX);
/**/

/*	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	Real *D_XX=AllOnDev(n_blocks);
	Zero(D_XX,n_blocks);
	Real *H_XX=(Real*) malloc(n_blocks*sizeof(Real));
	sum<<<n_blocks,block_size>>>(x,D_XX,M);
	cudaThreadSynchronize();
	TransferDataToHost(H_XX, D_XX,n_blocks);
	result=0; for (int i=0; i<n_blocks; i++) {result+=H_XX[i];}
	if (cudaSuccess != cudaGetLastError()) cout <<"Problem at Sum" << endl;
	for (int i=0; i<n_blocks; i++) if (isnan(H_XX[i])) cout <<" At " << i << " NaN" << endl;

	free(H_XX);
	cudaFree(D_XX);
*/
}
#else
void Sum(Real &result, Real *x,int M)   {
if (debug) cout <<"Sum in absence of cuda" << endl;
	result=0;
 	for (int i=0; i<M; i++) result +=x[i];
}
#endif

#ifdef CUDA
void AddTimes(Real *P, Real *A, Real *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addtimes<<<n_blocks,block_size>>>(P,A,B,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at AddTimes"<<endl;}
}
#else
void AddTimes(Real *P, Real *A, Real *B, int M)   {
	for (int i=0; i<M; i++) P[i]+=A[i]*B[i];
}
#endif

#ifdef CUDA
void Times(Real *P, Real *A, Real *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	times<<<n_blocks,block_size>>>(P,A,B,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Times"<<endl;}
}
#else
void Times(Real *P, Real *A, Real *B, int M)   {
	for (int i=0; i<M; i++) P[i]=A[i]*B[i];
}
#endif
#ifdef CUDA
void Times(Real *P, Real *A, int *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	times<<<n_blocks,block_size>>>(P,A,B,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Times"<<endl;}
}
#else
void Times(Real *P, Real *A, int *B, int M)   {
	for (int i=0; i<M; i++) P[i]=A[i]*B[i];
}
#endif

#ifdef CUDA
void Norm(Real *P, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	norm<<<n_blocks,block_size>>>(P,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Norm"<<endl;}
}
#else
void Norm(Real *P, Real C, int M)   {
	for (int i=0; i<M; i++) P[i] *= C;
}
#endif

#ifdef CUDA
void Composition(Real *phi, Real *Gf, Real *Gb, Real* G1, Real C, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	composition<<<n_blocks,block_size>>>(phi,Gf,Gb,G1,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Zero"<<endl;}
}
#else
void Composition(Real *phi, Real *Gf, Real *Gb, Real* G1, Real C, int M)   {
	for (int i=0; i<M; i++) if (G1[i]>0) phi[i]+=C*Gf[i]*Gb[i]/G1[i];
}
#endif

#ifdef CUDA
void Unity(Real* P, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	unity<<<n_blocks,block_size>>>(P,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Zero"<<endl;}
}
#else
void Unity(Real* P, int M)   {
	for (int i=0; i<M; i++) P[i] =1.0;
}
#endif

#ifdef CUDA
void Zero(Real* P, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	zero<<<n_blocks,block_size>>>(P,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Zero"<<endl;}
}
#else
void Zero(Real* P, int M)   {
	for (int i=0; i<M; i++) P[i] =0.0;
}
#endif

#ifdef CUDA
void Zero(int* P, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	zero<<<n_blocks,block_size>>>(P,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Zero"<<endl;}
}
#else
void Zero(int* P, int M)   {
	for (int i=0; i<M; i++) P[i] =0;
}
#endif

#ifdef CUDA
void Cp(Real *P,Real *A, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	cp<<<n_blocks,block_size>>>(P,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Cp"<<endl;}
}
#else
void Cp(Real *P,Real *A, int M)   {
	for (int i=0; i<M; i++) P[i] = A[i];
}
#endif

#ifdef CUDA
void Cp(Real *P,int *A, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	cp<<<n_blocks,block_size>>>(P,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Cp"<<endl;}
}
#else
void Cp(Real *P,int *A, int M)   {
	for (int i=0; i<M; i++) P[i] = 1.0*A[i];
}
#endif

#ifdef CUDA
void YisAminB(Real *Y, Real *A, Real *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaminb<<<n_blocks,block_size>>>(Y,A,B,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at YisAminB"<<endl;}
}
#else
  void YisAminB(Real *Y, Real *A, Real *B, int M)   {
	for (int i=0; i<M; i++) Y[i] = A[i]-B[i];
}
#endif
#ifdef CUDA
void YisAplusC(Real *Y, Real *A, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaplusc<<<n_blocks,block_size>>>(Y,A,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at YisAplusC"<<endl;}
}
#else
void YisAplusC(Real *Y, Real *A, Real C, int M)   {
	for (int i=0; i<M; i++) Y[i] = A[i]+C;
}
#endif
#ifdef CUDA
void YisAplusB(Real *Y, Real *A, Real *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaplusb<<<n_blocks,block_size>>>(Y,A,B,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at YisAplusB"<<endl;}
}
#else
void YisAplusB(Real *Y, Real *A, Real *B, int M)   {
	for (int i=0; i<M; i++) Y[i] = A[i]+B[i];
}
#endif

#ifdef CUDA
void YplusisCtimesX(Real *Y, Real *X, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yplusisctimesx<<<n_blocks,block_size>>>(Y,X,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Uplusisctimesx"<<endl;}
}
#else
void YplusisCtimesX(Real *Y, Real *X, Real C, int M)    {
	for (int i=0; i<M; i++) Y[i] += C*X[i];
}
#endif

#ifdef CUDA
void YisAplusCtimesB(Real *Y, Real *A, Real*B, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaplusctimesb<<<n_blocks,block_size>>>(Y,A,B,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Uplusisctimesx"<<endl;}
}
#else
void YisAplusCtimesB(Real *Y, Real *A, Real*B, Real C, int M)   {
	for (int i=0; i<M; i++) Y[i]=A[i]+C*B[i];
}
#endif

#ifdef CUDA
void UpdateAlpha(Real *Y, Real *X, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	updatealpha<<<n_blocks,block_size>>>(Y,X,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at UpdateAlpha"<<endl;}
}
#else
void UpdateAlpha(Real *Y, Real *X, Real C, int M)    {
	for (int i=0; i<M; i++) Y[i] += C*(X[i]-1.0);
}
#endif
#ifdef CUDA
  void Picard(Real *Y, Real *X, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	picard<<<n_blocks,block_size>>>(Y,X,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Picard"<<endl;}
}
#else
void Picard(Real *Y, Real *X, Real C, int M)    {
	for (int i=0; i<M; i++) Y[i] = C*Y[i]+(1.0-C)*X[i];
}
#endif

#ifdef CUDA
void Add(Real *P, Real *A, int M)    {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	add<<<n_blocks,block_size>>>(P,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Add"<<endl;}
}
#else
void Add(Real *P, Real *A, int M)   {
	for (int i=0; i<M; i++) P[i]+=A[i];
}
#endif
#ifdef CUDA
void Add(int *P, int *A, int M)    {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	add<<<n_blocks,block_size>>>(P,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Add"<<endl;}
}
#else
void Add(int *P, int *A, int M)   {
	for (int i=0; i<M; i++) P[i]+=A[i];
}
#endif

#ifdef CUDA
void Dubble(Real *P, Real *A, Real norm,int M)   {
       int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	dubble<<<n_blocks,block_size>>>(P,A,norm,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Dubble"<<endl;}
}
#else
void Dubble(Real *P, Real *A, Real norm,int M)   {
	for (int i=0; i<M; i++) P[i]*=norm/A[i];
}
#endif

#ifdef CUDA
void Boltzmann(Real *P, Real *A, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	boltzmann<<<n_blocks,block_size>>>(P,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Boltzmann"<<endl;}
}
#else
void Boltzmann(Real *P, Real *A, int M)   {
	for (int i=0; i<M; i++) P[i]=exp(-A[i]);
}
#endif
#ifdef CUDA
void OverwriteC(Real *P, int *Mask, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	overwritec<<<n_blocks,block_size>>>(P,Mask,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Overwrite"<<endl;}
}
#else
void OverwriteC(Real *P, int *Mask,Real C,int M) {
	for (int i=0; i<M; i++) if (Mask[i]==1) P[i]=C; else P[i]=0;
}
#endif
#ifdef CUDA
void OverwriteA(Real *P, int *Mask, Real* A, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	overwritea<<<n_blocks,block_size>>>(P,Mask,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Overwrite"<<endl;}
}
#else
void OverwriteA(Real *P, int *Mask,Real* A,int M) {
	for (int i=0; i<M; i++) if (Mask[i]==1) P[i]=A[i]; else P[i]=0;
}
#endif

#ifdef CUDA
void UpPsi(Real* g, Real* psi, Real* X, Real* eps, int JX, int JY, Real C, int* Mask, int M)  {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	uppsi<<<n_blocks,block_size>>>(g,psi,X,eps,JX,JY,C,Mask,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at UpPsi"<<endl;}
}
#else
void UpPsi(Real* g, Real* psi, Real* X, Real* eps, int JX, int JY, Real C, int* Mask, int M) {
	cout <<"UpPsi: program should not get here, because when no CUDA-flag is used, the routine in Lattice is used. " << endl;
}
#endif

#ifdef CUDA
void UpQ(Real* g, Real* q, Real* psi, Real* eps, int JX, int JY, Real C, int* Mask, int M)  {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	upq<<<n_blocks,block_size>>>(g,q,psi,eps,JX,JY,C,Mask,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at UpQ"<<endl;}
}
#else
void UpQ(Real* g, Real* q, Real* psi, Real* eps, int JX, int JY, Real C, int* Mask, int M) {
	cout <<"UpQ: program should not get here, becauase when no CUDA-flag is used, the routine in Lattice is used. " << endl;
}
#endif


#ifdef CUDA
void Invert(Real *KSAM, Real *MASK, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	invert<<<n_blocks,block_size>>>(KSAM,MASK,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Invert"<<endl;}
}
#else
void Invert(Real *KSAM, Real *MASK, int M)   {
	for (int i=0; i<M; i++) {if (MASK[i]==0) KSAM[i]=1.0; else KSAM[i]=0.0;}
}
#endif
#ifdef CUDA
void Invert(int *KSAM, int *MASK, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	invert<<<n_blocks,block_size>>>(KSAM,MASK,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Invert"<<endl;}
}
#else
void Invert(int *KSAM, int *MASK, int M)   {
	for (int i=0; i<M; i++) {if (MASK[i]==0) KSAM[i]=1.0; else KSAM[i]=0.0;}
}
#endif
#ifdef CUDA
void AddGradSquare(Real* EE, Real* X, Real* Y, Real* Z, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addgradsquare<<<n_blocks,block_size>>>(EE,X,Y,Z,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at AddGradSquare"<<endl;}
}
#else
void AddGradSquare(Real* EE, Real* X, Real* Y, Real* Z, int M)    {
	for (int i=0; i<M; i++) EE[i] += pow(X[i]-Y[i],2)+pow(Y[i]-Z[i],2);
}
#endif

#ifdef CUDA
void PutAlpha(Real *g, Real *phitot, Real *phi_side, Real chi, Real phibulk, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	putalpha<<<n_blocks,block_size>>>(g,phitot,phi_side,chi,phibulk,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at PutAlpha"<<endl;}
}
#else
void PutAlpha(Real *g, Real *phitot, Real *phi_side, Real chi, Real phibulk, int M)   {
	for (int i=0; i<M; i++) if (phitot[i]>0) g[i] = g[i] - chi*(phi_side[i]/phitot[i]-phibulk);
}
#endif

#ifdef CUDA
void PutAlpha(Real *g, Real *phi_side, Real chi, Real phibulk, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	putalpha<<<n_blocks,block_size>>>(g,phi_side,chi,phibulk,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at PutAlpha"<<endl;}
}
#else
void PutAlpha(Real *g, Real *phi_side, Real chi, Real phibulk, int M)   {
	for (int i=0; i<M; i++) g[i] = g[i] - chi*(phi_side[i]-phibulk);
}
#endif

#ifdef CUDA
void Div(Real *P, Real *A, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	div<<<n_blocks,block_size>>>(P,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at AddDiv"<<endl;}
}
#else
void Div(Real *P, Real *A, int M)   {
	for (int i=0; i<M; i++) if (A[i]!=0) P[i]/=A[i]; else P[i]=0;
}
#endif

#ifdef CUDA
void AddG(Real *g, Real *phitot, Real *alpha, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addg<<<n_blocks,block_size>>>(g,phitot,alpha,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at AddG"<<endl;}
}
#else
void AddG(Real *g, Real *phitot, Real *alpha, int M)   {
	for (int i=0; i<M; i++) if (phitot[i]>0)  g[i]= g[i] -alpha[i] +1/phitot[i]-1.0; else g[i]=0;
}
#endif

#ifdef CUDA
void OneMinusPhitot(Real *g, Real *phitot, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	oneminusphitot<<<n_blocks,block_size>>>(g,phitot,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at OneMinusPhitot"<<endl;}
}
#else
void OneMinusPhitot(Real *g, Real *phitot, int M)   {
	for (int i=0; i<M; i++) g[i]= 1/phitot[i]-1;
}
#endif


#ifdef CUDA
void DisG1(Real* G1, Real* g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	//int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
	//distributeg1<<<n_blocks,block_size>>>(G1,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
        dim3 blocks(ceil(Mx/8.0),ceil(My/8.0),ceil(Mz/8.0));
        dim3 blockdim(8,8,8);
        distributeg1<<<blocks,blockdim>>>(G1,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);

}
#else
void DisG1(Real* G1, Real* g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int pos_l=-M;
	int pos_x,pos_y,pos_z;
	int Bxp,Byp,Bzp;
	int ii=0,jj=0,kk=0;

	for (int p=0; p<n_box; p++) { pos_l +=M; ii=0; Bxp=Bx[p]; Byp=By[p]; Bzp=Bz[p];
		for (int i=1; i<Mx+1; i++) { ii+=jx; jj=0; if (Bxp+i>MX) pos_x=(Bxp+i-MX)*JX; else pos_x = (Bxp+i)*JX;
			for (int j=1; j<My+1; j++) {jj+=jy;  kk=0; if (Byp+j>MY) pos_y=(Byp+j-MY)*JY; else pos_y = (Byp+j)*JY;
				for (int k=1; k<Mz+1; k++) { kk++; if (Bzp+k>MZ) pos_z=(Bzp+k-MZ); else pos_z = (Bzp+k);
					g1[pos_l+ii+jj+kk]=G1[pos_x+pos_y+pos_z];
				}
			}
		}
	}
}
#endif

#ifdef CUDA
void ColPhi(Real* phi, Real* GN, Real* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	//int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
	//collectphi<<<n_blocks,block_size>>>(phi,GN,rho,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
 	 dim3 blocks(ceil(Mx/8.0),ceil(My/8.0),ceil(Mz/8.0));
        dim3 blockdim(8,8,8);
        collectphi<<<blocks,blockdim>>>(phi,GN,rho,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
}
#else
void ColPhi(Real* phi, Real* GN, Real* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int pos_l=-M;
	int pos_x,pos_y,pos_z;
	int Bxp,Byp,Bzp;
	Real Inv_H_GNp;
	int ii=0,jj=0,kk=0;
	for (int p=0; p<n_box; p++) {pos_l +=M; ii=0; Bxp=Bx[p]; Byp=By[p]; Bzp=Bz[p]; Inv_H_GNp=1.0/GN[p];
		for (int i=1; i<Mx+1; i++) {ii+=jx; jj=0;  if (Bxp+i>MX) pos_x=(Bxp+i-MX)*JX; else pos_x = (Bxp+i)*JX;
			for (int j=1; j<My+1; j++) {jj+=jy;  kk=0; if (Byp+j>MY) pos_y=(Byp+j-MY)*JY; else pos_y = (Byp+j)*JY;
				for (int k=1; k<Mz+1; k++) { kk++; if (Bzp+k>MZ) pos_z=(Bzp+k-MZ); else pos_z = (Bzp+k);
					phi[pos_x+pos_y+pos_z]+=rho[pos_l+ii+jj+kk]*Inv_H_GNp;
				}
			}
		}
	}
}

#endif


//#ifdef CUDA
//void ComputeGN(Real *GN, Real *G, int M, int n_box)   {
//	int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
//	computegn<<<n_blocks,block_size>>>(GN,G,M,n_box);
//}
//#else
//void ComputeGN(Real* GN, Real* G, int M, int n_box)    {
//	for (int p=0; p<n_box; p++) GN[p]=Sum(G+p*M,M);
//}
//#endif

void H_Zero(float* H_P, int M)   {//this precedure should act on a host P.
	for (int i=0; i<M; i++) H_P[i] = 0;
}
void H_Zero(Real* H_P, int M)   {//this precedure should act on a host P.
	for (int i=0; i<M; i++) H_P[i] = 0;
}
void H_Zero(int* H_P, int M)   {//this precedure should act on a host P.
	for (int i=0; i<M; i++) H_P[i] = 0;
}

Real H_Sum(Real* H, int M){
	Real Sum=0;
	for (int i=0; i<M; i++) Sum+=H[i];
	return Sum;
}
Real H_Dot(Real* A, Real *B, int M){
	Real Sum=0;
	for (int i=0; i<M; i++) Sum+=A[i]*B[i];
	return Sum;
}

void H_Invert(int* KSAM, int *MASK, int M)   { //only necessary on CPU
	for (int z=0; z<M; z++) {if (MASK[z]==0) KSAM[z]=1.0; else KSAM[z]=0.0;}
}

#ifdef CUDA
void SetBoundaries(Real *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)   {
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	bx<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at SetBoundaries"<<endl;}
}
#else
void SetBoundaries(Real *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)    {
	bx(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#endif
#ifdef CUDA
void SetBoundaries(int *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)   {
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	bx<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at SetBoundaries"<<endl;}
}
#else
void SetBoundaries(int *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)    {
	bx(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#endif


#ifdef CUDA
void RemoveBoundaries(Real *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)   {
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	b_x<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at RemoveBoundaries"<<endl;}
}
#else
void RemoveBoundaries(Real *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)    {
	b_x(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#endif
#ifdef CUDA
void RemoveBoundaries(int *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)   {
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	b_x<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at RemoveBoundaries"<<endl;}
}
#else
void RemoveBoundaries(int *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)    {
	b_x(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#endif

#define MAX_ITER 100
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(x,y) ((x)>(y)?(x):(y))

Real PYTHAG(Real a, Real b)
{
    Real at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}

int svdcmp(Real** a, int m, int n, Real *w, Real** v)
{
    int flag, i, its, j, jj, k, l, nm;
    Real c, f, h, s, x, y, z;
    Real anorm = 0.0, g = 0.0, scale = 0.0;
    Real *rv1;

    if (m < n)
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }

    rv1 = (Real *)malloc((unsigned int) n*sizeof(Real));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs((Real)a[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    a[k][i] = (Real)((Real)a[k][i]/scale);
                    s += ((Real)a[k][i] * (Real)a[k][i]);
                }
                f = (Real)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (Real)(f - g);
                if (i != n - 1)
                {
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += ((Real)a[k][i] * (Real)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            a[k][j] += (Real)(f * (Real)a[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    a[k][i] = (Real)((Real)a[k][i]*scale);
            }
        }
        w[i] = (Real)(scale * g);

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1)
        {
            for (k = l; k < n; k++)
                scale += fabs((Real)a[i][k]);
            if (scale)
            {
                for (k = l; k < n; k++)
                {
                    a[i][k] = (Real)((Real)a[i][k]/scale);
                    s += ((Real)a[i][k] * (Real)a[i][k]);
                }
                f = (Real)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (Real)(f - g);
                for (k = l; k < n; k++)
                    rv1[k] = (Real)a[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < n; k++)
                            s += ((Real)a[j][k] * (Real)a[i][k]);
                        for (k = l; k < n; k++)
                            a[j][k] += (Real)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++)
                    a[i][k] = (Real)((Real)a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((Real)w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        if (i < n - 1)
        {
            if (g)
            {
                for (j = l; j < n; j++)
                    v[j][i] = (Real)(((Real)a[i][j] / (Real)a[i][l]) / g);
                    /* Real division to avoid underflow */
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < n; k++)
                        s += ((Real)a[i][k] * (Real)v[k][j]);
                    for (k = l; k < n; k++)
                        v[k][j] += (Real)(s * (Real)v[k][i]);
                }
            }
            for (j = l; j < n; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        l = i + 1;
        g = (Real)w[i];
        if (i < n - 1)
            for (j = l; j < n; j++)
                a[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != n - 1)
            {
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += ((Real)a[k][i] * (Real)a[k][j]);
                    f = (s / (Real)a[i][i]) * g;
                    for (k = i; k < m; k++)
                        a[k][j] += (Real)(f * (Real)a[k][i]);
                }
            }
            for (j = i; j < m; j++)
                a[j][i] = (Real)((Real)a[j][i]*g);
        }
        else
        {
            for (j = i; j < m; j++)
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs((Real)w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = (Real)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (Real)h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = (Real)a[j][nm];
                            z = (Real)a[j][i];
                            a[j][nm] = (Real)(y * c + z * s);
                            a[j][i] = (Real)(z * c - y * s);
                        }
                    }
                }
            }
            z = (Real)w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (Real)(-z);
                    for (j = 0; j < n; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30! iterations \n");
                return(0);
            }

            /* shift from bottom 2 x 2 minor */
            x = (Real)w[l];
            nm = k - 1;
            y = (Real)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = (Real)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++)
                {
                    x = (Real)v[jj][j];
                    z = (Real)v[jj][i];
                    v[jj][j] = (Real)(x * c + z * s);
                    v[jj][i] = (Real)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (Real)z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = (Real)a[jj][j];
                    z = (Real)a[jj][i];
                    a[jj][j] = (Real)(y * c + z * s);
                    a[jj][i] = (Real)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (Real)x;
        }
    }
    free((void*) rv1);
    return(1);
}

int modern_svdcmp(Real** a, int m, int n, Real *w, Real** v)
{
    if (debug) cout << "modern_svdcmp in tools" << endl;
    int flag, i, its, j, jj, k, l, nm;
    Real c, f, h, s, x, y, z;
    Real anorm = 0.0, g = 0.0, scale = 0.0;

    if (m < n)
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }

    vector<Real> rv1(n);

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1.push_back( scale * g );

        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs((Real)a[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    a[k][i] = (Real)((Real)a[k][i]/scale);
                    s += ((Real)a[k][i] * (Real)a[k][i]);
                }
                f = (Real)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (Real)(f - g);
                if (i != n - 1)
                {
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += ((Real)a[k][i] * (Real)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            a[k][j] += (Real)(f * (Real)a[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    a[k][i] = (Real)((Real)a[k][i]*scale);
            }
        }
        w[i] = (Real)(scale * g);

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1)
        {
            for (k = l; k < n; k++)
                scale += fabs((Real)a[i][k]);
            if (scale)
            {
                for (k = l; k < n; k++)
                {
                    a[i][k] = (Real)((Real)a[i][k]/scale);
                    s += ((Real)a[i][k] * (Real)a[i][k]);
                }
                f = (Real)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (Real)(f - g);
                for (k = l; k < n; k++)
                    rv1[k] = (Real)a[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < n; k++)
                            s += ((Real)a[j][k] * (Real)a[i][k]);
                        for (k = l; k < n; k++)
                            a[j][k] += (Real)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++)
                    a[i][k] = (Real)((Real)a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((Real)w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        if (i < n - 1)
        {
            if (g)
            {
                for (j = l; j < n; j++){
                    v[j][i] = (Real)(((Real)a[i][j] / (Real)a[i][l]) / g);
                  }
                    /* Real division to avoid underflow */
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < n; k++)
                        s += ((Real)a[i][k] * (Real)v[k][j]);
                    for (k = l; k < n; k++)
                        v[k][j] += (Real)(s * (Real)v[k][i]);
                }
            }
            for (j = l; j < n; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        l = i + 1;
        g = (Real)w[i];
        if (i < n - 1)
            for (j = l; j < n; j++)
                a[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != n - 1)
            {
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += ((Real)a[k][i] * (Real)a[k][j]);
                    f = (s / (Real)a[i][i]) * g;
                    for (k = i; k < m; k++)
                        a[k][j] += (Real)(f * (Real)a[k][i]);
                }
            }
            for (j = i; j < m; j++)
                a[j][i] = (Real)((Real)a[j][i]*g);
        }
        else
        {
            for (j = i; j < m; j++)
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs((Real)w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = (Real)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (Real)h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = (Real)a[j][nm];
                            z = (Real)a[j][i];
                            a[j][nm] = (Real)(y * c + z * s);
                            a[j][i] = (Real)(z * c - y * s);
                        }
                    }
                }
            }
            z = (Real)w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (Real)(-z);
                    for (j = 0; j < n; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                fprintf(stderr, "No convergence after 30! iterations \n");
                return 0;
            }

            /* shift from bottom 2 x 2 minor */
            x = (Real)w[l];
            nm = k - 1;
            y = (Real)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = (Real)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++)
                {
                    x = (Real)v[jj][j];
                    z = (Real)v[jj][i];
                    v[jj][j] = (Real)(x * c + z * s);
                    v[jj][i] = (Real)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (Real)z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = (Real)a[jj][j];
                    z = (Real)a[jj][i];
                    a[jj][j] = (Real)(y * c + z * s);
                    a[jj][i] = (Real)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (Real)x;
        }
    }
    return 0;
}
