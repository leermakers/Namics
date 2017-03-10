#include "tools.h"  
#include "namics.h" 
#define MAX_BLOCK_SZ 512

   
#ifdef CUDA
//cublasHandle_t handle;
//cublasStatus_t stat=cublasCreate(&handle);
const int block_size = 256;
/*
__global__ void distributeg1(rene *G1, rene *g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY)   {
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

__global__ void collectphi(rene* phi, rene* GN,rene* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY)   {
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	int j = blockIdx.y*blockDim.y+threadIdx.y;
	int k = blockIdx.z*blockDim.z+threadIdx.z;
	int pos_r=(JX+JY+1);//This needs af fix**
	int pM=jx+jy+1+jx*i+jy*j+k;
	int ii,jj,kk;//This needs af fix**
	int MXm1 = MX-1;
	int MYm1 = MY-1;
	int MZm1 = MZ-1;
	if (k < Mz && j < My && i < Mx){
		for (int p = 0;p<n_box;++p){
			if (Bx[p]+i > MXm1)  ii=(Bx[p]+i-MX)*JX; else ii=(Bx[p]+i)*JX;
			if (By[p]+j > MYm1)  jj=(By[p]+j-MY)*JY; else jj=(By[p]+j)*JY;
			if (Bz[p]+k > MZm1)  kk=(Bz[p]+k-MZ); else kk=(Bz[p]+k);

			atomicAdd(&phi[pos_r+ii+jj+kk], rho[pM]/GN[p]); //This needs af fix*****************************
			pM+=M;
		}
	}
}
*/

__global__ void dot(rene *X, rene *Y, rene *Z, int M)   {
   __shared__ rene tmp[MAX_BLOCK_SZ];
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

__global__ void sum(rene *X, rene *Z, int M)   {
   __shared__ rene tmp[block_size];
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


__global__ void composition(rene *phi, rene *Gf, rene *Gb, rene* G1, rene C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (G1[idx]>0) phi[idx] += C*Gf[idx]*Gb[idx]/G1[idx];
}
__global__ void times(rene *P, rene *A, rene *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void times(rene *P, rene *A, int *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void addtimes(rene *P, rene *A, rene *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx]*B[idx];
}
__global__ void norm(rene *P, rene C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] *= C;
}
__global__ void zero(rene *P, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0.0;
}
__global__ void unity(rene *P, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 1.0;
}
__global__ void zero(int *P, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0;
}
__global__ void cp (rene *P, rene *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = A[idx];
}
__global__ void cp (rene *P, int *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 1.0*A[idx];
}

__global__ void yisaplusctimesb(rene *Y, rene *A,rene *B, rene C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]+ C * B[idx];
}
__global__ void yisaminb(rene *Y, rene *A,rene *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]-B[idx];
}
__global__ void yisaplusc(rene *Y, rene *A, rene C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]+C;
}
__global__ void yisaplusb(rene *Y, rene *A,rene *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]+B[idx];
}
__global__ void yplusisctimesx(rene *Y, rene *X, rene C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*X[idx];
}
__global__ void updatealpha(rene *Y, rene *X, rene C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*(X[idx]-1.0);
}
__global__ void picard(rene *Y, rene *X, rene C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = C*Y[idx]+(1-C)*X[idx];
}
__global__ void add(rene *P, rene *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx];
}
__global__ void add(int *P, int *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx];
}
__global__ void dubble(rene *P, rene *A, rene norm, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]*=norm/A[idx];
}
__global__ void boltzmann(rene *P, rene *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=exp(-A[idx]);
}
__global__ void invert(int *SKAM, int *MASK, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) SKAM[idx]=(MASK[idx]-1)*(MASK[idx]-1);
}
__global__ void invert(rene *SKAM, rene *MASK, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) SKAM[idx]=(MASK[idx]-1)*(MASK[idx]-1);
}
__global__ void putalpha(rene *g,rene *phitot,rene *phi_side,rene chi,rene phibulk,int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (phitot[idx]>0) g[idx] = g[idx] - chi*(phi_side[idx]/phitot[idx]-phibulk);
}
__global__ void div(rene *P,rene *A,int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (A[idx]>0) P[idx] /=A[idx];
}
__global__ void oneminusphitot(rene *g, rene *phitot, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) g[idx]= 1/phitot[idx]-1;
}
__global__ void addg(rene *g, rene *phitot, rene *alpha, int M)    {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {
		if (phitot[idx]>0) g[idx]= g[idx] -alpha[idx] +1/phitot[idx]-1; else g[idx]=0;
	}
}

//__global__ void computegn(rene *GN, rene* G, int M, int n_box)   {
//	int p= blockIdx.x*blockDim.x+threadIdx.x;
//	if (p<n_box){
//		cublasDasum(handle,M,G+p*M,1,GN+p); //in case of float we need the cublasSasum etcetera.
//	}
//}

void TransferDataToHost(rene *H, rene *D, int M)    {
	cudaMemcpy(H, D, sizeof(rene)*M,cudaMemcpyDeviceToHost);
}
void TransferDataToDevice(rene *H, rene *D, int M)    {
	cudaMemcpy(D, H, sizeof(rene)*M,cudaMemcpyHostToDevice);
}
void TransferIntDataToHost(int *H, int *D, int M)   {
	cudaMemcpy(H, D, sizeof(int)*M,cudaMemcpyDeviceToHost);
}
void TransferIntDataToDevice(int *H, int *D, int M)     {
	cudaMemcpy(D, H, sizeof(int)*M,cudaMemcpyHostToDevice);
}
#endif

#ifdef CUDA
__global__ void bx(rene *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx, jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=P[bx1_jx+idx];
		P[jx_mmx+idx]=P[jx_bxm+idx];
	}
}
__global__ void b_x(rene *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=0;
		P[jx_mmx+idx]=0;
	}
}
__global__ void by(rene *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=P[jy_by1+idx];
		P[jy_mmy+idx]=P[jy_bym+idx];
	}
}
__global__ void b_y(rene *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=0;
		P[jy_mmy+idx]=0;
	}
}
__global__ void bz(rene *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=P[idx+bz1];
		P[idx+mmz]=P[idx+bzm];
	}
}
__global__ void b_z(rene *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
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
void bx(rene *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
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
void b_x(rene *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}
void by(rene *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
	}
}
void b_y(rene *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}
void bz(rene *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
	}
}
void b_z(rene *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
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
rene *AllOnDev(int N)    {
	rene *X;
	if (cudaSuccess != cudaMalloc((void **) &X, sizeof(rene)*N))
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
rene *AllOnDev(int N)    {
	rene *X=NULL;
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
void Dot(rene &result, rene *x,rene *y, int M)   {
/**/
	rene *H_XXX=(rene*) malloc(M*sizeof(rene)); 
	rene *H_YYY=(rene*) malloc(M*sizeof(rene)); 
	TransferDataToHost(H_XXX, x, M);
	TransferDataToHost(H_YYY, y, M);
	result=H_Dot(H_XXX,H_YYY,M);
	free(H_XXX);
	free(H_YYY);
/**/


/*	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	rene *D_XX=AllOnDev(n_blocks);
	Zero(D_XX,n_blocks);
	rene *H_XX=(rene*) malloc(n_blocks*sizeof(rene)); 
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
void Dot(rene &result, rene *x,rene *y, int M)   {
	result=0.0;
 	for (int i=0; i<M; i++) result +=x[i]*y[i];
}
#endif

#ifdef CUDA
void Sum(rene &result, rene *x, int M)   {

/**/	
	rene *H_XXX=(rene*) malloc(M*sizeof(rene)); 
	TransferDataToHost(H_XXX, x, M);
	result=H_Sum(H_XXX,M);
	if (debug) cout <<"Host sum =" << result << endl;	
	for (int i=0; i<M; i++) if (isnan(H_XXX[i])) cout <<" At "  << i << " NaN" << endl; 
 	free(H_XXX); 
/**/

/*	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	rene *D_XX=AllOnDev(n_blocks);
	Zero(D_XX,n_blocks); 
	rene *H_XX=(rene*) malloc(n_blocks*sizeof(rene)); 
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
void Sum(rene &result, rene *x,int M)   {
if (debug) cout <<"Sum in absence of cuda" << endl; 
	result=0;
 	for (int i=0; i<M; i++) result +=x[i];
}
#endif

#ifdef CUDA
void AddTimes(rene *P, rene *A, rene *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addtimes<<<n_blocks,block_size>>>(P,A,B,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at AddTimes"<<endl;}
}
#else
void AddTimes(rene *P, rene *A, rene *B, int M)   {
	for (int i=0; i<M; i++) P[i]+=A[i]*B[i];
}
#endif

#ifdef CUDA
void Times(rene *P, rene *A, rene *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	times<<<n_blocks,block_size>>>(P,A,B,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Times"<<endl;}
}
#else
void Times(rene *P, rene *A, rene *B, int M)   {
	for (int i=0; i<M; i++) P[i]=A[i]*B[i];
}
#endif
#ifdef CUDA
void Times(rene *P, rene *A, int *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	times<<<n_blocks,block_size>>>(P,A,B,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Times"<<endl;}
}
#else
void Times(rene *P, rene *A, int *B, int M)   {
	for (int i=0; i<M; i++) P[i]=A[i]*B[i];
}
#endif

#ifdef CUDA
void Norm(rene *P, rene C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	norm<<<n_blocks,block_size>>>(P,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Norm"<<endl;}
}
#else
void Norm(rene *P, rene C, int M)   {
	for (int i=0; i<M; i++) P[i] *= C;
}
#endif

#ifdef CUDA
void Composition(rene *phi, rene *Gf, rene *Gb, rene* G1, rene C, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	composition<<<n_blocks,block_size>>>(phi,Gf,Gb,G1,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Zero"<<endl;}
}
#else
void Composition(rene *phi, rene *Gf, rene *Gb, rene* G1, rene C, int M)   {
	for (int i=0; i<M; i++) if (G1[i]>0) phi[i]+=C*Gf[i]*Gb[i]/G1[i];
}
#endif

#ifdef CUDA
void Unity(rene* P, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	unity<<<n_blocks,block_size>>>(P,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Zero"<<endl;}
}
#else
void Unity(rene* P, int M)   {
	for (int i=0; i<M; i++) P[i] =1.0;
}
#endif

#ifdef CUDA
void Zero(rene* P, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	zero<<<n_blocks,block_size>>>(P,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Zero"<<endl;}
}
#else
void Zero(rene* P, int M)   {
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
void Cp(rene *P,rene *A, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	cp<<<n_blocks,block_size>>>(P,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Cp"<<endl;}
}
#else
void Cp(rene *P,rene *A, int M)   {
	for (int i=0; i<M; i++) P[i] = A[i];
}
#endif

#ifdef CUDA
void Cp(rene *P,int *A, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	cp<<<n_blocks,block_size>>>(P,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Cp"<<endl;}
}
#else
void Cp(rene *P,int *A, int M)   {
	for (int i=0; i<M; i++) P[i] = 1.0*A[i];
}
#endif

#ifdef CUDA
void YisAminB(rene *Y, rene *A, rene *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaminb<<<n_blocks,block_size>>>(Y,A,B,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at YisAminB"<<endl;}
}
#else
  void YisAminB(rene *Y, rene *A, rene *B, int M)   {
	for (int i=0; i<M; i++) Y[i] = A[i]-B[i];
}
#endif
#ifdef CUDA
void YisAplusC(rene *Y, rene *A, rene C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaplusc<<<n_blocks,block_size>>>(Y,A,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at YisAplusC"<<endl;}
}
#else
void YisAplusC(rene *Y, rene *A, rene C, int M)   {
	for (int i=0; i<M; i++) Y[i] = A[i]+C;
}
#endif
#ifdef CUDA
void YisAplusB(rene *Y, rene *A, rene *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaplusb<<<n_blocks,block_size>>>(Y,A,B,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at YisAplusB"<<endl;}
}
#else
void YisAplusB(rene *Y, rene *A, rene *B, int M)   {
	for (int i=0; i<M; i++) Y[i] = A[i]+B[i];
}
#endif

#ifdef CUDA
void YplusisCtimesX(rene *Y, rene *X, rene C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yplusisctimesx<<<n_blocks,block_size>>>(Y,X,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Uplusisctimesx"<<endl;}
}
#else
void YplusisCtimesX(rene *Y, rene *X, rene C, int M)    {
	for (int i=0; i<M; i++) Y[i] += C*X[i];
}
#endif

#ifdef CUDA
void YisAplusCtimesB(rene *Y, rene *A, rene*B, rene C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaplusctimesb<<<n_blocks,block_size>>>(Y,A,B,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Uplusisctimesx"<<endl;}
}
#else
void YisAplusCtimesB(rene *Y, rene *A, rene*B, rene C, int M)   {
	for (int i=0; i<M; i++) Y[i]=A[i]+C*B[i];
}
#endif

#ifdef CUDA
void UpdateAlpha(rene *Y, rene *X, rene C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	updatealpha<<<n_blocks,block_size>>>(Y,X,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at UpdateAlpha"<<endl;}
}
#else
void UpdateAlpha(rene *Y, rene *X, rene C, int M)    {
	for (int i=0; i<M; i++) Y[i] += C*(X[i]-1.0);
}
#endif
#ifdef CUDA
  void Picard(rene *Y, rene *X, rene C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	picard<<<n_blocks,block_size>>>(Y,X,C,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Picard"<<endl;}
}
#else
void Picard(rene *Y, rene *X, rene C, int M)    {
	for (int i=0; i<M; i++) Y[i] = C*Y[i]+(1.0-C)*X[i];
}
#endif

#ifdef CUDA
void Add(rene *P, rene *A, int M)    {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	add<<<n_blocks,block_size>>>(P,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Add"<<endl;}
}
#else
void Add(rene *P, rene *A, int M)   {
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
void Dubble(rene *P, rene *A, rene norm,int M)   {
       int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	dubble<<<n_blocks,block_size>>>(P,A,norm,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Dubble"<<endl;}
}
#else
void Dubble(rene *P, rene *A, rene norm,int M)   {
	for (int i=0; i<M; i++) P[i]*=norm/A[i];
}
#endif

#ifdef CUDA
void Boltzmann(rene *P, rene *A, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	boltzmann<<<n_blocks,block_size>>>(P,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Boltzmann"<<endl;}
}
#else
void Boltzmann(rene *P, rene *A, int M)   {
	for (int i=0; i<M; i++) P[i]=exp(-A[i]);
}
#endif
#ifdef CUDA
void Invert(rene *KSAM, rene *MASK, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	invert<<<n_blocks,block_size>>>(KSAM,MASK,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at Invert"<<endl;}
}
#else
void Invert(rene *KSAM, rene *MASK, int M)   {
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
void PutAlpha(rene *g, rene *phitot, rene *phi_side, rene chi, rene phibulk, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	putalpha<<<n_blocks,block_size>>>(g,phitot,phi_side,chi,phibulk,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at PutAlpha"<<endl;}
}
#else
void PutAlpha(rene *g, rene *phitot, rene *phi_side, rene chi, rene phibulk, int M)   {
	for (int i=0; i<M; i++) if (phitot[i]>0) g[i] = g[i] - chi*(phi_side[i]/phitot[i]-phibulk);
}
#endif

#ifdef CUDA
void Div(rene *P, rene *A, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	div<<<n_blocks,block_size>>>(P,A,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at AddDiv"<<endl;}
}
#else
void Div(rene *P, rene *A, int M)   {
	for (int i=0; i<M; i++) if (A[i]!=0) P[i]/=A[i]; else P[i]=0;
}
#endif

#ifdef CUDA
void AddG(rene *g, rene *phitot, rene *alpha, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addg<<<n_blocks,block_size>>>(g,phitot,alpha,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at AddG"<<endl;}
}
#else
void AddG(rene *g, rene *phitot, rene *alpha, int M)   {
	for (int i=0; i<M; i++) if (phitot[i]>0)  g[i]= g[i] -alpha[i] +1/phitot[i]-1.0; else g[i]=0;
}
#endif

#ifdef CUDA
void OneMinusPhitot(rene *g, rene *phitot, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	oneminusphitot<<<n_blocks,block_size>>>(g,phitot,M);
	if (cudaSuccess != cudaGetLastError()) {cout <<"problem at OneMinusPhitot"<<endl;}
}
#else
void OneMinusPhitot(rene *g, rene *phitot, int M)   {
	for (int i=0; i<M; i++) g[i]= 1/phitot[i]-1;
}
#endif

//#ifdef CUDA
//void ComputeGN(rene *GN, rene *G, int M, int n_box)   {
//	int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
//	computegn<<<n_blocks,block_size>>>(GN,G,M,n_box);
//}
//#else
//void ComputeGN(rene* GN, rene* G, int M, int n_box)    {
//	for (int p=0; p<n_box; p++) GN[p]=Sum(G+p*M,M); 
//}
//#endif

void H_Zero(float* H_P, int M)   {//this precedure should act on a host P.
	for (int i=0; i<M; i++) H_P[i] = 0;
}
void H_Zero(rene* H_P, int M)   {//this precedure should act on a host P.
	for (int i=0; i<M; i++) H_P[i] = 0;
}
void H_Zero(int* H_P, int M)   {//this precedure should act on a host P.
	for (int i=0; i<M; i++) H_P[i] = 0;
}

rene H_Sum(rene* H, int M){
	rene Sum=0;
	for (int i=0; i<M; i++) Sum+=H[i];
	return Sum;
}
rene H_Dot(rene* A, rene *B, int M){
	rene Sum=0;
	for (int i=0; i<M; i++) Sum+=A[i]*B[i];
	return Sum;
}

void H_Invert(int* KSAM, int *MASK, int M)   { //only necessary on CPU
	for (int z=0; z<M; z++) {if (MASK[z]==0) KSAM[z]=1.0; else KSAM[z]=0.0;}
}

#ifdef CUDA
void SetBoundaries(rene *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)   {
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
void SetBoundaries(rene *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)    {
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
void RemoveBoundaries(rene *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)   {
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
void RemoveBoundaries(rene *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)    {
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
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static rene maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))
static rene sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

rene pythag(rene a, rene b)
//Computes (a 2 + b 2 ) 1/2 without destructive underflow or overflow.
{
    rene absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


void svdcmp(rene **a, int m, int n, rene *w, rene **v)
/**Given a matrix a[1..m][1..n] , this routine computes its singular value decomposition, A =
U ·W ·V T . The matrix U replaces a on output. The diagonal matrix of singular values W is out-
put as a vector w[1..n] . The matrix V (not the transpose V T ) is output as v[1..n][1..n] .
*/
{
    //rene pythag(rene a, rene b);
    
    int flag,i,its,j,jj,k,l,nm;
    rene anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  
    rv1 = new rene[n+1]; //vector(1,n);
    
    g=scale=anorm=0.0;
    //Householder reduction to bidiagonal form.
    for (i=1;i<=n;i++) {
        l=i+1;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i <= m) {
            for (k=i;k<=m;k++) scale += fabs(a[k][i]);
            if (scale) {
                for (k=i;k<=m;k++) {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f=a[i][i];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;
                for (j=l;j<=n;j++) {
                    for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
                    f=s/h;
                    for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                }
                for (k=i;k<=m;k++) a[k][i] *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
    
        if (i <= m && i != n) {
            for (k=l;k<=n;k++) scale += fabs(a[i][k]);
            if (scale) {
                for (k=l;k<=n;k++) {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f=a[i][l];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][l]=f-g;
                for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
                for (j=l;j<=m;j++) {
                    for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
                    for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
                }
                for (k=l;k<=n;k++) a[i][k] *= scale;
            }
        }
        anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n;i>=1;i--) { //Accumulation of right-hand transformations.
        if (i < n) {
            if (g) {
                for (j=l;j<=n;j++)
                    //Double division to avoid possible underflow.
                    v[j][i]=(a[i][j]/a[i][l])/g;
                for (j=l;j<=n;j++) {
                    for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
                    for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
                }
            }
            for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
        }
        v[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
        
    for (i=IMIN(m,n);i>=1;i--) { //Accumulation of left-hand transformations.
        l=i+1;
        g=w[i];
        for (j=l;j<=n;j++) a[i][j]=0.0;
        if (g) {
            g=1.0/g;
            for (j=l;j<=n;j++) {
                for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
                f=(s/a[i][i])*g;
                for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
            }
            for (j=i;j<=m;j++) a[j][i] *= g;
        } else for (j=i;j<=m;j++) a[j][i]=0.0;
        ++a[i][i];
    }     
    
    for (k=n;k>=1;k--) { //Diagonalization of the bidiagonal form: Loop over
        for (its=1;its<=MAX_ITER;its++) { //singular values, and over allowed iterations.
            flag=1; //Test for splitting.
            for (l=k;l>=1;l--) {
                nm=l-1; //Note that rv1[1] is always zero.
                if ((rene)(fabs(rv1[l])+anorm) == anorm) {
                    flag=0;
                    break;
                }
                if ((rene)(fabs(w[nm])+anorm) == anorm) break;
            }
            if (flag) {
                c=0.0; //Cancellation of rv1[l], if l > 1.
                s=1.0;
                for (i=l;i<=k;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if ((rene)(fabs(f)+anorm) == anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=1;j<=m;j++) {
                        y=a[j][nm];
                        z=a[j][i];
                        a[j][nm]=y*c+z*s;
                        a[j][i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k) { //Convergence.
                if (z < 0.0) { //Singular value is made nonnegative.
                    w[k] = -z;
                    for (j=1;j<=n;j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == MAX_ITER) printf("no convergence in %d svdcmp iterations\n", MAX_ITER);
                              //nrerror("no convergence in 30 svdcmp iterations");
            x=w[l]; //Shift from bottom 2-by-2 minor.
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x; 
            c=s=1.0; //Next QR transformation:   
        
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;  
                s=h/z;
                f=x*c+g*s;
                g = g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=1;jj<=n;jj++) {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z;
                //Rotation can be arbitrary if z = 0.
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=1;jj<=m;jj++) {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }
        
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    
    }
    //free_vector(rv1,1,n);
    
    delete [] rv1;
}



