

#ifdef CUDA

float *BlasResult;
cublasStatus_t stat;
cublasHandle_t handle;
int block_size=256;

__global__ void times(float *P, float *A, float *B, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void addtimes(float *P, float *A, float *B, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx]*B[idx];
}
__global__ void norm(float *P, float C, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] *= C;
}
__global__ void zero(float *P, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0.0;
}
__global__ void cp (float *P, float *A, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = A[idx];
}
__global__ void yisaminb(float *Y, float *A,float *B, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]-B[idx];
}
__global__ void yplusisctimesx(float *Y, float *X, float C, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*X[idx];
}
__global__ void add(float *P, float *A, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx];
}
__global__ void dubble(float *P, float *A, float norm, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]*=norm/A[idx];
}
__global__ void boltzmann(float *P, float *A, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=exp(-A[idx]);
}
__global__ void putalpha(float *g,float *phitot,float *phi_side,float chi,float phibulk,int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) g[idx] = g[idx] - chi*(phi_side[idx]/phitot[idx]-phibulk);
}
__global__ void div(float *P,float *A,int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (A[idx]>0) P[idx] /=A[idx];
}
__global__ void oneminusphitot(float *g, float *phitot, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) g[idx]= 1/phitot[idx]-1;
}
__global__ void addg(float *g, float *phitot, float *alpha, int M) {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {
		g[idx]= g[idx] -alpha[idx] +1/phitot[idx]-1;
	}
}

__global__ void computegn(float *GN, float* G, int M, int n_box){
	int p= blockIdx.x*blockDim.x+threadIdx.x;
	if (p<n_box){
		cublasSasum(handle,M,G+p*M,1,GN[p]);
	}
}

void TransferDataToHost(float *H, float *D, int M) {
	cudaMemcpy(H, D, sizeof(float)*M,cudaMemcpyDeviceToHost);
}
void TransferDataToDevice(float *H, float *D, int M) {
	cudaMemcpy(D, H, sizeof(float)*M,cudaMemcpyHostToDevice);
}
void TransferIntDataToHost(int *H, int *D, int M) {
	cudaMemcpy(H, D, sizeof(int)*M,cudaMemcpyDeviceToHost);
}
void TransferIntDataToDevice(int *H, int *D, int M) {
	cudaMemcpy(D, H, sizeof(int)*M,cudaMemcpyHostToDevice);
}
#endif

#ifdef CUDA
__global__ void bx(float *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
	int idx, jx_mmx=jx*mmx, jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=P[bx1_jx+idx];
		P[jx_mmx+idx]=P[jx_bxm+idx];
	}
}
__global__ void b_x(float *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
	int idx, jx_mmx=jx*mmx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=0;
		P[jx_mmx+idx]=0;
	}
}
__global__ void by(float *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int idx, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=P[jy_by1+idx];
		P[jy_mmy+idx]=P[jy_bym+idx];
	}
}
__global__ void b_y(float *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int idx, jy_mmy=jy*mmy;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=0;
		P[jy_mmy+idx]=0;
	}
}
__global__ void bz(float *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=P[idx+bz1];
		P[idx+mmz]=P[idx+bzm];
	}
}
__global__ void b_z(float *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=0;
		P[idx+mmz]=0;
	}
}


#else
void bx(float *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
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
void b_x(float *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}
void by(float *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
	}
}
void b_y(float *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}
void bz(float *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
	}
}
void b_z(float *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
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
float Dot(float *x,float *y,int M){
	float result;
 	cublasSdot(handle,M,x,1,y,1,BlasResult);
	cudaMemcpy(&result,BlasResult,sizeof(float),cudaMemcpyDeviceToHost);
	return result;
}
#else
float Dot(float *x,float *y,int M){
	float result=0;
 	for (int i=0; i<M; i++) result +=x[i]*y[i];
	return result;
}
#endif

#ifdef CUDA
float Sum(float *x,int M){
	float result;
 	cublasSasum(handle,M,x,1,BlasResult);
	cudaMemcpy(&result,BlasResult,sizeof(float),cudaMemcpyDeviceToHost);
	return result;
}
#else
float Sum(float *x,int M){
	float result=0;
 	for (int i=0; i<M; i++) result +=x[i];
	return result;
}
#endif

#ifdef CUDA
bool GPU_present() {
int deviceCount =0; cuDeviceGetCount(&deviceCount);
    if (deviceCount ==0) printf("There is no device supporting Cuda.\n");
    else cudaDeviceReset();
	return deviceCount > 0;
}
#else

bool GPU_present() {
	return false;
}
#endif


#ifdef CUDA
float *AllOnDev(int N) {
	float *X;
	if (cudaSuccess != cudaMalloc((void **) &X, sizeof(float)*N))
	printf("Memory allocation on GPU failed.\n Please reduce size of lattice and/or chain length(s) or turn on CUDA\n");
	return X;
}
int *AllIntOnDev(int N) {
	int *X;
	if (cudaSuccess != cudaMalloc((void **) &X, sizeof(int)*N))
	printf("Memory allocation on GPU failed.\n Please reduce size of lattice and/or chain length(s) or turn on CUDA\n");
	return X;
}
#else
float *AllOnDev(int N) {
	float *X=NULL;
	printf("Turn on CUDA\n");
	return X;
}
int *AllIntOnDev(int N) {
	int *X=NULL;
	printf("Turn on CUDA\n");
	return X;
}
#endif

#ifdef CUDA
void AddTimes(float *P, float *A, float *B, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addtimes<<<n_blocks,block_size>>>(P,A,B,M);
}
#else
void AddTimes(float *P, float *A, float *B, int M){
	for (int i=0; i<M; i++) P[i]+=A[i]*B[i];
}
#endif

#ifdef CUDA
void Times(float *P, float *A, float *B, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	times<<<n_blocks,block_size>>>(P,A,B,M);
}
#else
void Times(float *P, float *A, float *B, int M){
	for (int i=0; i<M; i++) P[i]=A[i]*B[i];
}
#endif

#ifdef CUDA
void Norm(float *P, float C, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	norm<<<n_blocks,block_size>>>(P,C,M);
}
#else
void Norm(float *P, float C, int M){
	for (int i=0; i<M; i++) P[i] *= C;
}
#endif

#ifdef CUDA
void Zero(float* P, int M){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	zero<<<n_blocks,block_size>>>(P,M);
}
#else
void Zero(float* P, int M){
	for (int i=0; i<M; i++) P[i] =0;
}
#endif

#ifdef CUDA
void Cp(float *P,float *A, int M){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	cp<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Cp(float *P,float *A, int M){
	for (int i=0; i<M; i++) P[i] = A[i];
}
#endif

#ifdef CUDA
void YisAminB(float *Y, float *A, float *B, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaminb<<<n_blocks,block_size>>>(Y,A,B,M);
}
#else
void YisAminB(float *Y, float *A, float *B, int M){
	for (int i=0; i<M; i++) Y[i] = A[i]-B[i];
}
#endif

#ifdef CUDA
void YplusisCtimesX(float *Y, float *X, float C, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yplusisctimesx<<<n_blocks,block_size>>>(Y,X,C,M);
}
#else
void YplusisCtimesX(float *Y, float *X, float C, int M) {
	for (int i=0; i<M; i++) Y[i] += C*X[i];
}
#endif

#ifdef CUDA
void Add(float *P, float *A, int M) {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	add<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Add(float *P, float *A, int M){
	for (int i=0; i<M; i++) P[i]+=A[i];
}
#endif

#ifdef CUDA
void Dubble(float *P, float *A, float norm,int M){
       int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	dubble<<<n_blocks,block_size>>>(P,A,norm,M);
}
#else
void Dubble(float *P, float *A, float norm,int M){
	for (int i=0; i<M; i++) P[i]*=norm/A[i];
}
#endif

#ifdef CUDA
void Boltzmann(float *P, float *A, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	boltzmann<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Boltzmann(float *P, float *A, int M){
	for (int i=0; i<M; i++) P[i]=exp(-A[i]);
}
#endif


#ifdef CUDA
void PutAlpha(float *g, float *phitot, float *phi_side, float chi, float phibulk, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	putalpha<<<n_blocks,block_size>>>(g,phitot,phi_side,chi,phibulk,M);
}
#else
void PutAlpha(float *g, float *phitot, float *phi_side, float chi, float phibulk, int M){
	for (int i=0; i<M; i++) g[i] = g[i] - chi*(phi_side[i]/phitot[i]-phibulk);
}
#endif

#ifdef CUDA
void Div(float *P, float *A, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	div<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Div(float *P, float *A, int M){
	for (int i=0; i<M; i++) if (A[i]>0) P[i]/=A[i]; else P[i]=0;
}
#endif

#ifdef CUDA
void AddG(float *g, float *phitot, float *alpha, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addg<<<n_blocks,block_size>>>(g,phitot,alpha,M);
}
#else
void AddG(float *g, float *phitot, float *alpha, int M){
	for (int i=0; i<M; i++) g[i]= g[i] -alpha[i] +1/phitot[i]-1;
}
#endif

#ifdef CUDA
void OneMinusPhitot(float *g, float *phitot, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	oneminusphitot<<<n_blocks,block_size>>>(g,phitot,M);
}
#else
void OneMinusPhitot(float *g, float *phitot, int M){
	for (int i=0; i<M; i++) g[i]= 1/phitot[i]-1;
}
#endif

#ifdef CUDA
void ComputeGN(float *GN, float *G, int M, int n_box){
	int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
	computegn<<<n_blocks,block_size>>>(GN,G,M,n_box);
}
#else
void ComputeGN(float* GN, float* G, int M, int n_box) {
	for (int p=0; p<n_box; p++) GN[p]=Sum(G+p*M,M); 
}
#endif

void H_Zero(float* H_P, int M){//this precedure should act on a host P.
	for (int i=0; i<M; i++) H_P[i] = 0;
}

void invert(float* KSAM, float *MASK, int M){ //only necessary in CPU
	for (int z=0; z<M; z++) if (MASK[z]==0) KSAM[z]=1.0; else KSAM[z]=0.0;
}

#ifdef CUDA
void SetBoundaries(float *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz){
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	bx<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#else
void SetBoundaries(float *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
	bx(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#endif


#ifdef CUDA
void RemoveBoundaries(float *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	b_x<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#else
void RemoveBoundaries(float *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
	b_x(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#endif
