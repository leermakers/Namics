

#ifdef CUDA

double *BlasResult;
cublasStatus_t stat;
cublasHandle_t handle;
int block_size=256;

__global__ void times(double *P, double *A, double *B, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void addtimes(double *P, double *A, double *B, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx]*B[idx];
}
__global__ void norm(double *P, double C, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] *= C;
}
__global__ void zero(double *P, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0.0;
}
__global__ void cp (double *P, double *A, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = A[idx];
}
__global__ void yisaminb(double *Y, double *A,double *B, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]-B[idx];
}
__global__ void yplusisctimesx(double *Y, double *X, double C, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*X[idx];
}
__global__ void add(double *P, double *A, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx];
}
__global__ void dubble(double *P, double *A, double norm, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]*=norm/A[idx];
}
__global__ void boltzmann(double *P, double *A, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=exp(-A[idx]);
}
__global__ void boltzmann(double *P, double *A, double *B, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=exp(-(A[idx]+B[idx]));
}
__global__ void putalpha(double *g,double *phitot,double *phi_side,double chi,double phibulk,int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) g[idx] = g[idx] - chi*(phi_side[idx]/phitot[idx]-phibulk);
}
__global__ void div(double *P,double *A,int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (A[idx]>0) P[idx] /=A[idx];
}
__global__ void oneminusphitot(double *g, double *phitot, int M){
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) g[idx]= 1/phitot[idx]-1;
}
__global__ void addg(double *g, double *phitot, double *alpha, int M) {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {
		g[idx]= g[idx] -alpha[idx] +1/phitot[idx]-1;
	}
}

__global__ void computegn(double *GN, double* G, int M, int n_box){
	int p= blockIdx.x*blockDim.x+threadIdx.x;
	if (p<n_box){
		cublasSasum(handle,M,G+p*M,1,GN[p]);
	}
}

void TransferDataToHost(double *H, double *D, int M) {
	cudaMemcpy(H, D, sizeof(double)*M,cudaMemcpyDeviceToHost);
}
void TransferDataToDevice(double *H, double *D, int M) {
	cudaMemcpy(D, H, sizeof(double)*M,cudaMemcpyHostToDevice);
}
void TransferIntDataToHost(int *H, int *D, int M) {
	cudaMemcpy(H, D, sizeof(int)*M,cudaMemcpyDeviceToHost);
}
void TransferIntDataToDevice(int *H, int *D, int M) {
	cudaMemcpy(D, H, sizeof(int)*M,cudaMemcpyHostToDevice);
}
#endif

#ifdef CUDA
__global__ void bx(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
	int idx, jx_mmx=jx*mmx, jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=P[bx1_jx+idx];
		P[jx_mmx+idx]=P[jx_bxm+idx];
	}
}
__global__ void b_x(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
	int idx, jx_mmx=jx*mmx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=0;
		P[jx_mmx+idx]=0;
	}
}
__global__ void by(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int idx, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=P[jy_by1+idx];
		P[jy_mmy+idx]=P[jy_bym+idx];
	}
}
__global__ void b_y(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int idx, jy_mmy=jy*mmy;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=0;
		P[jy_mmy+idx]=0;
	}
}
__global__ void bz(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=P[idx+bz1];
		P[idx+mmz]=P[idx+bzm];
	}
}
__global__ void b_z(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=0;
		P[idx+mmz]=0;
	}
}


#else
void bx(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
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
void b_x(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy){
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}
void by(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
	}
}
void b_y(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy){
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}
void bz(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
	}
}
void b_z(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy){
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
double Dot(double *x,double *y,int M){
	double result;
 	cublasSdot(handle,M,x,1,y,1,BlasResult);
	cudaMemcpy(&result,BlasResult,sizeof(double),cudaMemcpyDeviceToHost);
	return result;
}
#else
double Dot(double *x,double *y,int M){
	double result=0;
 	for (int i=0; i<M; i++) result +=x[i]*y[i];
	return result;
}
#endif

#ifdef CUDA
double Sum(double *x,int M){
	double result;
 	cublasSasum(handle,M,x,1,BlasResult);
	cudaMemcpy(&result,BlasResult,sizeof(double),cudaMemcpyDeviceToHost);
	return result;
}
#else
double Sum(double *x,int M){
	double result=0;
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
double *AllOnDev(int N) {
	double *X;
	if (cudaSuccess != cudaMalloc((void **) &X, sizeof(double)*N))
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
double *AllOnDev(int N) {
	double *X=NULL;
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
void AddTimes(double *P, double *A, double *B, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addtimes<<<n_blocks,block_size>>>(P,A,B,M);
}
#else
void AddTimes(double *P, double *A, double *B, int M){
	for (int i=0; i<M; i++) P[i]+=A[i]*B[i];
}
#endif

#ifdef CUDA
void Times(double *P, double *A, double *B, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	times<<<n_blocks,block_size>>>(P,A,B,M);
}
#else
void Times(double *P, double *A, double *B, int M){
	for (int i=0; i<M; i++) P[i]=A[i]*B[i];
}
#endif

#ifdef CUDA
void Norm(double *P, double C, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	norm<<<n_blocks,block_size>>>(P,C,M);
}
#else
void Norm(double *P, double C, int M){
	for (int i=0; i<M; i++) P[i] *= C;
}
#endif

#ifdef CUDA
void Zero(double* P, int M){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	zero<<<n_blocks,block_size>>>(P,M);
}
#else
void Zero(double* P, int M){
	for (int i=0; i<M; i++) P[i] =0;
}
#endif

#ifdef CUDA
void Cp(double *P,double *A, int M){int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	cp<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Cp(double *P,double *A, int M){
	for (int i=0; i<M; i++) P[i] = A[i];
}
#endif

#ifdef CUDA
void YisAminB(double *Y, double *A, double *B, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaminb<<<n_blocks,block_size>>>(Y,A,B,M);
}
#else
void YisAminB(double *Y, double *A, double *B, int M){
	for (int i=0; i<M; i++) Y[i] = A[i]-B[i];
}
#endif

#ifdef CUDA
void YplusisCtimesX(double *Y, double *X, double C, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yplusisctimesx<<<n_blocks,block_size>>>(Y,X,C,M);
}
#else
void YplusisCtimesX(double *Y, double *X, double C, int M) {
	for (int i=0; i<M; i++) Y[i] += C*X[i];
}
#endif

#ifdef CUDA
void Add(double *P, double *A, int M) {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	add<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Add(double *P, double *A, int M){
	for (int i=0; i<M; i++) P[i]+=A[i];
}
#endif

#ifdef CUDA
void Dubble(double *P, double *A, double norm,int M){
       int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	dubble<<<n_blocks,block_size>>>(P,A,norm,M);
}
#else
void Dubble(double *P, double *A, double norm,int M){
	for (int i=0; i<M; i++) P[i]*=norm/A[i];
}
#endif

#ifdef CUDA
void Boltzmann(double *P, double *A, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	boltzmann<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Boltzmann(double *P, double *A, int M){
	for (int i=0; i<M; i++) P[i]=exp(-A[i]);
}
#endif

#ifdef CUDA
void Boltzmann(double *P, double *A, double *B, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	boltzmann<<<n_blocks,block_size>>>(P,A,B,M);
}
#else
void Boltzmann(double *P, double *A, double *B, int M){
	for (int i=0; i<M; i++) P[i]=exp(-(A[i]+B[i]));
}
#endif


#ifdef CUDA
void PutAlpha(double *g, double *phitot, double *phi_side, double chi, double phibulk, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	putalpha<<<n_blocks,block_size>>>(g,phitot,phi_side,chi,phibulk,M);
}
#else
void PutAlpha(double *g, double *phitot, double *phi_side, double chi, double phibulk, int M){
	for (int i=0; i<M; i++) g[i] = g[i] - chi*(phi_side[i]/phitot[i]-phibulk);
}
#endif

#ifdef CUDA
void Div(double *P, double *A, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	div<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Div(double *P, double *A, int M){
	for (int i=0; i<M; i++) if (A[i]>0) P[i]/=A[i]; else P[i]=0;
}
#endif

#ifdef CUDA
void AddG(double *g, double *phitot, double *alpha, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addg<<<n_blocks,block_size>>>(g,phitot,alpha,M);
}
#else
void AddG(double *g, double *phitot, double *alpha, int M){
	for (int i=0; i<M; i++) g[i]= -g[i] +alpha[i] +1/phitot[i]+phitot[i];
}
#endif

#ifdef CUDA
void OneMinusPhitot(double *g, double *phitot, int M){
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	oneminusphitot<<<n_blocks,block_size>>>(g,phitot,M);
}
#else
void OneMinusPhitot(double *g, double *phitot, int M){
	for (int i=0; i<M; i++) g[i]= 1/phitot[i]-1;
}
#endif

#ifdef CUDA
void ComputeGN(double *GN, double *G, int M, int n_box){
	int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
	computegn<<<n_blocks,block_size>>>(GN,G,M,n_box);
}
#else
void ComputeGN(double* GN, double* G, int M, int n_box) {
	for (int p=0; p<n_box; p++) GN[p]=Sum(G+p*M,M); 
}
#endif

void H_Zero(double* H_P, int M){//this precedure should act on a host P.
	for (int i=0; i<M; i++) H_P[i] = 0;
}

void invert(double* KSAM, double *MASK, int M){ //only necessary in CPU
	for (int z=0; z<M; z++) if (MASK[z]==0) KSAM[z]=1.0; else KSAM[z]=0.0;
}

#ifdef CUDA
void SetBoundaries(double *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz){
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	bx<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#else
void SetBoundaries(double *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
	bx(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#endif


#ifdef CUDA
void RemoveBoundaries(double *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	b_x<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#else
void RemoveBoundaries(double *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
	b_x(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#endif



