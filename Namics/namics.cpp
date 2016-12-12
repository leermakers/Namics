#include "namics.h"
#include "tools.cpp" 
#include "input.cpp"
#include "segment.cpp"
#include "molecule.cpp"
#include "output.cpp"
#include "lattice.cpp"
#include "system.cpp"
#include "engine.cpp"
#include "newton.cpp"

#ifdef CUDA
__global__ void distributeg1(float *G1, float *g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
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

__global__ void collectphi(float* phi, float* GN,float* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
        int i = blockIdx.x*blockDim.x+threadIdx.x;
        int j = blockIdx.y*blockDim.y+threadIdx.y;
        int k = blockIdx.z*blockDim.z+threadIdx.z;
        int pos_r=(JX+JY+1);
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
#endif

void subdomain(int *fBx,int *fBy,int *fBz, int *fPx,int *fPy,int *fPz, int zloc){
    int loop = 0;
    for (int xcord=1; xcord<=10; xcord++){
		for(int ycord=1; ycord<=10; ycord++){	
		fPx[loop] = (xcord)*5-2 ; fPy[loop] = (ycord)*5-2 ; fPz[loop] = zloc ;
		fBx[loop] = (xcord*5-2)-10; fBy[loop] = (ycord*5-2)-10 ; fBz[loop] = zloc-10 ;	
		loop = loop+1;
		}
	}
	for (int diag=1; diag<=12; diag++){
		fPx[loop] = diag*4 ; fPy[loop] = diag*4 ; fPz[loop] = 26 ;
		fBx[loop] = (diag*4)-10; fBy[loop] = (diag*4)-10 ; fBz[loop] = 16 ;
		loop = loop+1;
	}
}

void Side(float *X_side, float *X, int M) {
	Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);

	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Add(X_side+JX,X+JY,MM-JX-JY); Add(X_side+JY,X+JX,MM-JY-JX);
	Add(X_side+1,X+JY,MM-JY-1); Add(X_side+JY,X+1,MM-JY-1);
	Add(X_side+1,X+JX,MM-JX-1); Add(X_side+JX,X+1,MM-JX-1);
	Norm(X_side,1.0/12.0,M);
}



void SideCUBIC(float *X_side, float *X, int M) {
	Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Norm(X_side,1.0/6.0,M);
}

void advanced_average(float *X_side, float *X, int M){
        Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);

	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);

	Add(X_side+JX+JY,X,M-JX-JY); Add(X_side, X+JX+JY, M-JX-JY);
	Add(X_side+JY+1,X,M-JY-1); Add(X_side, X+JY+1, M-JY-1);
	Add(X_side+JX+1,X,M-JX-1); Add(X_side, X+JX+1, M-JX-1);
        Add(X_side+JX-JY,X,M-JX+JY); Add(X_side, X+JX-JY, M-JX+JY);
	Add(X_side+JY-1,X,M-JY+1); Add(X_side, X+JY-1, M-JY+1);
	Add(X_side+JX-1,X,M-JX+1); Add(X_side, X+JX-1, M-JX+1);

	Add(X_side+JX+JY+1,X,M-JX-JY-1); Add(X_side, X+JX+JY+1, M-JX-JY-1);
	Add(X_side+JX+JY-1,X,M-JX-JY+1); Add(X_side, X+JX+JY-1, M-JX-JY+1);
	Add(X_side+JX-JY+1,X,M-JX+JY-1); Add(X_side, X+JX-JY+1, M-JX+JY-1);
	Add(X_side-JX+JY+1,X,M+JX-JY-1); Add(X_side, X-JX+JY+1, M+JX-JY-1);

	Norm(X_side,1.0/26.0,M);
}

void Propagate(float* G, float* G1, int s_from, int s_to) { //on small boxes
	int MMM=M*n_box;
	float *gs = G+MMM*s_to, *gs_1 = G+MMM*s_from, *g = G1;
	Zero(gs,MMM);
	for (int p=0; p<n_box; p++) SetBoundaries(gs_1+M*p,jx,jy,bx1,bxm,by1,bym,bz1,bzm,Mx,My,Mz);
	
	Add(gs+jx,gs_1,MMM-jx); Add(gs,gs_1+jx,MMM-jx);
	Add(gs+jy,gs_1,MMM-jy); Add(gs,gs_1+jy,MMM-jy);
	Add(gs+1,gs_1,MMM-1);  Add(gs,gs_1+1, MMM-1);	
	Add(gs+jx,gs_1+jy,MMM-jx-jy); Add(gs+jy,gs_1+jx,MMM-jx-jy);
	Add(gs+1,gs_1+jy,MMM-jy-1); Add(gs+jy,gs_1+1,MMM-jy-1);
	Add(gs+1,gs_1+jx,MMM-jx-1); Add(gs+jx,gs_1+1,MMM-jx-1);
	Norm(gs,1.0/12.0,MMM); Times(gs,gs,g,MMM);
}

void PROPAGATE(float *G, float *G1, int s_from, int s_to) { //on big box
	float *gs = G+MM*(s_to), *gs_1 = G+MM*(s_from), *g = G1;
	Zero(gs,MM);
	SetBoundaries(gs_1,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(gs+JX,gs_1,MM-JX); Add(gs,gs_1+JX,MM-JX);
	Add(gs+JY,gs_1,MM-JY); Add(gs,gs_1+JY,MM-JY);
	Add(gs+1,gs_1,MM-1);  Add(gs,gs_1+1, MM-1);
	Add(gs+JX,gs_1+JY,MM-JX-JY); Add(gs+JY,gs_1+JX,MM-JX-JY);
	Add(gs+1,gs_1+JY,MM-JY-1); Add(gs+JY,gs_1+1,MM-JY-1);
	Add(gs+1,gs_1+JX,MM-JX-1); Add(gs+JX,gs_1+1,MM-JX-1);
	Norm(gs,1.0/12.0,MM); Times(gs,gs,g,MM);
}



void PropagateCubic(float* G, float* G1, int s_from, int s_to) { //on small boxes
	int MMM=M*n_box;
	float *gs = G+MMM*s_to, *gs_1 = G+MMM*s_from, *g = G1;
	Zero(gs,MMM);
	for (int p=0; p<n_box; p++) SetBoundaries(gs_1+M*p,jx,jy,bx1,bxm,by1,bym,bz1,bzm,Mx,My,Mz);
	Add(gs+jx,gs_1,MMM-jx); Add(gs,gs_1+jx,MMM-jx);
	Add(gs+jy,gs_1,MMM-jy); Add(gs,gs_1+jy,MMM-jy);
	Add(gs+1,gs_1,MMM-1);  Add(gs,gs_1+1, MMM-1);
	Norm(gs,1.0/6.0,MMM); Times(gs,gs,g,MMM);
}

void PROPAGATECUBIC(float *G, float *G1, int s_from, int s_to) { //on big box
	float *gs = G+MM*(s_to), *gs_1 = G+MM*(s_from), *g = G1;
	Zero(gs,MM);
	SetBoundaries(gs_1,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(gs+JX,gs_1,MM-JX); Add(gs,gs_1+JX,MM-JX);
	Add(gs+JY,gs_1,MM-JY); Add(gs,gs_1+JY,MM-JY);
	Add(gs+1,gs_1,MM-1);  Add(gs,gs_1+1, MM-1);
	Norm(gs,1.0/6.0,MM); Times(gs,gs,g,MM);
}


#ifdef CUDA
void DistributeG1(float *G1, float *g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	//int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
	//distributeg1<<<n_blocks,block_size>>>(G1,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
        dim3 blocks(ceil(Mx/8.0),ceil(My/8.0),ceil(Mz/8.0));
        dim3 blockdim(8,8,8);
        distributeg1<<<blocks,blockdim>>>(G1,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);

}
#else
void DistributeG1(float *G1, float *g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int pos_l=-M;
	int pos_x,pos_y,pos_z;
	int Bxp,Byp,Bzp;
	int ii=0,jj=0,kk=0;

	for (int p=0; p<n_box; p++) { pos_l +=M; ii=0; Bxp=H_Bx[p]; Byp=H_By[p]; Bzp=H_Bz[p];
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
void CollectPhi(float* phi, float* GN, float* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
 	 dim3 blocks(ceil(Mx/8.0),ceil(My/8.0),ceil(Mz/8.0));
        dim3 blockdim(8,8,8);
        collectphi<<<blocks,blockdim>>>(phi,GN,rho,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
}
#else
void CollectPhi(float* phi, float* GN, float* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int pos_l=-M;
	int pos_x,pos_y,pos_z;
	int Bxp,Byp,Bzp;
	float Inv_H_GNp;
	int ii=0,jj=0,kk=0;

	for (int p=0; p<n_box; p++) { pos_l +=M; ii=0; Bxp=Bx[p]; Byp=By[p]; Bzp=Bz[p]; Inv_H_GNp=1.0/GN[p];
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


//void ComputePhis(){ //monomers are not connected to the subboxes and are distributed over the overall box (that's why MM is used)
//	if (charges) {
//		Cp(alpha,u,MM); Boltzmann(phi,alpha,MM);
//		Add(alpha,psi,MM); Boltzmann(phi_na,alpha,MM); Norm(phi_na,phib_s,MM);
//		YisAminB(alpha,u,psi,MM); Boltzmann(phi_cl,alpha,MM); Norm(phi_cl,phib_s,MM);
//		Cp(alpha,u,MM); YplusisCtimesX(alpha,psi,alpha_seg,MM);
//		Boltzmann(G1,alpha,MM);
//		Times(G1,G1,KSAM,MM);
//		Times(phi_cl,phi_cl,KSAM,MM);
//		Times(phi_na,phi_na,KSAM,MM);
//		Times(phi,phi,KSAM,MM);
//	} else {
//		Boltzmann(G1,u,MM); SetBoundaries(G1,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
//		Times(G1,G1,KSAM,MM);//monomers can not enter the sites taken by nodes/particles.
//		Cp(phi,G1,MM);
//	}
//}

void ComputePhi(){ //compute all densities.
//both block are equally long! 
	Boltzmann(G1,u,2*MM); 
	SetBoundaries(G1,   JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	SetBoundaries(G1+MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Times(G1,G1,KSAM,MM);Times(G1+MM,G1+MM,KSAM,MM);
	Cp(GG_F,G1,MM);
 
	for(int s=1; s<N_A; s++) PROPAGATE(GG_F,G1,s-1,s);
	RemoveBoundaries(GG_F+(N_A-1)*MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); 
	GNA=Sum(GG_F+(N_A-1)*MM,MM); phib[0]=Theta_A/GNA; phib[1]=1.0-phib[0];
	for(int s=0; s<N_A/2; s++) AddTimes(phi,GG_F+s*MM,GG_F+(N_A-(s+1))*MM,MM); 
	if (N_A>1) Norm(phi,2.0,MM);
	if ((N_A%2)==1) AddTimes(phi,GG_F+(N_A/2)*MM,GG_F+(N_A/2)*MM,MM);
	Div(phi,G1,MM); 
	Norm(phi,Theta_A/GNA/N_A,MM);

	Cp(GG_F,G1+MM,MM);
	for(int s=1; s<N_B; s++) PROPAGATE(GG_F,G1+MM,s-1,s);
	RemoveBoundaries(GG_F+(N_B-1)*MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); 
	for(int s=0; s<N_B/2; s++) AddTimes(phi+MM,GG_F+s*MM,GG_F+(N_B-(s+1))*MM,MM); 
	if (N_B>1) Norm(phi+MM,2.0,MM);
	if ((N_B%2)==1) AddTimes(phi+MM,GG_F+(N_B/2)*MM,GG_F+(N_B/2)*MM,MM);
	Div(phi+MM,G1+MM,MM); Norm(phi+MM,phib[1]/N_B,MM);

	Cp(Gg_f,mask,M*n_box);

	Zero(rho,M*n_box);
	DistributeG1(G1,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
	
	for (int s=1; s<=N; s++) Propagate(Gg_f,g1,s-1,s);
	Cp(Gg_b+(N%2)*M*n_box,g1,M*n_box);
       AddTimes(rho,Gg_f+N*M*n_box,Gg_b+(N%2)*M*n_box,M*n_box);

	for (int s=N-1; s>0; s--) {
		Propagate(Gg_b,g1,(s+1)%2,s%2);
		AddTimes(rho,Gg_f+s*M*n_box,Gg_b+(s%2)*M*n_box,M*n_box);
 	}
	RemoveBoundaries(Gg_f+N*n_box*M,jx,jy,bx1,bxm,by1,bym,bz1,bzm,Mx,My,Mz);
	ComputeGN(GN_A,Gg_f+N*n_box*M,M,n_box);
	CollectPhi(phi+2*MM,GN_A,rho,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);

	Zero(rho,M*n_box);
	DistributeG1(G1+MM,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);

	for (int s=1; s<=N; s++) Propagate(Gg_f,g1,s-1,s);
	Cp(Gg_b+(N%2)*M*n_box,g1,M*n_box);
       AddTimes(rho,Gg_f+N*M*n_box,Gg_b+(N%2)*M*n_box,M*n_box);

	for (int s=N-1; s>0; s--) {
		Propagate(Gg_b,g1,(s+1)%2,s%2);
		AddTimes(rho,Gg_f+s*M*n_box,Gg_b+(s%2)*M*n_box,M*n_box);
 	}
	RemoveBoundaries(Gg_f+N*n_box*M,jx,jy,bx1,bxm,by1,bym,bz1,bzm,Mx,My,Mz);
	ComputeGN(GN_B,Gg_f+N*n_box*M,M,n_box);
	CollectPhi(phi+3*MM,GN_B,rho,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
	/*RemoveBoundaries(phi,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); 
		RemoveBoundaries(phi+MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
			RemoveBoundaries(phi+2*MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
				RemoveBoundaries(phi+3*MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);*/

	Div(phi+2*MM,G1,2*MM);
		/*cout << Sum(phi+2*MM,MM)/N << " and " << Sum(phi+3*MM,MM)/N << endl; 
		cout << Sum(phi+0*MM,MM)/N_A << " and " << Sum(phi+1*MM,MM)/N_B << endl; */
	Cp(PHI,phi,2*MM); Add(PHI,phi+2*MM,2*MM);
}


void ComputeG(){
//	ComputePhis();
#ifdef CUDA
	Zero(phi,4*MM);
#endif
	H_Zero(H_phi,4*MM); //initialize phi for polymer;
	ComputePhi(); //run this for calculation the chain densities between nodes.
	Cp(phitot,phi,MM); Add(phitot,phi+MM,MM); Add(phitot,phi+2*MM,MM);Add(phitot,phi+3*MM,MM);Add(phitot,MASK,MM);
	Side(phi_side,PHI,MM);Side(phi_side+MM,PHI+MM,MM);
	Cp(g,u,iv);
	PutAlpha(g,phitot,phi_side+MM,CHI,phib[1],MM);
	PutAlpha(g+MM,phitot,phi_side,CHI,phib[0],MM);
	Cp(alpha,g,MM); Add(alpha,g+MM,MM); Norm(alpha,1.0/n_seg,MM);
	for (int i=0; i<n_seg; i++) {
		AddG(g+i*MM,phitot,alpha,MM);
		RemoveBoundaries(g+i*MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	}
//	if (charges) {
//		Cp(psi_0,psi,MM);
//		Add(phitot,phi_na,MM);
//		Add(phitot,phi_cl,MM);
//		YisAminB(q,phi_na,phi_cl,MM);
//		YplusisCtimesX(q,phi+MM,alpha_seg,MM);
//		Norm(q,factor/6,MM);
//		Side(psi_side,psi,MM);
//		Cp(psi,psi_side,MM); Add(psi,q,MM);
//		YisAminB(g+MM,psi_0,psi,MM);
//		RemoveBoundaries(g+MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
//	}

}

void Ax(float* A, float* X, int N){//From Ax_B; below B is not used: it is assumed to contain a row of unities.
	float* U = new float[N*N];
	float* S = new float[N];
	float* VT = new float[N*N];
	integer MM = (integer)N, NN = (integer)N;
	integer LDA=MM, LDU=MM, LDVT=NN, INFO, LWORK;
	int lwork;
	float WKOPT;
	float* WORK;
	char JOBU='S'; //'S' is nodig om alleen de eerste N colommen in U te schrijven.
	char JOBVT='A';

	LWORK = -1; //grootte hulpgeheugen aanvragen
	sgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT, &WKOPT, &LWORK, &INFO );
	lwork = (int)WKOPT;
	WORK = (float*)malloc( lwork*sizeof(float) );
	LWORK = (integer)lwork; //nu uitrekenen.
	sgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT,WORK, &LWORK, &INFO );
	if (INFO >0) { //error message genereren
	};
	free(WORK);
	for (int i=0; i<N; i++) X[i]=0;
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += U[i*N + j];//*B[j];
	for (int i=0; i<N; i++) {S[i] = X[i]/S[i]; X[i]=0;} //S is use decause it is no longer needed.
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += VT[i*N + j]*S[j];
	delete U;
	delete S;
	delete VT;
}
void DIIS(float* x, float* x_x0, float* xR, float* Aij, float* Apij,float* Ci, int k, int m, int iv) {
	float normC=0; int posi;
	if (k_diis>m) { k_diis =m;
		for (int i=1; i<m; i++) for (int j=1; j<m; j++)
		Aij[m*(i-1)+j-1]=Aij[m*i+j]; //remove oldest elements
	}
	for (int i=0; i<k_diis; i++) {posi = k-k_diis+1+i; if (posi<0) posi +=m;
		Aij[i+m*(k_diis-1)] = Aij[k_diis-1+m*i] = Dot(x_x0+posi*iv, x_x0+k*iv,iv);	}
		// write to (compressed) matrix Apij
	for (int i=0; i<k_diis; i++) for (int j=0; j<k_diis; j++) {
		Apij[j+k_diis*i] = Aij[j+m*i];
	}
	Ax(Apij,Ci,k_diis);
	for (int i=0; i<k_diis; i++) normC +=Ci[i];
	for (int i=0; i<k_diis; i++) {Ci[i] =Ci[i]/normC; }
	Zero(x,iv);
	posi = k-k_diis+1; if (posi<0) posi +=m;

	YplusisCtimesX(x,xR+posi*iv,Ci[0],iv); //pv = Ci[0]*xR[0];
	for (int i=1; i<k_diis; i++) {
		posi = k-k_diis+1+i; if (posi<0) posi +=m;
		YplusisCtimesX(x,xR+posi*iv,Ci[i],iv);
	}
}
float GrandPotential(){
	float GP=0;
	float volume = 1.0*MX*MY*MZ;
	RemoveBoundaries(phi,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	RemoveBoundaries(phi+MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	RemoveBoundaries(phi+2*MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	RemoveBoundaries(phi+3*MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	GP -= (Sum(phi,MM)-volume*phib[0])/N_A;
	GP -= (Sum(phi+MM,MM)-volume*phib[1])/N_B;
	GP -= Sum(phi+2*MM,MM)/N;
	GP -= Sum(phi+3*MM,MM)/N;
	GP -= Dot(alpha,MASK,MM);
		
	return GP;	
}
float Helmholtz(){
	float F_Helmholtz, entropy, energy;
	entropy = 0;
	energy = 0;

	RemoveBoundaries(u,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	RemoveBoundaries(u+MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	RemoveBoundaries(PHI,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	RemoveBoundaries(PHI+MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	entropy = Dot(PHI,u,MM); 
	

	//SetBoundaries(PHI,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	SetBoundaries(PHI+MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Side(phi_side,PHI,MM);
	energy = CHI*Dot(PHI+MM,phi_side,MM);

#ifdef CUDA
TransferDataToHost(H_GN_A,GN_A,n_box);
TransferDataToHost(H_GN_B,GN_B,n_box);
#endif
	for (int p=0; p<n_box; p++) {entropy +=log(H_GN_A[p]/N); entropy +=log(H_GN_B[p]/N);}
	RemoveBoundaries(phi,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	RemoveBoundaries(phi+MM,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	entropy -=Sum(phi,MM)/N_A*log(phib[0]);
	entropy -=Sum(phi+MM,MM)/N_B*log(phib[1]);
	F_Helmholtz = energy - entropy;
	return F_Helmholtz;
}

float SCF() {
	for (int i=0; i<MX+2; i++) for (int j=0; j<MY+2; j++) for (int k=0; k<MZ+2; k++) {
		if (k<MZ/2) {
			H_u[JX*i+JY*j+k]= -0.38335;
			H_u[JX*i+JY*j+k+MM]= 0.38335;
		} else {
			H_u[JX*i+JY*j+k]=0;
			H_u[JX*i+JY*j+k+MM]=0;
		}
	}

#ifdef CUDA
	TransferDataToDevice(H_mask, mask, M*n_box);
	TransferDataToDevice(H_MASK, MASK, MM);
	TransferDataToDevice(H_KSAM, KSAM, MM);
	TransferDataToDevice(H_u, u, MM*n_seg);
	TransferIntDataToDevice(H_Bx, Bx, n_box);
	TransferIntDataToDevice(H_By, By, n_box);
	TransferIntDataToDevice(H_Bz, Bz, n_box);
	TransferIntDataToDevice(H_Px, Px, n_box);
	TransferIntDataToDevice(H_Py, Py, n_box);
	TransferIntDataToDevice(H_Pz, Pz, n_box);
//	if (charges) TransferDataToDevice(H_psi,psi,MM);
#else
	Cp(u,H_u,MM*n_seg);
//	if (charges) Cp(psi,H_psi,MM);
#endif
	Zero(x0,iv);
	it=0; k_diis=1; k=0;
	ComputeG();
	YplusisCtimesX(x,g,-eta,iv);
	YisAminB(x_x0,x,x0,iv);
	Cp(xR,x,iv);
	error = sqrt(Dot(g,g,iv));
	printf("DIIS has been notified\n");
	printf("Your guess = %1e \n",error);
	while (error > tolerance && it < iterations) {
		it++;
		Cp(x0,x,iv); ComputeG();
		k=it % m; k_diis++; //plek voor laatste opslag
		YplusisCtimesX(x,g,-eta,iv);
		Cp(xR+k*iv,x,iv); YisAminB(x_x0+k*iv,x,x0,iv);
		DIIS(x,x_x0,xR,Aij,Apij,Ci,k,m,iv);
		error = sqrt(Dot(g,g,iv));
		if(it%10 == 0){
			printf("it = %i error = %1e \n",it,error);
		}
	}
	return Helmholtz();
}



void vtk_output(string filename, float* X) {
	FILE * fp;
	fp = fopen(filename.c_str(),"w+");
	fprintf(fp, "# vtk DataFile Version 3.0 \nvtk output \nASCII \nDATASET STRUCTURED_POINTS \nDIMENSIONS %i %i %i \n", MX, MY, MX);
	fprintf(fp,"SPACING 1 1 1 \nORIGIN 0 0 0 \nPOINT_DATA %i \n", MX*MY*MZ);
	fprintf(fp,"SCALARS Box_profile float\nLOOKUP_TABLE default \n");

	for (int i=1; i<MX+1; i++) for (int j=1; j<MY+1; j++) for (int k=1; k<MZ+1; k++)
	fprintf(fp," %f \n", X[i*JX+j*JY+k]);
	fclose(fp);
}


int main(int argc, char *argv[]) {
	string fname;
	string filename;
	string line;
	ofstream out_file;


	if (argc == 2) fname = argv[1]; else {printf("Use: namics filename -without extension- \n"); return 1;}
	filename = fname + ".in";
	vector<Input*> In; In.push_back(new Input(filename)); if (In[0]->Input_error) {return 0;}
	vector<Lattice*> Lat; Lat.push_back(new Lattice(In,In[0]->LatList[0])); if (!Lat[0]->CheckInput()) {return 0;}
	
	n_seg=In[0]->MonList.size();

	vector<Segment*> Seg; 
	for (int i=0; i<n_seg; i++) Seg.push_back(new Segment(In,In[0]->MonList[i],i,n_seg,Lat[0]->MX,Lat[0]->MY,Lat[0]->MZ));
	for (int i=0; i<n_seg; i++) {  
		for (int k=0; k<n_seg; k++) Seg[i]->PutChiKEY(Seg[k]->name); if (!Seg[i]->CheckInput()) return 0;   
	}

	n_mol = In[0]->MolList.size(); 
 
	vector<Molecule*> Mol; for (int i=0; i<n_mol; i++) Mol.push_back(new Molecule(In,Seg,In[0]->MolList[i],i,n_mol));
	for (int i=0; i<n_mol; i++) {if (!Mol[i]->CheckInput()) return 0;}
	vector<System*> Sys; Sys.push_back(new System(In,Lat,Seg,Mol,In[0]->SysList[0])); if (!Sys[0]->CheckInput()) {return 0;}
	if (!Sys[0]->CheckChi_values(n_seg)) return 0;
	if (In[0]->NewtonList.size()>0) {
		vector<Newton*> New; New.push_back(new Newton(In,Sys,In[0]->NewtonList[0])); if (!New[0]->CheckInput()) {return 0;}
	}
	if (In[0]->EngineList.size()>0) {
		vector<Engine*> Eng; Eng.push_back(new Engine(In,Sys,In[0]->EngineList[0])); if (!Eng[0]->CheckInput()) {return 0;}
	}
	n_out = In[0]->OutputList.size(); 
	vector<Output*> Out; for (int i=0; i<n_out; i++) Out.push_back(new Output(In,In[0]->OutputList[i],i,n_out));
	for (int i=0; i<n_out; i++) {
		if (!Out[i]->CheckOutInput()) {return 0;}
		string template_ = Out[i]->GetValue("template");
		if (!Out[i]->Load(template_)) return 0;  
	}
	cal_types.push_back("micro-emulsion"); 
	cal_types.push_back("membrane"); 

	if (!In[0]->Get_string(Sys[0]->GetValue("calculation_type"),cal_type,cal_types,"In 'sys' the parameter 'calculation_type' is required")) {return 0;} 
bool MEmulsion=cal_type=="micro-emulsion"; if (MEmulsion) cout << "micro-emulsion" << endl; 
bool Membrane=cal_type=="membrane"; if (Membrane) cout<< "micro-emulsion" << endl;

	MX=Lat[0]->MX; MY=Lat[0]->MY;MZ=Lat[0]->MZ;  


	Aij = new float[m*m]; for (int i=0; i<m*m; i++) Aij[i]=0; //needed for SVD which is done on cpu.
	Ci = new float[m]; for (int i=0; i<m; i++) Ci[i]=0;
	Apij = new float[m*m]; for (int i=0; i<m*m; i++) Apij[i]=0;

#ifdef CUDA
	cudaDeviceReset();
	stat = cublasCreate(&handle); if (stat !=CUBLAS_STATUS_SUCCESS) {printf("CUBLAS failed \n");}
	cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
#endif
	n_seg=2; //number of segment types.




	if (MEmulsion) { //Read info from 'in' file for micro-emulsion model. 

	n_cosol = 0;
	n_sol = 2; 
	N = 20; 
	N_A = 4;
	N_B = 4;
	CHI = 0.6;

		charges = false;

	} else { //Read info forom 'in' file for 'membrane' model
		charges = false;
	}


	JX = (MY+2)*(MZ+2); JY = MZ+2; MM=(MX+2)*(MY+2)*(MZ+2);
	BX1=MX; BXM=1; BY1=MY; BYM=1; BZ1=1; BZM=MZ; // Large BOX periodic in X and Y, but reflecting in Z.
	iv = MM*n_seg;
	Theta_A = 1.0*MX*MY*MZ/2.0;	//Fill half the Big box with solvent A.

		phib = new float[2];
		H_MASK = new float[MM]; H_Zero(H_MASK,MM);
		H_KSAM = new float[MM];
		H_u = new float[MM*n_seg];
		H_G1= new float[MM*n_seg]; //there are just two segment types A and B
		H_phi = new float[MM*4];  //we will collect phi(A4), phi(B4), phi(A20) and phi(B20) in phi vector.
		H_PHI = new float[MM*n_seg]; //will collect phi(A), phi(B).
		//H_psi = new float[MM];

#ifdef CUDA
		BlasResult= (float*)AllOnDev(1);
		x = (float*)AllOnDev(iv); u=x; Zero(x,iv);
		x0 = (float*)AllOnDev(iv);
		g= (float*)AllOnDev(iv);
		xR= (float*)AllOnDev(m*iv);
		x_x0= (float*)AllOnDev(m*iv);
		phi= (float*)AllOnDev(MM*4);	   //NA, NB, block A and Bloc B of copolymer.
		PHI = (float*)AllOnDev(MM*n_seg); //only A and B types
		phi_side = (float*)AllOnDev(MM*n_seg);
		phitot= (float*)AllOnDev(MM);
		G1= (float*)AllOnDev(MM*n_seg);
		alpha= (float*)AllOnDev(MM);
		MASK=(float*)AllOnDev(MM); Zero(MASK,MM);
		KSAM =(float*)AllOnDev(MM);
//		if (charges) {
//			phi_na =(float*)AllOnDev(MM);
//			phi_cl =(float*)AllOnDev(MM);
//			psi_0 =(float*)AllOnDev(MM);
//			psi_side =(float*)AllOnDev(MM);
//			psi=x+MM;
//			q =(float*)AllOnDev(MM);
//		}

#else
		x = new float[iv]; u=x;  Zero(x,iv);
		x0 = new float[iv];
		g = new float[iv];
		xR = new float[m*iv];
		x_x0 =new float[m*iv];
		phi = H_phi;
		PHI = H_PHI;	  //collects densities per segment type.
		phi_side = new float[MM*n_seg];
		phitot = new float[MM];
		G1 = H_G1;
		alpha = new float[MM];
		MASK=H_MASK;
		KSAM=H_KSAM;
//		if (charges) {
//			psi=x+MM;
//			phi_na = new float[MM];
//			phi_cl = new float[MM];
//			q = new float[MM];
//			psi_0 = new float[MM];
//			psi_side = new float[MM];
//		}
#endif
//		if (!Full) for (int i=1; i<MX+1; i++) for (int j=1; j<MY+1; j++) for (int k=1; k<MZ+1; k++)
//		in_file >> H_u[JX*i+JY*j+k];//only read the potentials when needed. In Full SCF compute u's from scratch!
//		in_file.close();
//	} else cout << "file " << filename << " is not available." << endl;

//	if (charges && !Full) {
//		filename = fname + "_psi.vtk";
//		in_file.open(filename.c_str());
//		if (in_file.is_open()) {
//			getline(in_file,line_); getline(in_file,line_);
//			getline(in_file,line_); getline(in_file,line_);
//			getline(in_file,line_); getline(in_file,line_);
//			getline(in_file,line_); getline(in_file,line_);
//			getline(in_file,line_); getline(in_file,line_);
//			if (!Full) for (int i=1; i<MX+1; i++) for (int j=1; j<MY+1; j++) for (int k=1; k<MZ+1; k++)
//			in_file >> H_psi[JX*i+JY*j+k];//only read the potentials when needed. In Full SCF compute u's from scratch!
//
//		} else cout << "file " << filename << " is not available." << endl;
//	}

//	int n_lines=0;
//	filename = fname + ".sub";  //sub-boxes
//	in_file.open(filename.c_str());
//	if (in_file.is_open()) {
//		in_file >> line_ >>  N >> line_ >> Mx >> My >> Mz;
//		while (getline(in_file,line_)) {n_lines++;
//		}
//		in_file.close();
//	} else cout << "file " << filename << " is not available." << endl;
//	if (n_lines%4 >0) cout << ".sub file does not contain a multiple of 4 lines" << endl;
//	n_box = n_lines/4;



	n_box = 112;
	Mx=My=Mz=20; 

	jx = (My+2)*(Mz+2); jy = Mz+2; M=(Mx+2)*(My+2)*(Mz+2);
	bx1=1; bxm=Mx; by1=1; bym=My; bz1=1; bzm=Mz; //reflecting for small boxes

	H_mask = new float[M*n_box];
#ifdef CUDA
	mask= (float*)AllOnDev(M*n_box);
#else
	mask=H_mask;
#endif

	H_Zero(H_mask,M*n_box);
	H_Bx  = new int[n_box];  H_By  = new int[n_box];   H_Bz  = new int[n_box];
	H_Px = new int[n_box];  H_Py = new int[n_box];   H_Pz = new int[n_box];
#ifdef CUDA
	Bx=(int*)AllIntOnDev(n_box);
	By=(int*)AllIntOnDev(n_box);
	Bz=(int*)AllIntOnDev(n_box);
	Px=(int*)AllIntOnDev(n_box);
	Py=(int*)AllIntOnDev(n_box);
	Pz=(int*)AllIntOnDev(n_box);

#else
	Bx=H_Bx; By=H_By; Bz=H_Bz;
	Px=H_Px; Py=H_Py; Pz=H_Pz;

#endif

//	filename = fname + ".sub";  //sub-boxes
//	in_file.open(filename.c_str());
//	if (in_file.is_open()) {
//		for (int p=0; p<n_box; p++){
//			getline(in_file,line_);
//			in_file >> H_Bx[p] >> H_By[p] >> H_Bz[p];
//			in_file >> H_Px1[p] >> H_Py1[p] >> H_Pz1[p];
//			in_file >> H_Px2[p] >> H_Py2[p] >> H_Pz2[p];
//			getline(in_file,line_);
//
//			if (H_Bx[p]<1) {H_Bx[p] +=MX; H_Px1[p] +=MX; H_Px2[p] +=MX;}
//			if (H_By[p]<1) {H_By[p] +=MY; H_Py1[p] +=MY; H_Py2[p] +=MY;}
//			if (H_Bz[p]<1) {H_Bz[p] +=MZ; H_Pz1[p] +=MZ; H_Pz2[p] +=MZ;}
//		}
//		in_file.close();
//	} else cout << "file " << filename << " is not available." << endl;
	float Free_energy;
	int variation = 1;
	
	if (variation == 1)
	{
		
		ofstream varfile;
		varfile.open ("variation.dat");
		
		for (int var=0; var<6; var++){
		printf("\n");
		printf("Problem %d of 10 \n",var+1);
		printf("\n");
		int interface = 26+var;
		//Theta_A = Theta_A + MZ*MZ/10.0;
		subdomain(H_Bx,H_By,H_Bz,H_Px,H_Py,H_Pz, interface);
		for (int p=0; p<n_box; p++){
			H_mask[p*M + jx*(H_Px[p]-H_Bx[p])+jy*(H_Py[p]-H_By[p])+(H_Pz[p]-H_Bz[p])]=1;
			H_MASK[((H_Px[p]-1)%MX+1)*JX + ((H_Py[p]-1)%MY+1)*JY + (H_Pz[p]-1)%MZ+1]=1;
		}

		invert(H_KSAM,H_MASK,MM);
		H_g1= new float[M*n_box]; H_Zero(H_g1,M*n_box);
		H_rho = new float[M*n_box];
		H_GN_A = new float[n_box];
		H_GN_B = new float[n_box];

	#ifdef CUDA
		Gg_f=(float*)AllOnDev(M*N*n_box);
		GG_F=(float*)AllOnDev(MM*N_A); //Let N_A be the largest N of the solvents!!!!
		Gg_b=(float*)AllOnDev(M*2*n_box);
		GN_A=(float*)AllOnDev(n_box);
		GN_B=(float*)AllOnDev(n_box);
		rho =(float*)AllOnDev(M*n_box);
		g1=(float*)AllOnDev(M*n_box);
	#else
		Gg_f = new float[M*(N+1)*n_box]; Gg_b = new float[M*2*n_box];
		GG_F = new float[MM*N_A]; //Let N_A be the largest N of the solvents!!!!
		GN_A= H_GN_A;
		GN_B= H_GN_B;
		rho = H_rho;
		g1 = H_g1;
	#endif
		Free_energy = SCF(); 
		Free_energy=Helmholtz();
		printf("Free energy : %1f \n", Free_energy);
		int loc = 0;
		for (int iw=1; iw<MZ; iw++){
			if(H_PHI[25*JX+25*JY+ iw]>=0.5){
				loc = iw;
			}
		}
		varfile << Theta_A <<"     " << loc <<  "     " << interface << endl;
		}
		varfile.close();
	}
	else{

		int interface = 26;
		subdomain(H_Bx,H_By,H_Bz,H_Px,H_Py,H_Pz,interface);

	/*H_Bx[0]=15;
	H_By[0]=15;
	H_Bz[0]=15;
	H_Px[0]=26;
	H_Py[0]=26;
	H_Pz[0]=26;
	H_Bx[1]=10;
	H_By[1]=10;
	H_Bz[1]=15;
	H_Px[1]=20;
	H_Py[1]=20;
	H_Pz[1]=26;*/

		for (int p=0; p<n_box; p++){
			H_mask[p*M + jx*(H_Px[p]-H_Bx[p])+jy*(H_Py[p]-H_By[p])+(H_Pz[p]-H_Bz[p])]=1;
			H_MASK[((H_Px[p]-1)%MX+1)*JX + ((H_Py[p]-1)%MY+1)*JY + (H_Pz[p]-1)%MZ+1]=1;
		}

		invert(H_KSAM,H_MASK,MM);
		H_g1= new float[M*n_box]; H_Zero(H_g1,M*n_box);
		H_rho = new float[M*n_box];
		H_GN_A = new float[n_box];
		H_GN_B = new float[n_box];

	#ifdef CUDA
		Gg_f=(float*)AllOnDev(M*N*n_box);
		GG_F=(float*)AllOnDev(MM*N_A); //Let N_A be the largest N of the solvents!!!!
		Gg_b=(float*)AllOnDev(M*2*n_box);
		GN_A=(float*)AllOnDev(n_box);
		GN_B=(float*)AllOnDev(n_box);
		rho =(float*)AllOnDev(M*n_box);
		g1=(float*)AllOnDev(M*n_box);
	#else
		Gg_f = new float[M*(N+1)*n_box]; Gg_b = new float[M*2*n_box];
		GG_F = new float[MM*N_A]; //Let N_A be the largest N of the solvents!!!!
		GN_A= H_GN_A;
		GN_B= H_GN_B;
		rho = H_rho;
		g1 = H_g1;
	#endif	
		
		Free_energy = SCF(); 
		Free_energy=Helmholtz();
		printf("Free energy : %1f \n", Free_energy);	
	}



#ifdef CUDA
	TransferDataToHost(H_u,u,MM);
#endif
	vtk_output(fname+"_pot.vtk",H_u);
#ifdef CUDA
	TransferDataToHost(H_phi,phi,4*MM);
#endif

	Cp(PHI,H_phi+2*MM,MM); Add(PHI,phi+3*MM,MM); 
	vtk_output(fname+"_phi.vtk",PHI);

/*	ofstream varfile;
	varfile.open ("surf_profile.dat");
	for (int z=1; z<MZ; z++){	
		varfile << z <<"\t" << PHI[10*JX+40*JY+ z] << endl;
		}
	varfile.close();
	
	ofstream varfile2;
	varfile2.open ("surf_profile_part_1.dat");
	for (int z=1; z<MZ; z++){	
		varfile2 << z <<"\t" << PHI[15*JX+35*JY+ z] << endl;
		}
	varfile2.close();
*/	
	
	
	
	Cp(PHI,H_phi,MM);
	vtk_output(fname+"A_phi.vtk",PHI);
	
/*	ofstream varfile3;
	varfile3.open ("mol_profile.dat");
	for (int z=1; z<MZ; z++){	
		varfile3 << z <<"\t" << PHI[10*JX+40*JY+ z] << endl;
		}
	varfile3.close();
	
	ofstream varfile4;
	varfile4.open ("mol_profile_part_1.dat");
	for (int z=1; z<MZ; z++){	
		varfile4 << z <<"\t" << PHI[15*JX+35*JY+ z] << endl;
		}
	varfile4.close();
*/	
	
	Cp(PHI,H_phi+MM,MM);
	vtk_output(fname+"B_phi.vtk",PHI);
	filename=fname+".out";
	out_file.open(filename.c_str());
	out_file << "Free_energy : " << Free_energy << endl;
	out_file.close();

	return 0;
}

