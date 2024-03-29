#include <iostream>
#include "lattice.h"
#include "LG1Planar.h"

LG1Planar::LG1Planar(vector<Input*> In_,string name_): LGrad1(In_,name_) {}

LG1Planar::~LG1Planar() {
if (debug) cout <<"LG1Planar destructor " << endl;
}

void LG1Planar:: ComputeLambdas() {
	for (int i=1; i<MX+1; i++) L[i]=1;

	if (fjc>1) {
		for (int i = 0; i < M; i++) {
			L[i] = 1.0/fjc;
			LAMBDA[i] = 1.0/(2*(FJC-1));
			LAMBDA[i + (FJC-1)*M] = 1.0/(2*(FJC-1));
			LAMBDA[i+(FJC-1)/2*M] = 1.0/(FJC-1);
			for (int j = 1; j < FJC/2; j++) {
				LAMBDA[i+j*M] = 1.0/(FJC-1);
				LAMBDA[i+(FJC-j-1)*M] = 1.0/(FJC-1);
			}
		}
	}
}

void LG1Planar::Side(Real *X_side, Real *X, int M) {
if (debug) cout <<" Side in LG1Planar " << endl;
	if (ignore_sites) {
		Cp(X_side,X,M); return;
	}
	Zero(X_side,M); //set_bounds(X);
	int kk;

	if (fcc_sites) {
		Add(X_side+1,X,M-1);
		Add(X_side,X+1,M-1);
		Add(X_side,X,M);
		Norm(X_side,1.0/3.0,M);
	} else {
		if (fjc==1) {
			YplusisCtimesX(X_side+1,X,lambda,M-1);
			YplusisCtimesX(X_side,X+1,lambda,M-1);
			YplusisCtimesX(X_side,X,1.0-2.0*lambda,M);
		} else {
			for (int j = 0; j < FJC/2; j++) {
				kk = (FJC-1)/2-j;
				AddTimes(X_side+kk, X, LAMBDA+j*M+kk, M-kk);
				AddTimes(X_side, X+kk, LAMBDA+(FJC-j-1)*M, M-kk);
			}
			AddTimes(X_side, X, LAMBDA+(FJC-1)/2*M, M);
		}
	}
}

void LG1Planar::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
if (debug) cout <<" propagateF in LG1Planar " << endl;

	Real *gs = G+M*FJC*(s_to), *gs_1 = G+M*FJC*(s_from);
	Real *g = G1;

	Zero (gs,M*FJC);
	for (int k=0; k<(FJC-1)/2; k++) set_bounds(gs_1+k*M,gs_1+(FJC-k-1)*M);
	set_bounds(gs_1+(FJC-1)/2*M);

	if (lattice_type==simple_cubic) {
		Real *gz0 = gs_1, *gz1 = gs_1+M, *gz2 = gs_1+2*M;
		Real *gx0 = gs, *gx1 = gs+M, *gx2 = gs+2*M;

		YplusisCtimesX(gx0+1,gz0,P[0],     M-1);
		YplusisCtimesX(gx0+1,gz1,4*P[1],   M-1);

		YplusisCtimesX(gx1,gz0,P[1],       M);
		YplusisCtimesX(gx1,gz1,2*P[1]+P[0],M);
		YplusisCtimesX(gx1,gz2,P[1],       M);

		YplusisCtimesX(gx2,gz1+1,4*P[1],   M-1);
		YplusisCtimesX(gx2,gz2+1,P[0],     M-1);

		for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
	} else {
		int a,b; Real c;

		for (int p=0; p<FJC; p++){
			a=p-fjc; if (a<0) {b=0; a=-a; } else {b=a; a=0;}
			for (int q=0; q<FJC; q++) {
				c=P[abs(-p+q)];
				if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
				if (c!=0) YplusisCtimesX(gs+p*M+a,gs_1+q*M+b,c,M-a-b);
			}
		}
		for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
	}
}

void LG1Planar::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
if (debug) cout <<" propagateB in LG1Planar " << endl;

	Real *gs = G+M*FJC*(s_to), *gs_1 = G+M*FJC*(s_from);
	Real *g = G1;

	Zero (gs,M*FJC);
	for (int k=0; k<(FJC-1)/2; k++) set_bounds(gs_1+k*M,gs_1+(FJC-k-1)*M);
	set_bounds(gs_1+(FJC-1)/2*M);

	if (lattice_type==simple_cubic) {
		Real *gz0 = gs_1, *gz1 = gs_1+M, *gz2 = gs_1+2*M;
		Real *gx0 = gs,   *gx1 = gs+M,   *gx2 = gs+2*M;

		YplusisCtimesX(gx0,  gz0+1,P[0],       M-1);
		YplusisCtimesX(gx0,  gz1,  4*P[1],     M);

		YplusisCtimesX(gx1,  gz0+1,P[1],       M-1);
		YplusisCtimesX(gx1,  gz1,  2*P[1]+P[0],M);
		YplusisCtimesX(gx1+1,gz2,  P[1],       M-1);

		YplusisCtimesX(gx2,  gz1,  4*P[1],     M);
		YplusisCtimesX(gx2+1,gz2,  P[0],       M-1);

		for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
	} else {
		int a,b; Real c;

		for (int q=FJC-1; q>-1; q--){
			a=q-fjc; if (a>0) {b=0;} else {b=-a; a=0;}
			for (int p=FJC-1; p>-1; p--) {
				c=P[abs(-p+q)];
				if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
				if (c!=0) YplusisCtimesX(gs+p*M+a,gs_1+q*M+b,c,M-a-b);
			}
		}
		for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
	}

}

void LG1Planar::propagate(Real *G, Real *G1, int s_from, int s_to,int M) {
if (debug) cout <<" propagate in LG1Planar " << endl;
	Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	int kk;
	int j;
	Zero (gs,M); set_bounds(gs_1);

	if (fjc==1) {
		YplusisCtimesX(gs+1,gs_1,lambda,M-1);
		YplusisCtimesX(gs,gs_1,1.0-2.0*lambda,M);
		YplusisCtimesX(gs,gs_1+1,lambda,M-1);
		Times(gs,gs,G1,M);
	} else {
		for (j = 0; j < FJC/2; j++) {
			kk = (FJC-1)/2-j;
			AddTimes(gs+kk, gs_1, LAMBDA+j*M+kk, M-kk);
			AddTimes(gs, gs_1+kk, LAMBDA+(FJC-j-1)*M, M-kk);
		}
		AddTimes(gs, gs_1, LAMBDA+(FJC-1)/2*M, M);
		Times(gs, gs, G1, M);
	}
}


void LG1Planar::UpdateEE(Real* EE, Real* psi, Real* E) {
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_M_bounds(psi);
	Zero(EE,M);
	Real Exmin,Explus;

	pf =pf/2.0*fjc*fjc;
	Explus=psi[fjc-1]-psi[fjc];
	Explus *=Explus;

	for (int x=fjc; x<MX+fjc; x++) {
		Exmin=Explus;
		Explus=psi[x]-psi[x+1];
		Explus *=Explus;
		EE[x]=pf*(Exmin+Explus);
	}
}


void LG1Planar::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, int* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	Real a,b,c,a_,b_,c_;
	Real epsXplus, epsXmin;
	//set_M_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);

   if (!fixedPsi0) {
	C=C*2.0/fjc/fjc;
	epsXplus=eps[fjc-1]+eps[fjc];
	a=0; b=psi[fjc-1]; c=psi[fjc];
	for (int x=fjc; x<MX+fjc; x++) {
		epsXmin=epsXplus;
		epsXplus=eps[x]+eps[x+1];
		//a=b; b=c; c=psi[x+1];
		//X[x]=(epsXmin*a  +C*q[x] + epsXplus*c)/(epsXmin+epsXplus);
		if (x==fjc) a=psi[fjc-1]; else a=X[x-1]; //upwind
		X[x]=(epsXmin*a  +C*q[x] + epsXplus*psi[x+1])/(epsXmin+epsXplus);
	}
	YisAminB(g,g,X,M);
   } else { //fixedPsi0 is true;
	a=0; b=psi[fjc-1]; c=psi[fjc];
	for (int x=fjc; x<MX+fjc; x++) {
		a=b; b=c; c=psi[x+1];
		if (Mask[x] == 0) psi[x]=0.5*(a+c)+q[x]*C/eps[x];
	}
	if (grad_epsilon) {
		a=0; b=psi[fjc-1]; c=psi[fjc];a_=0; b_=eps[fjc-1]; c_=eps[fjc];
		for (int x=fjc; x<MX+fjc; x++) {//for all geometries
			a=b; b=c; c=psi[x+1]; a_=b_; b_=c_; c_=eps[x+1];
			if (Mask[x] == 0) {
				psi[x]+=0.25*(c_-a_)*(c-a)/eps[x]*fjc*fjc;
			}
		}
	}
	for (int x=fjc; x<MX+fjc; x++)
	if (Mask[x] == 0) {
		g[x]-=psi[x];
	}
   }
}


void LG1Planar::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, int* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
	Real a,b,c,a_,b_,c_;

	Real C = -e*e/(eps0*k_BT*bond_length);
	a=0; b=psi[fjc-1]; c=psi[fjc];
	for (int x=fjc; x<MX+fjc; x++) { //for all geometries
		a=b; b=c; c=psi[x+1];
		if (Mask[x] == 1) q[x] = -0.5*(a-2*b+c)*fjc*fjc*eps[x]/C;
	}

	if (grad_epsilon) {
		a=0; b=psi[fjc-1]; c=psi[fjc]; a_=0; b_=eps[fjc-1]; c_=eps[fjc];
		for (int x=fjc; x<MX+fjc; x++) {//for all geometries
			a=b; b=c; c=psi[x+1]; a_=b_; b_=c_; c_=eps[x+1];
			if (Mask[x] == 1) q[x]-=0.25*(c_-a_)*(c-a)*fjc*fjc/C;
		}
	}
	for (int x=fjc; x<MX+fjc; x++)
	if (Mask[x] == 1) {
		g[x]=-q[x];
	}
}
bool LG1Planar:: PutMask(int* MASK,vector<int>px,vector<int>py,vector<int>pz,int R){
	bool success=true;
	cout <<"PutMask does not make sence in planar 1 gradient system " << endl;
	return success;
}

Real LG1Planar::DphiDt(Real *g, Real* B_phitot, Real* phiA, Real* phiB, Real* alphaA, Real* alphaB, Real B_A, Real B_B) {
	Real AverageJ=0;
	Real a,b,c,Ma,Mb,Mc;
	g[1]=phiA[0]/phiA[1]-1.0;
	b=phiA[1]*phiB[1]*B_B/B_phitot[1];
	c=phiA[2]*phiB[2]*B_B/B_phitot[2];
	Mb=alphaA[1]-alphaB[1];
	Mc=alphaA[2]-alphaB[2];
	for (int z=2; z<M-2; z++) {
		a=b; b=c; c=phiA[z+1]*phiB[z+1]*B_B/B_phitot[z+1];
		Ma=Mb; Mb=Mc; Mc=alphaA[z+1]-alphaB[z+1];
		g[z] = g[z]  + B_A*((a+b)*(Mb-Ma)-(b+c)*(Mc-Mb));
		AverageJ+=(a+b)*(Mb-Ma);
	}
	AverageJ/=(M-4);

        //for (int z=2; z<M-2; z++) g[z]-=AverageJ;

	g[M-2]=phiA[M-1]/phiA[M-2]-1.0;
 	//cout <<"J " << -B_A*AverageJ/(2*(M-4)) << endl;

	return -B_A*AverageJ/2;// /(2*(M-4));
}


/* //Forward propagator as is was before the compactation.
	switch (fjc) {
		case 1:
			if (lattice_type ==hexagonal) {

				int a,b; Real c;
				for (int p=0; p<FJC; p++){
					a=p-fjc; if (a<0) {b=0; a=-a; } else {b=a; a=0;}
					for (int q=0; q<FJC; q++) {
						c=P[abs(-p+q)];
						if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
						if (c!=0) YplusisCtimesX(gs+p*M+a,gs_1+q*M+b,c,M-a-b);
					}}

						for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);

//				YplusisCtimesX(gx0+1,gz0,P[0],  M-1);
//				YplusisCtimesX(gx0+1,gz1,2*P[1],M-1);

//				YplusisCtimesX(gx1,  gz0,P[1],  M);
//				YplusisCtimesX(gx1,  gz1,P[0],  M);
//				YplusisCtimesX(gx1,  gz2,P[1],  M);

//				YplusisCtimesX(gx2,gz1+1,2*P[1],M-1);
//				YplusisCtimesX(gx2,gz2+1,P[0],  M-1);



			} else { //simple cubic

				YplusisCtimesX(gx0+1,gz0,P[0],     M-1);
				YplusisCtimesX(gx0+1,gz1,4*P[1],   M-1);

				YplusisCtimesX(gx1,gz0,P[1],       M);
				YplusisCtimesX(gx1,gz1,2*P[1]+P[0],M);
				YplusisCtimesX(gx1,gz2,P[1],       M);

				YplusisCtimesX(gx2,gz1+1,4*P[1],   M-1);
				YplusisCtimesX(gx2,gz2+1,P[0],     M-1);

				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
			}
			break;

		case 2:
			if (lattice_type==hexagonal) {
				//Real *gx3 = gs+3*M,   *gx4 = gs+4*M;
				//Real *gz3 = gs_1+3*M, *gz4 = gs_1+4*M;

				int a,b; Real c;
				for (int p=0; p<FJC; p++){
					a=p-fjc; if (a<0) {b=0; a=-a; } else {b=a; a=0;}
					for (int q=0; q<FJC; q++) {
						c=P[abs(-p+q)];
						if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
						if (c!=0) YplusisCtimesX(gs+p*M+a,gs_1+q*M+b,c,M-a-b);
					}
				}

				YplusisCtimesX(gx0+2,gz0,P[0],     M-2);
				YplusisCtimesX(gx0+2,gz1,2*P[1],   M-2);
				YplusisCtimesX(gx0+2,gz2,2*P[2],   M-2);
				YplusisCtimesX(gx0+2,gz3,2*P[3],   M-2);

				YplusisCtimesX(gx1+1,gz0,P[1],     M-1);
				YplusisCtimesX(gx1+1,gz1,P[0]+P[2],M-1);
				YplusisCtimesX(gx1+1,gz2,P[1]+P[3],M-1);
				YplusisCtimesX(gx1+1,gz3,P[2],     M-1);
				YplusisCtimesX(gx1+1,gz4,P[3],     M-1);

				YplusisCtimesX(gx2,  gz0,P[2],     M);
				YplusisCtimesX(gx2,  gz1,P[1]+P[3],M);
				YplusisCtimesX(gx2,  gz2,P[0],     M);
				YplusisCtimesX(gx2,  gz3,P[1]+P[3],M);
				YplusisCtimesX(gx2,  gz4,P[2],     M);

				YplusisCtimesX(gx3,gz0+1,P[3],     M-1);
				YplusisCtimesX(gx3,gz1+1,P[2],     M-1);
				YplusisCtimesX(gx3,gz2+1,P[1]+P[3],M-1);
				YplusisCtimesX(gx3,gz3+1,P[0]+P[2],M-1);
				YplusisCtimesX(gx3,gz4+1,P[1],     M-1);

				YplusisCtimesX(gx4,gz1+2,2*P[3],   M-2);
				YplusisCtimesX(gx4,gz2+2,2*P[2],   M-2);
				YplusisCtimesX(gx4,gz3+2,2*P[1],   M-2);
				YplusisCtimesX(gx4,gz4+2,P[0],     M-2);

				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
			} else {
				cout <<"cubic lattice and fjc=2 Markov 2 not implemented " << endl;
			}
			break;
		case 3:
			if (lattice_type==hexagonal) {
				//Real *gx3 = gs+3*M,   *gx4 = gs+4*M,   *gx5 = gs+5*M,   *gx6 = gs+6*M;
				//Real *gz3 = gs_1+3*M, *gz4 = gs_1+4*M, *gz5 = gs_1+5*M, *gz6 = gs_1+6*M;

				int a,b; Real c;

				for (int p=0; p<FJC; p++){
					a=p-fjc; if (a<0) {b=0; a=-a; } else {b=a; a=0;}
					for (int q=0; q<FJC; q++) {
						c=P[abs(-p+q)];
						if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
						if (c!=0) YplusisCtimesX(gs+p*M+a,gs_1+q*M+b,c,M-a-b);
					}
				}

				YplusisCtimesX(gx0+3,gz0,P[0],     M-3);
				YplusisCtimesX(gx0+3,gz1,2*P[1],   M-3);
				YplusisCtimesX(gx0+3,gz2,2*P[2],   M-3);
				YplusisCtimesX(gx0+3,gz3,2*P[3],   M-3);
				YplusisCtimesX(gx0+3,gz4,2*P[4],   M-3);
				YplusisCtimesX(gx0+3,gz5,2*P[5],   M-3);

				YplusisCtimesX(gx1+2,gz0,P[1],     M-2);
				YplusisCtimesX(gx1+2,gz1,P[0]+P[2],M-2);
				YplusisCtimesX(gx1+2,gz2,P[1]+P[3],M-2);
				YplusisCtimesX(gx1+2,gz3,P[2]+P[4],M-2);
				YplusisCtimesX(gx1+2,gz4,P[3]+P[5],M-2);
				YplusisCtimesX(gx1+2,gz5,P[4],     M-2);
				YplusisCtimesX(gx1+2,gz6,P[5],     M-2);

				YplusisCtimesX(gx2+1,gz0,P[2],     M-1);
				YplusisCtimesX(gx2+1,gz1,P[1]+P[3],M-1);
				YplusisCtimesX(gx2+1,gz2,P[0]+P[4],M-1);
				YplusisCtimesX(gx2+1,gz3,P[1]+P[5],M-1);
				YplusisCtimesX(gx2+1,gz4,P[2],     M-1);
				YplusisCtimesX(gx2+1,gz5,P[3]+P[5],M-1);
				YplusisCtimesX(gx2+1,gz6,P[4],     M-1);

				YplusisCtimesX(gx3,  gz0,P[3],     M);
				YplusisCtimesX(gx3,  gz1,P[2]+P[4],M);
				YplusisCtimesX(gx3,  gz2,P[1]+P[5],M);
				YplusisCtimesX(gx3,  gz3,P[0]     ,M);
				YplusisCtimesX(gx3,  gz4,P[1]+P[5],M);
				YplusisCtimesX(gx3,  gz5,P[2]+P[4],M);
				YplusisCtimesX(gx3,  gz6,P[3],     M);

				YplusisCtimesX(gx4,gz0+1,P[4],     M-1);
				YplusisCtimesX(gx4,gz1+1,P[3]+P[5],M-1);
				YplusisCtimesX(gx4,gz2+1,P[2],     M-1);
				YplusisCtimesX(gx4,gz3+1,P[1]+P[5],M-1);
				YplusisCtimesX(gx4,gz4+1,P[0]+P[4],M-1);
				YplusisCtimesX(gx4,gz5+1,P[1]+P[3],M-1);
				YplusisCtimesX(gx4,gz6+1,P[2],     M-1);

				YplusisCtimesX(gx5,gz0+2,P[5],     M-2);
				YplusisCtimesX(gx5,gz1+2,P[4],     M-2);
				YplusisCtimesX(gx5,gz2+2,P[3]+P[5],M-2);
				YplusisCtimesX(gx5,gz3+2,P[2]+P[4],M-2);
				YplusisCtimesX(gx5,gz4+2,P[1]+P[3],M-2);
				YplusisCtimesX(gx5,gz5+2,P[0]+P[2],M-2);
				YplusisCtimesX(gx5,gz6+2,P[1],     M-2);

				YplusisCtimesX(gx6,gz1+3,2*P[5],   M-3);
				YplusisCtimesX(gx6,gz2+3,2*P[4],   M-3);
				YplusisCtimesX(gx6,gz3+3,2*P[3],   M-3);
				YplusisCtimesX(gx6,gz4+3,2*P[2],   M-3);
				YplusisCtimesX(gx6,gz5+3,2*P[1],   M-3);
				YplusisCtimesX(gx6,gz6+3,P[0],     M-3);

				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
			} else {
				cout <<"cubic lattice and FJC_choices=3 Markov 2 not implemented " << endl;
			}
			break;

		default:
			if (lattice_type==hexagonal) {
				int a,b; Real c;

				for (int p=0; p<FJC; p++){
					a=p-fjc; if (a<0) {b=0; a=-a; } else {b=a; a=0;}
					for (int q=0; q<FJC; q++) {
						c=P[abs(-p+q)];
						if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
						if (c!=0) YplusisCtimesX(gs+p*M+a,gs_1+q*M+b,c,M-a-b);
					}
				}
				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
			} else {
				cout <<"cubic lattice and FJC_choices=3 Markov 2 not implemented " << endl;
			}

			break;
	}
//Backward propagator before compactation.
	switch (fjc) {
		case 1:
			if (lattice_type ==hexagonal) {
				int a,b; Real c;

				for (int q=FJC-1; q>-1; q--){
					a=q-fjc; if (a>0) {b=0;} else {b=-a; a=0;}
					for (int p=FJC-1; p>-1; p--) {
						c=P[abs(-p+q)];
						if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
						if (c!=0) YplusisCtimesX(gs+p*M+a,gs_1+q*M+b,c,M-a-b);
					}
				}
				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);

				//YplusisCtimesX(gx0,  gz0+1,P[0],  M-1);
				//YplusisCtimesX(gx0,  gz1,  2*P[1],M);

				//YplusisCtimesX(gx1,  gz0+1,P[1],  M-1);
				//YplusisCtimesX(gx1,  gz1,  P[0],  M);
				//YplusisCtimesX(gx1+1,gz2,  P[1],  M-1);

				//YplusisCtimesX(gx2,  gz1,  2*P[1],M);
				//YplusisCtimesX(gx2+1,gz2,  P[0],  M-1);

				//for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);

			} else { //simple cubic
				YplusisCtimesX(gx0,  gz0+1,P[0],       M-1);
				YplusisCtimesX(gx0,  gz1,  4*P[1],     M);

				YplusisCtimesX(gx1,  gz0+1,P[1],       M-1);
				YplusisCtimesX(gx1,  gz1,  2*P[1]+P[0],M);
				YplusisCtimesX(gx1+1,gz2,  P[1],       M-1);

				YplusisCtimesX(gx2,  gz1,  4*P[1],     M);
				YplusisCtimesX(gx2+1,gz2,  P[0],       M-1);

				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
			}
			break;

		case 2:
			if (lattice_type ==hexagonal) {
				//Real *gx3 = gs+3*M,   *gx4 = gs+4*M;
				//Real *gz3 = gs_1+3*M, *gz4 = gs_1+4*M;
				int a,b; Real c;

				for (int q=FJC-1; q>-1; q--){
					a=q-fjc; if (a>0) {b=0;} else {b=-a; a=0;}
					for (int p=FJC-1; p>-1; p--) {
						c=P[abs(-p+q)];
						if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
						if (c!=0) YplusisCtimesX(gs+p*M+a,gs_1+q*M+b,c,M-a-b);
					}
				}

				YplusisCtimesX(gx1+2,gz4,   P[3],     M-2);
				YplusisCtimesX(gx2+2,gz4,   P[2],     M-2);
				YplusisCtimesX(gx3+2,gz4,   P[1],     M-2);
				YplusisCtimesX(gx4+2,gz4,   P[0],     M-2);

				YplusisCtimesX(gx0+1,gz3,   2*P[3],   M-1);
				YplusisCtimesX(gx1+1,gz3,   P[2],     M-1);
				YplusisCtimesX(gx2+1,gz3,   P[1]+P[3],M-1);
				YplusisCtimesX(gx3+1,gz3,   P[0]+P[2],M-1);
				YplusisCtimesX(gx4+1,gz3,   2*P[1],   M-1);

				YplusisCtimesX(gx0,  gz2,   2*P[2],   M);
				YplusisCtimesX(gx1,  gz2,   P[1]+P[3],M);
				YplusisCtimesX(gx2,  gz2,   P[0],     M);
				YplusisCtimesX(gx3,  gz2,   P[1]+P[3],M);
				YplusisCtimesX(gx4,  gz2,   2*P[2],   M);

				YplusisCtimesX(gx0,  gz1+1, 2*P[1],   M-1);
				YplusisCtimesX(gx1,  gz1+1, P[0]+P[2],M-1);
				YplusisCtimesX(gx2,  gz1+1, P[1]+P[3],M-1);
				YplusisCtimesX(gx3,  gz1+1, P[2],     M-1);
				YplusisCtimesX(gx4,  gz1+1, 2*P[3],   M-1);

				YplusisCtimesX(gx0,  gz0+2, P[0],     M-2);
				YplusisCtimesX(gx1,  gz0+2, P[1],     M-2);
				YplusisCtimesX(gx2,  gz0+2, P[2],     M-2);
				YplusisCtimesX(gx3,  gz0+2, P[3],     M-2);

				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
			} else {
				cout <<"cubic lattice type in FJC_choices>3 not inplemented " << endl;
			}
			break;
		case 3:
			if (lattice_type ==hexagonal) {
				//Real *gx3 = gs+3*M,   *gx4 = gs+4*M,  *gx5 = gs+5*M,  *gx6 = gs+6*M;
				//Real *gz3 = gs_1+3*M, *gz4 = gs_1+4*M,*gz5 = gs_1+5*M,*gz6 = gs_1+6*M;

				int a,b; Real c;

				for (int q=FJC-1; q>-1; q--){
					a=q-fjc; if (a>0) {b=0;} else {b=-a; a=0;}
					for (int p=FJC-1; p>-1; p--) {
						c=P[abs(-p+q)];
						if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
						if (c!=0) YplusisCtimesX(gs+p*M+a,gs_1+q*M+b,c,M-a-b);
						//cout <<"p " << p << " q " << q << " a " << a << " b " << b << "abs(-p+q) " << abs(-p+q) << " FJC-1-abs(FJC-1-p-q) " << FJC-1-abs(FJC-1-p-q) << endl;
					}
				}

				YplusisCtimesX(gx6+3,gz6,   P[0],     M-3);
				YplusisCtimesX(gx5+3,gz6,   P[1],     M-3);
				YplusisCtimesX(gx4+3,gz6,   P[2],     M-3);
				YplusisCtimesX(gx3+3,gz6,   P[3],     M-3);
				YplusisCtimesX(gx2+3,gz6,   P[4],     M-3);
				YplusisCtimesX(gx1+3,gz6,   P[5],     M-3);

				YplusisCtimesX(gx6+2,gz5,   2*P[1],   M-2);
				YplusisCtimesX(gx5+2,gz5,   P[0]+P[2],M-2);
				YplusisCtimesX(gx4+2,gz5,   P[1]+P[3],M-2);
				YplusisCtimesX(gx3+2,gz5,   P[2]+P[4],M-2);
				YplusisCtimesX(gx2+2,gz5,   P[3]+P[5],M-2);
				YplusisCtimesX(gx1+2,gz5,   P[4],     M-2);
				YplusisCtimesX(gx0+2,gz5,   2*P[5],   M-2);

				YplusisCtimesX(gx6+1,gz4,   2*P[2],   M-1);
				YplusisCtimesX(gx5+1,gz4,   P[1]+P[3],M-1);
				YplusisCtimesX(gx4+1,gz4,   P[0]+P[4],M-1);
				YplusisCtimesX(gx3+1,gz4,   P[1]+P[5],M-1);
				YplusisCtimesX(gx2+1,gz4,   P[2],     M-1);
				YplusisCtimesX(gx1+1,gz4,   P[3]+P[5],M-1);
				YplusisCtimesX(gx0+1,gz4,   2*P[4],   M-1);

				YplusisCtimesX(gx6,  gz3,   2*P[3],   M);
				YplusisCtimesX(gx5,  gz3,   P[2]+P[4],M);
				YplusisCtimesX(gx4,  gz3,   P[1]+P[5],M);
				YplusisCtimesX(gx3,  gz3,   P[0],     M);
				YplusisCtimesX(gx2,  gz3,   P[1]+P[5],M);
				YplusisCtimesX(gx1,  gz3,   P[2]+P[4],M);
				YplusisCtimesX(gx0,  gz3,   2*P[3],   M);

				YplusisCtimesX(gx6,  gz2+1, 2*P[4],   M-1);
				YplusisCtimesX(gx5,  gz2+1, P[3]+P[5],M-1);
				YplusisCtimesX(gx4,  gz2+1, P[2],     M-1);
				YplusisCtimesX(gx3,  gz2+1, P[1]+P[5],M-1);
				YplusisCtimesX(gx2,  gz2+1, P[0]+P[4],M-1);
				YplusisCtimesX(gx1,  gz2+1, P[1]+P[3],M-1);
				YplusisCtimesX(gx0,  gz2+1, 2*P[2],   M-1);

				YplusisCtimesX(gx6,  gz1+2, 2*P[5],   M-2);
				YplusisCtimesX(gx5,  gz1+2, P[4],     M-2);
				YplusisCtimesX(gx4,  gz1+2, P[3]+P[5],M-2);
				YplusisCtimesX(gx3,  gz1+2, P[2]+P[4],M-2);
				YplusisCtimesX(gx2,  gz1+2, P[1]+P[3],M-2);
				YplusisCtimesX(gx1,  gz1+2, P[0]+P[2],M-2);
				YplusisCtimesX(gx0,  gz1+2, 2*P[1],   M-2);

				YplusisCtimesX(gx5,  gz0+3, P[5],     M-3);
				YplusisCtimesX(gx4,  gz0+3, P[4],     M-3);
				YplusisCtimesX(gx3,  gz0+3, P[3],     M-3);
				YplusisCtimesX(gx2,  gz0+3, P[2],     M-3);
				YplusisCtimesX(gx1,  gz0+3, P[1],     M-3);
				YplusisCtimesX(gx0,  gz0+3, P[0],     M-3);

				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
			} else {
				cout <<"cubic lattice type in FJC_choices>3 not inplemented " << endl;
			}
			break;

		default:
			if (lattice_type ==hexagonal) {
				int a,b; Real c;

				for (int q=FJC-1; q>-1; q--){
					a=q-fjc; if (a>0) {b=0;} else {b=-a; a=0;}
					for (int p=FJC-1; p>-1; p--) {
						c=P[abs(-p+q)];
						if (q>0 && q<FJC-1) c+= P[FJC-1-abs(FJC-1-p-q)];
						if (c!=0) YplusisCtimesX(gs+p*M+a,gs_1+q*M+b,c,M-a-b);
					}
				}
				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
			} else {
				cout <<"cubic lattice type in FJC_choices>3 not inplemented " << endl;
			}
			break;
	}
*/
