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


void LG1Planar::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
if (debug) cout <<" Side in LG1Planar " << endl;

	if (ignore_sites) {
		Cp(X_side,X,M); return;
	}
	Zero(X_side,M);set_bounds(X);
	int kk;

	if (fcc_sites) {
		YplusisCtimesX(X_side+1,X,1.0/3.0,M-1);
		YplusisCtimesX(X_side,X+1,1.0/3.0,M-1);
		YplusisCtimesX(X_side,X,1.0/3.0,M);
	} else {
		if (fjc==1) {
			YplusisCtimesX(X_side+1,X,l1,M-1);
			YplusisCtimesX(X_side,X+1,l1,M-1);
			YplusisCtimesX(X_side,X,l0,M);
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

void LG1Planar::propagate(Real *G, Real *G1, int s_from, int s_to,int M) { //this procedure should function on simple cubic lattice.
if (debug) cout <<" propagate in LG1Planar " << endl;
	Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	int kk;
	int j;

	if (fjc==1) {
		YplusisCtimesX(gs+1,gs_1,l1,M-1);
		YplusisCtimesX(gs,gs_1,l0,M);
		YplusisCtimesX(gs,gs_1+1,l1,M-1);
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
	set_bounds(psi);
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
	set_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);

   if (!fixedPsi0) {
	C=C*2.0/fjc/fjc;
	epsXplus=eps[fjc-1]+eps[fjc];
	a=0; b=psi[fjc-1]; c=psi[fjc];
	for (int x=fjc; x<MX+fjc; x++) {
		epsXmin=epsXplus;
		epsXplus=eps[x]+eps[x+1];
		a=b; b=c; c=psi[x+1];
		X[x]=(epsXmin*a  +C*q[x] + epsXplus*c)/(epsXmin+epsXplus);
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

