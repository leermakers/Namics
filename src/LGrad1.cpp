#include <iostream>
#include <string>
#include "lattice.h"
#include "LGrad1.h"

//planar geometry is in LG1Planar.cpp

LGrad1::LGrad1(vector<Input*> In_,string name_): Lattice(In_,name_) {
if (debug) cout <<"LGrad1 constructor " << endl;
}

LGrad1::~LGrad1() {
if (debug) cout <<"LGrad1 destructor " << endl;
}

void LGrad1:: ComputeLambdas() {
if (debug) cout <<"LGrad1 computeLambda's " << endl;

	Real r, VL, LS;
	Real rlow, rhigh;

	if (fcc_sites){
		if (geometry=="cylindrical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=PIE*(pow(r,2)-pow(r-1,2));
				lambda1[i]=2.0*PIE*r/L[i]/3.0;
				lambda_1[i]=2.0*PIE*(r-1)/L[i]/3.0;
				lambda0[i]=1.0/3.0;
			}
		}
		if (geometry=="spherical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=4.0/3.0*PIE*(pow(r,3)-pow(r-1,3));
				fcc_lambda1[i]=4.0*PIE*pow(r,2)/L[i]/3.0;
				fcc_lambda_1[i]=4.0*PIE*pow(r-1,2)/L[i]/3.0;
				fcc_lambda0[i]=1.0-fcc_lambda1[i]-fcc_lambda_1[i];
			}
		}
	}

	if (fjc==1) {
		if (geometry=="cylindrical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=PIE*(pow(r,2)-pow(r-1,2));
				lambda1[i]=2.0*PIE*r/L[i]*lambda;
				lambda_1[i]=2.0*PIE*(r-1)/L[i]*lambda;
				lambda0[i]=1.0-2.0*lambda;
			}
		}
		if (geometry=="spherical") {
			for (int i=1; i<MX+1; i++) {
				r=offset_first_layer + i;
				L[i]=4.0/3.0*PIE*(pow(r,3)-pow(r-1,3));
				lambda1[i]=4.0*PIE*pow(r,2)/L[i]*lambda;
				lambda_1[i]=4.0*PIE*pow(r-1,2)/L[i]*lambda;
				lambda0[i]=1.0-lambda1[i]-lambda_1[i];

			}
		}
		if (Markov==2) {
			for (int i=1; i<MX+1; i++) {
				l1[i]=lambda1[i]/lambda; l11[i]=1.0-l1[i];
				l_1[i]=lambda_1[i]/lambda; l_11[i]=1.0-l_1[i];
			}
		}

	}

	if (fjc>1) {
		if (geometry == "cylindrical") {
			for (int i = fjc; i < M - fjc; i++) {
				r = offset_first_layer+1.0*(i-fjc+1.0)/fjc;
				rlow = r - 0.5;
				rhigh = r + 0.5;
				L[i] = PIE * (2.0 * r) / fjc;
				VL = L[i] / PIE * fjc;
				if ((rlow - r) * 2 + r > 0.0)
					LAMBDA[i] += 1.0/(1.0*FJC-1.0)*rlow/VL;
				if ((rhigh - r) * 2 + r < 1.0*MX/fjc)
					LAMBDA[i+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
				else {
					if (2*rhigh-r-1.0*MX/fjc > -0.001 && 2 * rhigh-r-1.0*MX/fjc < 0.001) {
						LAMBDA[i+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
					}
					for (int j = 1; j <= fjc; j++) {
						if (2*rhigh-r-1.0*MX/fjc > 0.99*j/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*j/fjc) {
							LAMBDA[i+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*(rhigh-1.0*j/fjc)/VL;
						}
					}
				}
				for (int j = 1; j < fjc; j++) {
					rlow += 0.5/(fjc);
					rhigh -= 0.5/(fjc);
					if ((rlow-r)*2+r > 0.0)
						LAMBDA[i+j*M] += 1.0/(1.0*FJC-1.0)*2.0*rlow/VL;
					if ((rhigh-r)*2+r < offset_first_layer+1.0*MX/fjc)
						LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
					else {
						if (2 * rhigh-r-1.0*MX/fjc > -0.001 && 2*rhigh-r-1.0*MX/fjc < 0.001) {
							LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
						}
						for (int k = 1; k <= fjc; k++) {
							if (2 * rhigh-r-1.0*MX/fjc > 0.99*k/fjc && 2*rhigh-r-1.0*MX/fjc<1.01*k/fjc) {
								LAMBDA[i + (FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*(rhigh-1.0*k/fjc)/VL;
							}
						}
					}
				}
				LS = 0;
				for (int j = 0; j < FJC; j++)
					LS += LAMBDA[i+j*M];
				LAMBDA[i+(FJC/2)*M] += 1.0 - LS;
			}
		}

		if (geometry == "spherical") {
			for (int i = fjc; i < M - fjc; i++) {
				r = offset_first_layer+1.0*(1.0*i-1.0*fjc+1.0)/fjc;
				rlow = r-0.5;
				rhigh = r+0.5;
				L[i] = PIE*4.0/3.0*(rhigh*rhigh*rhigh-rlow*rlow*rlow)/fjc;
				VL = L[i] / PIE * fjc;
				if ((rlow-r)*2+r > 0.0)
					LAMBDA[i] += 0.5/(1.0*FJC-1.0)*4.0*rlow*rlow/VL;
				if ((rhigh -r)*2+r < 1.0*MX/fjc)
					LAMBDA[i+(FJC-1)*M] += 0.5/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
				else {
					if (2*rhigh-r-1.0*MX/fjc>-0.001 && 2*rhigh-r-1.0*MX/fjc<0.001) {
						LAMBDA[i+(FJC-1)*M] += 0.5/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
					}
					for (int j = 1; j <= fjc; j++) {
						if (2*rhigh-r-1.0*MX/fjc > 0.99*j/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*j/fjc) {
							LAMBDA[i+(FJC-1)*M] += 0.5/(1.0*FJC-1.0)*4.0*(rhigh-1.0*j/fjc)*(rhigh-1.0*j/fjc)/VL;
						}
					}
				}
				for (int j = 1; j < fjc; j++) {
					rlow += 0.5/(fjc);
					rhigh -= 0.5/(fjc);
					if ((rlow-r)*2+r > 0.0)
						LAMBDA[i+j*M] += 1.0/(1.0*FJC-1.0)*4.0*rlow*rlow/VL;
					if ((rhigh - r) * 2 + r < offset_first_layer + 1.0*MX/fjc)
						LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
					else {
						if (2*rhigh-r-1.0*MX/fjc > -0.001 && 2*rhigh-r-1.0*MX/fjc < 0.001) {
							LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
						}
						for (int k = 1; k <= fjc; k++) {
							if (2*rhigh-r-1.0*MX/fjc > 0.99*k/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*k/fjc) {
								LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*4.0*(rhigh-1.0*k/fjc)*(rhigh-1.0*k/fjc)/VL;
							}
						}
					}
				}
				LS = 0;
				for (int j = 0; j < FJC; j++)
					LS += LAMBDA[i+j*M];
				LAMBDA[i+(FJC/2)*M] += 1.0-LS;
			}
		}
		if (Markov==2) { //planar gemaakt voor debugging....
			for (int i = fjc; i < M - fjc; i++) {
				LABDA[i]=LAMBDA[i]*8.0;
				LABDA_1[i]=1.0-LABDA[i];
				for (int j=1; j<FJC-1; j++) {
					LABDA[i+j*M]=LAMBDA[i+j*M]*4.0;
					LABDA_1[i+j*M]=1.0-LABDA[i+j*M];
				}
				LABDA[i+(FJC-1)*M]=LAMBDA[i+(FJC-1)*M]*8.0;
				LABDA_1[i+(FJC-1)*M]=1.0-LABDA[i+(FJC-1)*M];
			}
		}
	}
}

bool LGrad1::PutM() {
if (debug) cout << "PutM in LGrad1 " << endl;
	bool success=true;

	JX=1; JY=0; JZ=0; M=MX+2*fjc;
	if (geometry=="planar") {volume = MX; }
	if (geometry=="spherical") {volume = 4/3*PIE*(pow(MX+offset_first_layer,3)-pow(offset_first_layer,3));}
	if (geometry=="cylindrical") {volume = PIE*(pow(MX+offset_first_layer,2)-pow(offset_first_layer,2));}

	Accesible_volume=volume;
	return success;
}

void LGrad1::TimesL(Real* X){
if (debug) cout << "TimesL in LGrad1 " << endl;
	if (geometry!="planar") Times(X,X,L,M);
}

void LGrad1::DivL(Real* X){
if (debug) cout << "DivL in LGrad1 " << endl;
	if (geometry!="planar") Div(X,L,M);
}

Real LGrad1:: Moment(Real* X,Real Xb, int n) {
if (debug) cout << "Moment in LGrad1 " << endl;
	Real Result=0;
	Real cor;
	remove_bounds(X);
	for (int i = fjc; i<M; i++) {
		cor = (i-fjc+0.5)/fjc;
		Result += pow(cor,n)*(X[i]-Xb)*L[i];
	}
	return Result/fjc;
}

Real LGrad1::WeightedSum(Real* X){
if (debug) cout << "weighted sum in LGrad1 " << endl;
	Real sum{0};
	remove_bounds(X);
	if (geometry=="planar") {
		Sum(sum,X,M); sum/=fjc;
	} else Dot(sum,X,L,M);
	return sum;
}

void LGrad1::vtk(string filename, Real* X, string id,bool writebounds) {
if (debug) cout << "vtk in LGrad1 " << endl;
	cout << "for system with one gradient there is no VTK output available " << endl;
}

void LGrad1::PutProfiles(FILE* pf,vector<Real*> X,bool writebounds){
if (debug) cout <<"PutProfiles in LGrad1 " << endl;
	int x,i;
	int length=X.size();
	int a;
	if (writebounds) a=0; else a = fjc;

	for (x=a; x<MX+2*fjc-a; x++){
		//fprintf(pf,"%e\t",offset_first_layer+1.0*x/fjc-1/(2.0*fjc));
#ifdef LongReal
		fprintf(pf,"%Le\t",offset_first_layer/fjc+1.0*(x-fjc+1)/fjc-0.5/fjc); //g - e
		for (i=0; i<length; i++)
		if (i<length-1) fprintf(pf,"%.20Lg\t",X[i][x]); else fprintf(pf,"%.20Lg",X[i][x]);
		fprintf(pf,"\n");
#else
		fprintf(pf,"%e\t",offset_first_layer/fjc+1.0*(x-fjc+1)/fjc-0.5/fjc); //g - e
		for (i=0; i<length; i++)
		if (i<length-1) fprintf(pf,"%.20g\t",X[i][x]); else fprintf(pf,"%.20g",X[i][x]);
		fprintf(pf,"\n");
#endif
	}
}

void LGrad1::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
if (debug) cout <<" Side in LGrad1 " << endl;

	if (ignore_sites) {
		Cp(X_side,X,M); return;
	}
	Zero(X_side,M);//set_bounds(X);
	int j, kk;

	if (fcc_sites) {
		AddTimes(X_side,X,fcc_lambda0,M);
		AddTimes(X_side+1,X,fcc_lambda_1+1,M-1);
		AddTimes(X_side,X+1,fcc_lambda1,M-1);

	} else {
		if (fjc==1) {
			AddTimes(X_side,X,lambda0,M);
			AddTimes(X_side+1,X,lambda_1+1,M-1);
			AddTimes(X_side,X+1,lambda1,M-1);
		} else {

			for (j = 0; j < FJC/2; j++) {
				kk = (FJC-1)/2-j;
				AddTimes(X_side+kk, X, LAMBDA+j*M+kk, M-kk);
				AddTimes(X_side, X+kk, LAMBDA+(FJC-j-1)*M, M-kk);
			}
			AddTimes(X_side, X, LAMBDA+(FJC-1)/2*M, M);

		}
	}
}

//void LGrad1::LReflect(Real *H, Real *P, Real *Q) {
//	Times   (H,   l_1 +1,  P,    M-1);
//	AddTimes(H,   l_11+1,  Q+1,  M-1);
//}

//void LGrad1::UReflect(Real *H, Real *P, Real* Q) {
//	Times   (H+1, l1,      P+1,  M-1);
//	AddTimes(H+1, l11,     Q,    M-1);
//}

void LGrad1::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	Real *gs=G+FJC*M*(s_to);
	Real *gs_1=G+FJC*M*(s_from);
	Real *gz0=gs_1;
	Real *gz1=gs_1+M;
	Real *gz2=gs_1+2*M;
	Real *gz3=gs_1+3*M;
	Real *gz4=gs_1+4*M;
	Real *gx0=gs;
	Real *gx1=gs+M;
	Real *gx2=gs+2*M;
	Real *gx3=gs+3*M;
	Real *gx4=gs+4*M;
	Real *g =G1;

	Zero (gs,M*FJC);
	for (int k=0; k<(FJC-1)/2; k++) set_bounds(gs_1+k*M,gs_1+(FJC-k-1)*M);
	set_bounds(gs_1+(FJC-1)/2*M);

	switch (fjc) {
		case 1:
			if (lattice_type==hexagonal) {
				//LReflect(H,gz0,gz2);
				Times   (H,   l_1 +1,  gz0,    M-1);
				AddTimes(H,   l_11+1,  gz2+1,  M-1);
				YplusisCtimesX(gx0+1,H,P[0],M-1);
				//LReflect(H,gz1,gz1);
				Times   (H,   l_1 +1,  gz1,    M-1);
				AddTimes(H,   l_11+1,  gz1+1,  M-1);
				YplusisCtimesX(gx0+1,H,2*P[1],M-1);

				YplusisCtimesX(gx1,gz0,P[1],M);
				YplusisCtimesX(gx1,gz1,P[0],M);
				YplusisCtimesX(gx1,gz2,P[1],M);

				//UReflect(H,gz1,gz1);
				Times   (H+1, l1,      gz1+1,  M-1);
				AddTimes(H+1, l11,     gz1,    M-1);
				YplusisCtimesX(gx2,H+1,2*P[1],M-1);
				//UReflect(H,gz2,gz0);
				Times   (H+1, l1,      gz2+1,  M-1);
				AddTimes(H+1, l11,     gz0,    M-1);
				YplusisCtimesX(gx2,H+1,P[0],M-1);

				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);

			} else {
				//LReflect(H,gz0,gz2);
				  Times   (H,   l_1 +1,  gz0,    M-1);
				  AddTimes(H,   l_11+1,  gz2+1,  M-1);
				YplusisCtimesX(gx0+1,H,P[0],M-1);

				//LReflect(H,gz1,gz1);
				  Times   (H,   l_1 +1,  gz1,    M-1);
				  AddTimes(H,   l_11+1,  gz1+1,  M-1);
				YplusisCtimesX(gx0+1,H,4*P[1],M-1);

				YplusisCtimesX(gx1,gz0,P[1],M);
				YplusisCtimesX(gx1,gz1,2*P[1]+P[0],M);
				YplusisCtimesX(gx1,gz2,P[1],M);

				//UReflect(H,gz1,gz1);
				  Times   (H+1, l1,      gz1+1,  M-1);
				  AddTimes(H+1, l11,     gz1,    M-1);
				YplusisCtimesX(gx2,H+1,4*P[1],M-1);
				//UReflect(H,gz2,gz0);
				  Times   (H+1, l1,      gz2+1,  M-1);
				  AddTimes(H+1, l11,     gz0,    M-1);
				YplusisCtimesX(gx2,H+1,P[0],M-1);
/*
				LReflect(H,gz0,gz2); YplusisCtimesX(gx0+1,H,P[0],M-1);
				LReflect(H,gz1,gz1); YplusisCtimesX(gx0+1,H,4*P[1],M-1);

				YplusisCtimesX(gx1,gz0,P[1],M);
				YplusisCtimesX(gx1,gz1,2*P[1]+P[0],M);
				YplusisCtimesX(gx1,gz2,P[1],M);

				UReflect(H,gz1,gz1); YplusisCtimesX(gx2,H+1,4*P[1],M-1);
				UReflect(H,gz2,gz0); YplusisCtimesX(gx2,H+1,P[0],M-1);
*/

				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
			}
			break;
		case 2:
			if (lattice_type==hexagonal) {
				  Times   (H,   LABDA+0*M+2,    gz0,    M-2);
				  AddTimes(H,   LABDA_1+0*M+2,  gz4+2,  M-2);
				YplusisCtimesX(gx0+2,H,P[0],     M-2);
				  Times   (H,   LABDA+0*M+2,    gz1,    M-2);
				  AddTimes(H,   LABDA_1+0*M+2,  gz3+2,  M-2);
				YplusisCtimesX(gx0+2,H,2*P[1],   M-2);
				  Times   (H,   LABDA+0*M+2,    gz2,    M-2);
				  AddTimes(H,   LABDA_1+0*M+2,  gz2+2,  M-2);
				YplusisCtimesX(gx0+2,H,2*P[2],   M-2);
				  Times   (H,   LABDA+0*M+2,    gz3,    M-2);
				  AddTimes(H,   LABDA_1+0*M+2,  gz1+2,  M-2);
				YplusisCtimesX(gx0+2,H,2*P[3],   M-2);

				  Times   (H,   LABDA+1*M+1,    gz0,    M-1);
				  AddTimes(H,   LABDA_1+1*M+1,  gz4+1,  M-1);
				YplusisCtimesX(gx1+1,H,P[1],     M-1);
				  Times   (H,   LABDA+1*M+1,    gz1,    M-1);
				  AddTimes(H,   LABDA_1+1*M+1,  gz3+1,  M-1);
				YplusisCtimesX(gx1+1,H,P[0]+P[2],M-1);
				  Times   (H,   LABDA+1*M+1,    gz2,    M-1);
				  AddTimes(H,   LABDA_1+1*M+1,  gz2+1,  M-1);
				YplusisCtimesX(gx1+1,H,P[1]+P[3],M-1);
				  Times   (H,   LABDA+1*M+1,    gz3,    M-1);
				  AddTimes(H,   LABDA_1+1*M+1,  gz1+1,  M-1);
				YplusisCtimesX(gx1+1,H,P[2],     M-1);
				  Times   (H,   LABDA+1*M+1,    gz4,    M-1);
				  AddTimes(H,   LABDA_1+1*M+1,  gz0+1,  M-1);
				YplusisCtimesX(gx1+1,H,P[3],     M-1);

				YplusisCtimesX(gx2,  gz0,P[2],     M);
				YplusisCtimesX(gx2,  gz1,P[1]+P[3],M);
				YplusisCtimesX(gx2,  gz2,P[0],     M);
				YplusisCtimesX(gx2,  gz3,P[1]+P[3],M);
				YplusisCtimesX(gx2,  gz4,P[2],     M);

				  Times   (H+1, LABDA+3*M,      gz0+1,  M-1);
				  AddTimes(H+1, LABDA_1+3*M,    gz4,   M-1);
				YplusisCtimesX(gx3,H+1,P[3],     M-1);
				  Times   (H+1, LABDA+3*M,      gz1+1,  M-1);
				  AddTimes(H+1, LABDA_1+3*M,    gz3,    M-1);
				YplusisCtimesX(gx3,H+1,P[2],     M-1);
				  Times   (H+1, LABDA+3*M,      gz2+1,  M-1);
				  AddTimes(H+1, LABDA_1+3*M,    gz2,    M-1);
				YplusisCtimesX(gx3,H+1,P[1]+P[3],M-1);
				  Times   (H+1, LABDA+3*M,      gz3+1,  M-1);
				  AddTimes(H+1, LABDA_1+3*M,    gz1,    M-1);
				YplusisCtimesX(gx3,H+1,P[0]+P[2],M-1);
				  Times   (H+1, LABDA+3*M,      gz4+1,  M-1);
				  AddTimes(H+1, LABDA_1+3*M,    gz0,    M-1);
				YplusisCtimesX(gx3,H+1,P[1],     M-1);

				  Times   (H+2, LABDA+4*M,      gz1+2,  M-2);
				  AddTimes(H+2, LABDA_1+4*M,    gz3,    M-2);
				YplusisCtimesX(gx4,H+2,2*P[3],   M-2);
				  Times   (H+2, LABDA+4*M,      gz2+2,  M-2);
				  AddTimes(H+2, LABDA_1+4*M,    gz2,    M-2);
				YplusisCtimesX(gx4,H+2,2*P[2],   M-2);
				  Times   (H+2, LABDA+4*M,      gz3+2,  M-2);
				  AddTimes(H+2, LABDA_1+4*M,    gz1,    M-2);
				YplusisCtimesX(gx4,H+2,2*P[1],   M-2);
				  Times   (H+2, LABDA+4*M,      gz4+2,  M-2);
				  AddTimes(H+2, LABDA_1+4*M,    gz0,    M-2);
				YplusisCtimesX(gx4,H+2,P[0],     M-2);

				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);

			} else {
				cout <<"cubic lattice and fjc=2 Markov 2 not implemented " << endl;
			}
			break;
		default:
			cout <<"Error, FJC_choices > 5 not implemented " << endl;
			break;
	}
}

void LGrad1::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	Real *gs=G+FJC*M*(s_to);
	Real *gs_1=G+FJC*M*(s_from);
	Real *gz0=gs_1;
	Real *gz1=gs_1+M;
	Real *gz2=gs_1+2*M;
	Real *gz3=gs_1+3*M;
	Real *gz4=gs_1+4*M;

	Real *gx0=gs;
	Real *gx1=gs+M;
	Real *gx2=gs+2*M;
	Real *gx3=gs+3*M;
	Real *gx4=gs+4*M;
	Real *g =G1;

	Zero (gs,M*FJC);
	for (int k=0; k<(FJC-1)/2; k++) set_bounds(gs_1+k*M,gs_1+(FJC-k-1)*M);
	set_bounds(gs_1+(FJC-1)/2*M);
	switch (fjc) {
		case 1:
			if (lattice_type==hexagonal) {
				//LReflect(H,gz2,gz0);
				  Times   (H,   l_1 +1,  gz2,    M-1);
				  AddTimes(H,   l_11+1,  gz0+1,  M-1);
				YplusisCtimesX(gx1+1,H,P[1],M-1);
				YplusisCtimesX(gx2+1,H,P[0],M-1);

				YplusisCtimesX(gx0,gz1,2*P[1],M);
				YplusisCtimesX(gx1,gz1,P[0],M);
				YplusisCtimesX(gx2,gz1,2*P[1],M);

				//UReflect(H,gz0,gz2);
				  Times   (H+1, l1,      gz0+1,  M-1);
				  AddTimes(H+1, l11,     gz2,    M-1);
				YplusisCtimesX(gx0,H+1,P[0],M-1);
				YplusisCtimesX(gx1,H+1,P[1],M-1);
				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
			} else {
				//LReflect(H,gz2,gz0);
				  Times   (H,   l_1 +1,  gz2,    M-1);
				  AddTimes(H,   l_11+1,  gz0+1,  M-1);
				YplusisCtimesX(gx1+1,H,P[1],M-1);
				YplusisCtimesX(gx2+1,H,P[0],M-1);

				YplusisCtimesX(gx0,gz1,4*P[1],M);
				YplusisCtimesX(gx1,gz1,2*P[1]+P[0],M);
				YplusisCtimesX(gx2,gz1,4*P[1],M);

				//UReflect(H,gz0,gz2);
				  Times   (H+1, l1,      gz0+1,  M-1);
				  AddTimes(H+1, l11,     gz2,    M-1);
				YplusisCtimesX(gx0,H+1,P[0],M-1);
				YplusisCtimesX(gx1,H+1,P[1],M-1);

/*
			LReflect(H,gz2,gz0);
			YplusisCtimesX(gx1+1,H,P[1],M);
			YplusisCtimesX(gx2+1,H,P[0],M-1);

			YplusisCtimesX(gx0,gz1,4*P[1],M);
			YplusisCtimesX(gx1,gz1,2*P[1]+P[0],M);
			YplusisCtimesX(gx2,gz1,4*P[1],M);

			UReflect(H,gz0,gz2);
			YplusisCtimesX(gx0,H+1,P[0],M-1);
			YplusisCtimesX(gx1,H+1,P[1],M-1);
			for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
*/

				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);
			}
			break;
		case 2:
           if (lattice_type ==hexagonal) {
				  Times   (H,   LABDA+0*M+2,    gz4,    M-2);
				  AddTimes(H,   LABDA_1+0*M+2,  gz0+2,  M-2);
				YplusisCtimesX(gx1+2,H,   P[3],     M-2);
				YplusisCtimesX(gx2+2,H,   P[2],     M-2);
				YplusisCtimesX(gx3+2,H,   P[1],     M-2);
				YplusisCtimesX(gx4+2,H,   P[0],     M-2);

				  Times   (H,   LABDA+1*M+1,    gz3,    M-1);
				  AddTimes(H,   LABDA_1+1*M+1,  gz1+1,  M-1);
				YplusisCtimesX(gx0+1,H,   2*P[3],   M-1);
				YplusisCtimesX(gx1+1,H,   P[2],     M-1);
				YplusisCtimesX(gx2+1,H,   P[1]+P[3],M-1);
				YplusisCtimesX(gx3+1,H,   P[0]+P[2],M-1);
				YplusisCtimesX(gx4+1,H,   2*P[1],   M-1);

				YplusisCtimesX(gx0,  gz2,   2*P[2],   M);
				YplusisCtimesX(gx1,  gz2,   P[1]+P[3],M);
				YplusisCtimesX(gx2,  gz2,   P[0],     M);
				YplusisCtimesX(gx3,  gz2,   P[1]+P[3],M);
				YplusisCtimesX(gx4,  gz2,   2*P[2],   M);

				  Times   (H+1, LABDA+3*M,      gz1+1,  M-1);
				  AddTimes(H+1, LABDA_1+3*M,    gz3,    M-1);
				YplusisCtimesX(gx0,  H+1, 2*P[1],   M-1);
				YplusisCtimesX(gx1,  H+1, P[0]+P[2],M-1);
				YplusisCtimesX(gx2,  H+1, P[1]+P[3],M-1);
				YplusisCtimesX(gx3,  H+1, P[2],     M-1);
				YplusisCtimesX(gx4,  H+1, 2*P[3],   M-1);

				  Times   (H+2, LABDA+4*M,      gz0+2,  M-2);
				  AddTimes(H+2, LABDA_1+4*M,    gz4,    M-2);;
				YplusisCtimesX(gx0,  H+2, P[0],     M-2);
				YplusisCtimesX(gx1,  H+2, P[1],     M-2);
				YplusisCtimesX(gx2,  H+2, P[2],     M-2);
				YplusisCtimesX(gx3,  H+2, P[3],     M-2);

				for (int k=0; k<FJC; k++) Times(gs+k*M,gs+k*M,g,M);

			} else {
				cout <<" cubic lattice, fjc=2, Markov=2 not implemented " << endl;
			}
			break;
		default:
			cout <<" fjc > 2 not implemented " << endl;
			break;

	}
}


void LGrad1::propagate(Real *G, Real *G1, int s_from, int s_to,int M) {
if (debug) cout <<" propagate in LGrad1 " << endl;
	Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	int kk;
	int j;
	Zero(gs,M); set_bounds(gs_1);

	if (fjc==1) {
		AddTimes(gs,gs_1,lambda0,M);
		AddTimes(gs+1,gs_1,lambda_1+1,M-1);
		AddTimes(gs,gs_1+1,lambda1,M-1);
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


bool LGrad1::ReadRange(int* r, int* H_p, int &n_pos, bool &block, string range, int var_pos, string seg_name, string range_type) {
if (debug) cout <<"ReadRange in LGrad1 " << endl;
	bool success=true;
	vector<string>set;
	vector<string>coor;
	vector<string>xyz;
	In[0]->split(range,';',set);

	coor.clear();
	block=true;
	In[0]->split(set[0],',',coor);
	if (coor.size()!=1) {cout << "In mon " + seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes must come as a single coordinate 'x'" << endl; r[0]=0; success=false;}
	else r[0]=In[0]->Get_int(coor[0],-1) ;

	coor.clear(); In[0]->split(set[1],',',coor);

	if (coor.size()!=1) {cout << "In mon " + seg_name+ ", for 'pos 2', in '" + range_type + "' the coordinates must come as a single coordinate 'x'" << endl; r[3]=0; success=false;}
	else r[3]=In[0]->Get_int(coor[0],-1);
	if (r[0] > r[3]) {cout << "In mon " + seg_name+ ", for 'pos 1', the x-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}

	return success;
}

bool LGrad1::ReadRangeFile(string filename,int* H_p, int &n_pos, string seg_name, string range_type) {
if (debug) cout <<"ReadRangeFile in LGrad1 " << endl;
	if (fjc>1) {
		cout << "Rangefile is not implemented for FJC-choices >3; contact FL. " << endl;
		return false;
	}
	bool success=true;
	string content;
	vector<string> lines;
	vector<string> sub;
	vector<string> xyz;
	string Infilename=In[0]->name;
	In[0]->split(Infilename,'.',sub);

	int length;
	int length_xyz;
	int px,p_i,x;
	int i=0;
	if (!In[0]->ReadFile(sub[0].append(".").append(filename),content)) {
		success=false;
		return success;
	}

	In[0]->split(content,'#',lines);
	length = lines.size();
	if (length == MX) { //expect to read 'mask file';
		if (n_pos==0) {
			for (i = 0 ; i < length ; ++i) {
				if (In[0]->Get_int(lines[i],0)==1) n_pos++;
			}
			if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any elements." << endl;}
		} else {
			p_i=0;
			for (x=1; x<MX+1; x++) {
				if (In[0]->Get_int(lines[x-1],0)==1) {H_p[p_i]=x; p_i++;}
			}
		}
	} else { //expect to read x only
		px=0; i=0;
		if (n_pos==0) n_pos=length;
		else {
			while (i<length) {
				xyz.clear();
				In[0]->split(lines[i],',',xyz);
				length_xyz=xyz.size();
				if (length_xyz!=1) {
					cout << "In mon " + seg_name + " " +range_type+"_filename  the expected 'single coordinate' 'x' was not found. " << endl;  success = false;
				} else {
					px=In[0]->Get_int(xyz[0],0);
					if (px < 1 || px > MX) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the x-coordinate in "+range_type+"_filename out of bounds: 1.." << MX << endl; success =false;}
				}
				H_p[i]=px;
				i++;
			}
		}
	}
	return success;
}

bool LGrad1::FillMask(int* Mask, vector<int>px, vector<int>py, vector<int>pz, string filename) {
	bool success=true;
	bool readfile=false;
	int length=0;
	int length_px = px.size();

	vector<string> lines;
	int p;
	if (px.size()==0) {
		readfile=true;
		string content;
		success=In[0]->ReadFile(filename,content);
		if (success) {
			In[0]->split(content,'#',lines);
			length = lines.size();
		}
	}
	if (readfile) {
		if (MX!=length) {success=false; cout <<"inputfile for filling delta_range has not expected length in x-direction" << endl;
		} else {
			for (int x=1; x<MX+1; x++) Mask[x]=In[0]->Get_int(lines[x],-1);
		}
	} else  {
		for (int i=0; i<length_px; i++) {
			p=px[i]; if (p<1 || p>MX) {success=false; cout<<" x-value in delta_range out of bounds; " << endl; }
			else Mask[fjc-1+px[i]]=1;
		}
	}
	for (int i=0; i<M; i++) if (!(Mask[i]==0 || Mask[i]==1)) {success =false; cout <<"Delta_range does not contain '0' or '1' values. Check delta_inputfile values"<<endl; }
	return success;
}

bool LGrad1::CreateMASK(int* H_MASK, int* r, int* H_P, int n_pos, bool block) {
if (debug) cout <<"CreateMask for LGrad1 " + name << endl;
	bool success=true;
	H_Zero(H_MASK,M);
	if (block) {
		for (int x=r[0]; x<r[3]+1; x++) {
			H_MASK[x]=1;
		}
	} else {
		for (int i = 0; i<n_pos; i++) H_MASK[H_P[i]]=1;
	}
	return success;
}


Real LGrad1::ComputeTheta(Real* phi) {
	Real result=0; remove_bounds(phi);
	if (geometry !="planar") Dot(result,phi,L,M);
	else {if (fjc==1) Sum(result,phi,M); else  Dot(result,phi,L,M);}
	return result/fjc;
}

void LGrad1::UpdateEE(Real* EE, Real* psi, Real* E) {
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_M_bounds(psi);
	Zero(EE,M);
	Real Exmin,Explus;
	int x;
	int r;

	if (geometry=="cylindrical" ) {
		r=offset_first_layer*fjc;
		pf=pf*PIE;
		for (x=fjc; x<MX+fjc; x++) {
			r++;
			Exmin=psi[x]-psi[x-1];
			Exmin*=(r-1)*Exmin;
			Explus=psi[x]-psi[x+1];
			Explus*=(r)*Explus;
			EE[x]=pf*(Exmin+Explus)/L[x];
		}


/*
		pf=pf*PIE;
		r=offset_first_layer*fjc-1.5;
		if (offset_first_layer <=0) {
			r++;
			Explus = psi[fjc+1]-psi[fjc+2];
			Explus *=(r+0.5)*Explus;
			EE[fjc+1]=pf*Explus/L[fjc+1];
			x=fjc+2;
		} else {
			Explus=psi[fjc]-psi[fjc+1];
			Explus *=(r+0.5)*Explus;
			x=fjc+1;
		}

		for (; x<MX+fjc; x++) {
			r +=1.0;
			Exmin=Explus;
			Explus=(psi[x]-psi[x+1]);
			Explus *= (r+0.5)*Explus;
			EE[x]=pf*(Exmin+Explus)/L[x];
		}
*/

	}
	if (geometry=="spherical" ) {
		pf=pf*PIE*2/fjc;
		r=offset_first_layer*fjc +1.0;
		Explus=r*(psi[fjc]-psi[fjc+1]);
		Explus *=Explus;
		EE[fjc]=pf*Explus/(L[fjc]);
		for (x=fjc+1; x<MX+fjc; x++) {
			r +=1.0;
			Exmin=Explus;
			Explus=r*(psi[x]-psi[x+1]);
			Explus *=Explus;
			EE[x]=pf*(Exmin+Explus)/L[x];
		}

	}
}


void LGrad1::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, int* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	int x;
	Real a,b,c,a_,b_,c_;
	Real r;
	Real epsXplus, epsXmin;
	//set_M_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);

   if (!fixedPsi0) {
	if (geometry=="cylindrical") {
		C=C/PIE;
		r=offset_first_layer*fjc;
		epsXplus=r*(eps[fjc-1]+eps[fjc]);
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			r++;
			epsXmin=epsXplus;
			epsXplus=r*(eps[x]+eps[x+1]);
			a=b; b=c; c=psi[x+1];
			X[x]=(epsXmin*a + C*q[x]*L[x] + epsXplus*c)/(epsXmin+epsXplus);
		 }
	}
	if (geometry=="spherical") {
		C=C/(2.0*PIE)*fjc;
		r=offset_first_layer*fjc;
		epsXplus=r*r*(eps[fjc-1]+eps[fjc]);
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			epsXmin=epsXplus;
			r++;
			epsXplus=r*r*(eps[x]+eps[x+1]);
			a=b; b=c; c=psi[x+1];
			X[x]=(epsXmin*a + C*q[x]*L[x] + epsXplus*c)/(epsXmin+epsXplus);
		 }
	}
	//Cp(psi,X,M);
	YisAminB(g,g,X,M);
   } else { //fixedPsi0 is true
	a=0; b=psi[fjc-1]; c=psi[fjc];
	for (x=fjc; x<MX+fjc; x++) {
		a=b; b=c; c=psi[x+1];
		if (Mask[x] == 0) psi[x]=0.5*(a+c)+q[x]*C/eps[x];
	}

	if (geometry=="cylindrical") {
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			a=b; b=c; c=psi[x+1];
			if (Mask[x] == 0) psi[x]+=(c-a)/(2.0*(offset_first_layer*fjc+x-fjc+0.5))*fjc;
		}
	}
	if (geometry=="spherial") {
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			a=b; b=c; c=psi[x+1];
			if (Mask[x] == 0) psi[x]+=(c-a)/(offset_first_layer*fjc+x-fjc+0.5)*fjc;
		}
	}
	if (grad_epsilon) {
		a=0; b=psi[fjc-1]; c=psi[fjc];a_=0; b_=eps[fjc-1]; c_=eps[fjc];
		for (x=fjc; x<MX+fjc; x++) {//for all geometries
			a=b; b=c; c=psi[x+1]; a_=b_; b_=c_; c_=eps[x+1];
			if (Mask[x] == 0) {
				psi[x]+=0.25*(c_-a_)*(c-a)/eps[x]*fjc*fjc;
			}
		}
	}
	for (x=fjc; x<MX+fjc; x++)
	if (Mask[x] == 0) {
		g[x]-=psi[x];
	}
   }
}


void LGrad1::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, int* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
	int x;
	Real a,b,c,a_,b_,c_;

	Real C = -e*e/(eps0*k_BT*bond_length);
	a=0; b=psi[fjc-1]; c=psi[fjc];
	for (x=fjc; x<MX+fjc; x++) { //for all geometries
		a=b; b=c; c=psi[x+1];
		if (Mask[x] == 1) q[x] = -0.5*(a-2*b+c)*fjc*fjc*eps[x]/C;
	}

	if (geometry=="cylindrical") {
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			a=b; b=c; c=psi[x+1];
			if (Mask[x] == 1) q[x]-=(c-a)/(2.0*(offset_first_layer*fjc+x-fjc+0.5))*fjc*eps[x]/C;
		}
	}
	if (geometry=="spherial") {
		a=0; b=psi[fjc-1]; c=psi[fjc];
		for (x=fjc; x<MX+fjc; x++) {
			a=b; b=c; c=psi[x+1];
			if (Mask[x] == 1) q[x]-=(c-a)/(offset_first_layer*fjc+x-fjc+0.5)*fjc*eps[x]/C;
		}
	}
	if (grad_epsilon) {
		a=0; b=psi[fjc-1]; c=psi[fjc]; a_=0; b_=eps[fjc-1]; c_=eps[fjc];
		for (x=fjc; x<MX+fjc; x++) {//for all geometries
			a=b; b=c; c=psi[x+1]; a_=b_; b_=c_; c_=eps[x+1];
			if (Mask[x] == 1) q[x]-=0.25*(c_-a_)*(c-a)*fjc*fjc/C;
		}
	}
	for (x=fjc; x<MX+fjc; x++)
	if (Mask[x] == 1) {
		g[x]=-q[x];
	}

}

void LGrad1::remove_bounds(Real *X){
if (debug) cout <<"remove_bounds in LGrad1 " << endl;
	int k;
	if (fjc==1) {
		X[0]=0;
		X[MX+1]=0;
	} else {
		for (k=0; k<fjc; k++) {
			X[k]=0;
			X[MX+fjc+k]=0;
		}
	}
}

void LGrad1::set_bounds(Real* X, Real* Y){
if (debug) cout <<"set_bounds in LGrad1 " << endl;
	int k;
	if (fjc==1) {
		X[0]=Y[BX1];
		X[MX+1]=Y[BXM];
		Y[0]=X[BX1];
		Y[MX+1]=X[BXM];

	} else {
		for (k=0; k<fjc; k++) {
			X[k]=Y[B_X1[k]];
			X[MX+fjc+k]=Y[B_XM[k]];
			Y[k]=X[B_X1[k]];
			Y[MX+fjc+k]=X[B_XM[k]];
		}
	}
}

void LGrad1::set_bounds(Real* X){
if (debug) cout <<"set_bounds in LGrad1 " << endl;
	int k=0;
	if (fjc==1) {
		X[0]=X[BX1];
		X[MX+1]=X[BXM];
	} else {
		for (k=0; k<fjc; k++) {
			X[k]=X[B_X1[k]];
			X[MX+fjc+k]=X[B_XM[k]];
		}
	}
}

void LGrad1::set_M_bounds(Real* X){
if (debug) cout <<"set_M_bounds in LGrad1 " << endl; //set mirror bounds
	int k=0;
	if (fjc==1) {
		X[0]=X[1];
		X[MX+1]=X[MX];
	} else {
		for (k=0; k<fjc; k++) {
			X[k]=X[2*fjc-k-1];
			X[MX+fjc+k]=X[MX+fjc-k-1];
		}
	}
}


void LGrad1::remove_bounds(int *X){
if (debug) cout <<"remove_bounds in LGrad1 " << endl;
int k;
	if (fjc==1) {
		X[0]=0;
		X[MX+1]=0;
	} else {
		for (k=0; k<fjc; k++) {
			X[k]=0;
			X[MX+fjc+k]=0;
		}
	}
}

void LGrad1::set_bounds(int* X){
if (debug) cout <<"set_bounds in LGrad1 " << endl;
	int k=0;
	if (fjc==1) {
		X[0]=X[BX1];
		X[MX+1]=X[BXM];
	} else {
		for (k=0; k<fjc; k++) {
			X[k]=X[B_X1[k]];
			X[MX+fjc+k]=X[B_XM[k]];
		}
	}
}

Real LGrad1::ComputeGN(Real* G,int Markov, int M){
	Real GN=0;
	if (Markov==2) {
		GN=WeightedSum(G);
		for (int k=1; k<FJC-1; k++) {
			if (lattice_type == hexagonal) GN += 2.0*WeightedSum(G+k*M); else GN +=4.0*WeightedSum(G+k*M);
		}
		GN+=WeightedSum(G+(FJC-1)*M);
		if (lattice_type == hexagonal) GN /= 4.0*fjc; else GN /= 6.0;
	} else GN=WeightedSum(G);
	return GN;
}

void LGrad1::AddPhiS(Real* phi,Real* Gf,Real* Gb,int Markov, int M){
	if (debug) cout <<"AddPhiS_markov " << endl;
	if (Markov==2) {
		if (lattice_type ==hexagonal) {
			if (fjc==1) {
				YplusisCtimesAtimesB(phi,Gf,    Gb,    1.0/4.0,M);
				YplusisCtimesAtimesB(phi,Gf+1*M,Gb+1*M,2.0/4.0,M);
				YplusisCtimesAtimesB(phi,Gf+2*M,Gb+2*M,1.0/4.0,M);
			} else {
				YplusisCtimesAtimesB(phi,Gf,    Gb,    0.5/(FJC-1.0),M);
				for (int k=1; k<FJC-1; k++) YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,1.0/(FJC-1.0),M);
				YplusisCtimesAtimesB(phi,Gf+(FJC-1)*M,Gb+(FJC-1)*M,0.5/(FJC-1.0),M);
			}
		} else { //markov=2 cubic fjc=1
			YplusisCtimesAtimesB(phi,Gf,    Gb,    1.0/6.0,M);
			YplusisCtimesAtimesB(phi,Gf+1*M,Gb+1*M,4.0/6.0,M);
			YplusisCtimesAtimesB(phi,Gf+2*M,Gb+2*M,1.0/6.0,M);
		}
	} else {
		AddTimes(phi,Gf,Gb,M);
	}
}

void LGrad1::AddPhiS(Real* phi,Real* Gf,Real* Gb, Real degeneracy, int Markov, int M){
if (debug) cout <<"AddPhiS_degeneracy markov " << endl;
	if (Markov==2) {
		if (lattice_type ==hexagonal) {
			YplusisCtimesAtimesB(phi,Gf,Gb,degeneracy*0.5/(FJC-1.0),M);
			for (int k=1; k<FJC-1; k++) YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,degeneracy/(FJC-1.0),M);
			YplusisCtimesAtimesB(phi,Gf+(FJC-1)*M,Gb+(FJC-1)*M,degeneracy*0.5/(FJC-1.0),M);
		} else {
			YplusisCtimesAtimesB(phi,Gf,Gb,degeneracy/6.0,M);
			YplusisCtimesAtimesB(phi,Gf+1*M,Gb+1*M,degeneracy*4.0/6.0,M);
			YplusisCtimesAtimesB(phi,Gf+2*M,Gb+2*M,degeneracy/6.0,M);
		}
	} else {
		YplusisCtimesAtimesB(phi,Gf,Gb,degeneracy,M);
	}
}


void LGrad1::AddPhiS(Real* phi,Real* Gf,Real* Gb,Real* G1, Real norm, int Markov, int M){
if (debug) cout <<"AddPhiS_norm_markov " << endl;
	if (Markov==2) {
		if (lattice_type ==hexagonal) {
			Composition (phi,Gf,Gb,G1,norm*0.5/(FJC-1.0),M);
			for (int k=1; k<FJC-1; k++) Composition (phi,Gf+k*M,Gb+k*M,G1,norm/(FJC-1.0),M);
			Composition (phi,Gf+(FJC-1)*M,Gb+(FJC-1)*M,G1,norm*0.5/(FJC-1.0),M);
		} else {
			//if (fjc>1) cout<<"AddPhiS in cubic lattice, Markov=2 FJC_choice >3 not implemented" << endl;
			Composition (phi,Gf,Gb,G1,norm/6.0,M);
			Composition (phi,Gf+1*M,Gb+1*M,G1,norm*4.0/6.0,M);
			Composition (phi,Gf+2*M,Gb+2*M,G1,norm/6.0,M);
		}
	} else {
		Composition (phi,Gf,Gb,G1,norm,M);
	}
}


void LGrad1::Initiate(Real* G,Real* Gz,int Markov, int M){
	if (Markov==2) {
		for (int k=0; k<FJC; k++) Cp(G+k*M,Gz,M);
	} else {
		Cp(G,Gz,M);
	}
}

void LGrad1::Terminate(Real* Gz ,Real* G, int Markov, int M){
if (!debug) cout <<"LGrad1::Terminite " << endl;
	Real one=1.0;
	if (Markov==2) {
		Zero(Gz,M);
		if (lattice_type == simple_cubic) {
			Add(Gz,G+M,M); Norm(Gz,4.0*one,M);
			Add(Gz,G,M); Add(Gz,G+2*M,M);
			Norm(Gz,1.0/6.0*one,M);
		} else {
			Add(Gz,G+M,M); Norm(Gz,2.0*one,M);
			Add(Gz,G,M); Add(Gz,G+2*M,M);
			Norm(Gz,1.0/4.0*one,M);

		}
	} else {
		Cp(Gz,G,M);
	}
}

bool LGrad1:: PutMask(int* MASK,vector<int>px,vector<int>py,vector<int>pz,int R){
	bool success=false;
	cout <<"PutMask does not make sence in 1 gradient system " << endl;
	return success;
}

Real LGrad1::DphiDt(Real* g, Real* B_phitot, Real* phiA, Real* phiB, Real* alphaA, Real* alphaB,Real B_A, Real B_B) {
	//cout <<"LGrad1 : DphiDt not implemented yet " << endl;

	Real AverageJ=0;
	//Real Jplus,Jmin;
	Real a,b,c,Ma,Mb,Mc;


	//g[1]=phiA[0]/phiA[1]-1.0;
	b=phiA[1]*phiB[1]*B_B/B_phitot[1];
	c=phiA[2]*phiB[2]*B_B/B_phitot[2];
	Mb=alphaA[1]-alphaB[1];
	Mc=alphaA[2]-alphaB[2];
	for (int z=2; z<M-2; z++) {
		a=b; b=c; c=phiA[z+1]*phiB[z+1]*B_B/B_phitot[z+1];
		Ma=Mb; Mb=Mc; Mc=alphaA[z+1]-alphaB[z+1];
		//Jmin =(a+b)*(Mb-Ma);
		//if (Jmin <0) Jmin=-log(-Jmin); else Jmin=log(Jmin);
		//Jplus = (b+c)*(Mc-Mb);
		//if (Jplus <0) Jplus=-log(-Jplus); else Jplus=log(Jplus);
		//g[z]=g[z]+ (Jmin-Jplus);
		//g[z]=g[z]+ (Jmin-Jplus)/abs(Jmin+Jplus);
		g[z] = g[z]  + (a+b)*(Mb-Ma)*lambda_1[z]-(b+c)*(Mc-Mb)*lambda1[z];///L[z];

		AverageJ+=lambda_1[z]*L[z]*(a+b)*(Mb-Ma);
	}
	//g[M-2]=phiA[M-1]/phiA[M-2]-1.0;

	return -B_A*AverageJ/(2*(M-4)*lambda);

}
