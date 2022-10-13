#include <iostream>
#include <string>
#include "lattice.h"
#include "LGrad2.h"

LGrad2::LGrad2(vector<Input*> In_,string name_): Lattice(In_,name_) {}

LGrad2::~LGrad2() {
if (debug) cout <<"LGrad2 destructor " << endl;
}


void LGrad2:: ComputeLambdas() {
	Real r, VL, LS;
	Real rlow, rhigh;


	//if (fcc_sites){
	//	for (int x=1; x<MX+1; x++)
	//	for (int y=1; y<MY+1; y++) {
	//		r=offset_first_layer + 1.0*x;
	//		L[P(x,y)]=PIE*(pow(r,2)-pow(r-1,2));
	//		fcc_lambda1[P(x,y)]=2.0*PIE*r/L[P(x,y)]/3.0;
	//		fcc_lambda_1[P(x,y)]=2.0*PIE*(r-1)/L[P(x,y)]/3.0;
	//		fcc_lambda0[P(x,y)]=1.0-2.0/3.0;
	//	}
	//
	//}
	if (fjc==1) {
		for (int x=1; x<MX+1; x++)
		for (int y=1; y<MY+1; y++) {
			r=offset_first_layer + 1.0*x;
			L[P(x,y)]=PIE*(pow(r,2)-pow(r-1,2));
			lambda1[P(x,y)]=2.0*PIE*r/L[P(x,y)]*lambda;
			lambda_1[P(x,y)]=2.0*PIE*(r-1)/L[P(x,y)]*lambda;
			lambda0[P(x,y)]=1.0-2.0*lambda;
			if (fcc_sites) {
				fcc_lambda1[P(x,y)]=2.0*PIE*r/L[P(x,y)]/3.0;
				fcc_lambda_1[P(x,y)]=2.0*PIE*(r-1)/L[P(x,y)]/3.0;
				fcc_lambda0[P(x,y)]=1.0-2.0/3.0;
			}
		}
		if (Markov ==2) {
			for (int i=0; i<M; i++) {
				l1[i]=lambda1[i]/lambda; l11[i]=1.0-l1[i];
				l_1[i]=lambda_1[i]/lambda; l_11[i]=1.0-l_1[i];
			}
		}

	}
	if (fjc>1) {
		for (int y=fjc; y<MY+fjc; y++) {
			for (int x = fjc; x < MX+fjc; x++) {
				r = offset_first_layer+1.0*(x-fjc+1.0)/fjc;
				rlow = r - 0.5;
				rhigh = r + 0.5;
				L[P(x,y)] = PIE * (2.0 * r) / fjc;
				VL = L[P(x,y)] / PIE * fjc;
				if ((rlow - r) * 2 + r > 0.0) {
					LAMBDA[P(x,y)] += 1.0/(1.0*FJC-1.0)*rlow/VL;
				}
				if ((rhigh - r) * 2 + r < 1.0*MX/fjc) {
					LAMBDA[P(x,y)+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
				} else {
					if (2*rhigh-r-1.0*MX/fjc > -0.001 && 2 * rhigh-r-1.0*MX/fjc < 0.001) {
						LAMBDA[P(x,y)+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
					}
					for (int j = 1; j <= fjc; j++) {
						if (2*rhigh-r-1.0*MX/fjc > 0.99*j/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*j/fjc) {
							LAMBDA[P(x,y)+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*(rhigh-1.0*j/fjc)/VL;
						}
					}								}
				for (int j = 1; j < fjc; j++) {
					rlow += 0.5/(fjc);
					rhigh -= 0.5/(fjc);
					if ((rlow-r)*2+r > 0.0)
					LAMBDA[P(x,y)+j*M] += 1.0/(1.0*FJC-1.0)*2.0*rlow/VL;
					if ((rhigh-r)*2+r < offset_first_layer+1.0*MX/fjc)
					LAMBDA[P(x,y)+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
					else {
						if (2 * rhigh-r-1.0*MX/fjc > -0.001 && 2*rhigh-r-1.0*MX/fjc < 0.001) {
							LAMBDA[P(x,y)+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
						}
						for (int k = 1; k <= fjc; k++) {
							if (2 * rhigh-r-1.0*MX/fjc > 0.99*k/fjc && 2*rhigh-r-1.0*MX/fjc<1.01*k/fjc) {
								LAMBDA[P(x,y) + (FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*(rhigh-1.0*k/fjc)/VL;
							}
						}
					}
				}
				LS = 0;
				for (int j = 0; j < FJC; j++)
				LS += LAMBDA[P(x,y)+j*M];
				LAMBDA[P(x,y)+(FJC/2)*M] += 1.0 - LS;
			}
		}
	}
	for (int k=0; k<FJC; k++) for (int y=1-fjc; y<MY+fjc; y++) for (int x=MX+1; x< MX+fjc; x++)
		LAMBDA[P(x,y)+k*M]=LAMBDA[P(2*MX-x+1,y)+(FJC-k-1)*M];
}

bool LGrad2::PutM() {
if (debug) cout << "PutM in LGrad2 " << endl;
	bool success=true;

	if (geometry=="cylindrical")
		volume = MY*PIE*(pow(MX+offset_first_layer,2)-pow(offset_first_layer,2));
	else volume = MX*MY;
	JX=MY+2*fjc; JY=1; JZ=0; M=(MX+2*fjc)*(MY+2*fjc);

	Accesible_volume=volume;
	return success;
}

void LGrad2::TimesL(Real* X){
if (debug) cout << "TimesL in LGrad2 " << endl;
	if (geometry!="planar") Times(X,X,L,M);
}

void LGrad2::DivL(Real* X){
if (debug) cout << "DivL in LGrad2 " << endl;
	if (geometry!="planar") Div(X,L,M);
}

Real LGrad2:: Moment(Real* X,Real Xb, int n) {
if (debug) cout << "Moment in LGrad2 " << endl;
	Real Result=0;
	int x,y;
	int zz;
	Real Nz;
	if (fjc==1) {
		zz=0;
		for (y=1; y<=MY ; y++) {
			Nz=0;
			for (x=1; x<=MX; x++) {
				if (X[P(x,y)]>0) Nz+=(X[P(x,y)]-Xb)*L[P(x,y)];
			}
			if (Nz>0) zz++;
			if (zz>0) Result+= pow(zz,n)*Nz;
		}
	} else {
		//cout <<"Moment analysis not yet implemented in LGrad" << endl;
	}
	return Result/fjc;
}

Real LGrad2::WeightedSum(Real* X){
if (debug) cout << "weighted sum in LGrad2 " << endl;
	Real sum{0};
	remove_bounds(X);
	if (geometry=="planar") {
		Sum(sum,X,M); sum = sum/(fjc*fjc);
	} else	{
		Dot(sum,X,L,M);
	}
	return sum;
}

void LGrad2::vtk(string filename, Real* X, string id,bool writebounds) {
if (debug) cout << "vtk in LGrad2 " << endl;
	FILE *fp;
	int i;
	fp = fopen(filename.c_str(),"w+");
	fprintf(fp,"# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\n",MY,MX,1);
	if (writebounds) {
		fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",(MX+2*fjc)*(MY+2*fjc));
	} else {
		fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",MX*MY);
	}
	fprintf(fp,"SCALARS %s double\nLOOKUP_TABLE default\n",id.c_str());
#ifdef LongReal
	if (writebounds) for (i=0; i<M; i++) fprintf(fp,"%Lf\n",X[i]);
	for (int x=fjc; x<MX+fjc; x++)
	for (int y=fjc; y<MY+fjc; y++)
	fprintf(fp,"%Lf\n",X[P(x,y)]);
#else
	if (writebounds) for (i=0; i<M; i++) fprintf(fp,"%f\n",X[i]);
	for (int x=fjc; x<MX+fjc; x++)
	for (int y=fjc; y<MY+fjc; y++)
	fprintf(fp,"%f\n",X[P(x,y)]);
#endif
	fclose(fp);
}

void LGrad2::PutProfiles(FILE* pf,vector<Real*> X,bool writebounds){
if (debug) cout <<"PutProfiles in LGrad2 " << endl;
	Real one=1.0;
	int x,y,i;
	int length=X.size();
	int a;
	if (writebounds) a=0; else a = fjc;
	for (x=a; x<MX+2*fjc-a; x++)
	for (y=a; y<MY+2*fjc-a; y++){
#ifdef LongReal
		fprintf(pf,"%Le\t%Le\t",offset_first_layer/fjc+one*(x-fjc+1)/fjc-0.5/fjc,one*(y-fjc+1)/fjc-0.5/fjc);
		for (i=0; i<length; i++) fprintf(pf,"%.20Le\t",X[i][P(x,y)]);
		fprintf(pf,"\n");
#else
		fprintf(pf,"%e\t%e\t",offset_first_layer/fjc+one*(x-fjc+1)/fjc-0.5/fjc,one*(y-fjc+1)/fjc-0.5/fjc);
		for (i=0; i<length; i++) fprintf(pf,"%.20e\t",X[i][P(x,y)]);
		fprintf(pf,"\n");
#endif
	}
}


void LGrad2::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
if (debug) cout <<" Side in LGrad2 " << endl;
	if (ignore_sites) {
		Cp(X_side,X,M); return;
	}
	Zero(X_side,M);//set_bounds(X);

	if (fcc_sites) {
		YplusisCtimesX(X_side,X,1.0/3.0,M);
		AddTimes(X_side+JX,X,fcc_lambda_1+JX,M-JX);
		AddTimes(X_side,X+JX,fcc_lambda1,    M-JX);
		YplusisCtimesX(X_side+1,X,1.0/3.0,M-1);
		YplusisCtimesX(X_side,X+1,1.0/3.0,M-1);
		AddTimes(X_side+JX+1,X,fcc_lambda_1+JX+1,M-JX-1);
		AddTimes(X_side+JX,X+1,fcc_lambda_1+JX,  M-JX-1);
		AddTimes(X_side+1,X+JX,fcc_lambda1+1,M-JX-1);
		AddTimes(X_side,X+JX+1,fcc_lambda1,  M-JX-1);
		Norm(X_side,1.0/3.0,M);

	} else {
		switch (fjc) {
			case 1:
				if (lattice_type ==simple_cubic) {
					YplusisCtimesX(X_side,X,4.0/6.0,M);
					AddTimes(X_side+JX,X,lambda_1+JX,M-JX);
					AddTimes(X_side,X+JX,lambda1,    M-JX);
					YplusisCtimesX(X_side+1,X,1.0/6.0,M-1);
					YplusisCtimesX(X_side,X+1,1.0/6.0,M-1);
					Norm(X_side,4.0,M);
					AddTimes(X_side+JX+1,X,lambda_1+JX+1,M-JX-1);
					AddTimes(X_side+JX,X+1,lambda_1+JX,  M-JX-1);
					AddTimes(X_side+1,X+JX,lambda1+1,M-JX-1);
					AddTimes(X_side,X+JX+1,lambda1,  M-JX-1);
					Norm(X_side,1.0/6.0,M);
				} else {
					YplusisCtimesX(X_side,X,2.0/4.0,M);
					AddTimes(X_side+JX,X,lambda_1+JX,M-JX);
					AddTimes(X_side,X+JX,lambda1,    M-JX);
					YplusisCtimesX(X_side+1,X,1.0/4.0,M-1);
					YplusisCtimesX(X_side,X+1,1.0/4.0,M-1);
					Norm(X_side,2.0,M);
					AddTimes(X_side+JX+1,X,lambda_1+JX+1,M-JX-1);
					AddTimes(X_side+JX,X+1,lambda_1+JX,  M-JX-1);
					AddTimes(X_side+1,X+JX,lambda1+1,    M-JX-1);
					AddTimes(X_side,X+JX+1,lambda1,      M-JX-1);
					Norm(X_side,3.0/12.0,M);
				}
				break;
			case 2:
				AddTimes(X_side,     X,         LAMBDA+2*M,     M);

				AddTimes(X_side+1,   X,         LAMBDA+2*M+1,   M-1);
				AddTimes(X_side,     X+1,       LAMBDA+2*M,     M-1);
				AddTimes(X_side+JX,  X,         LAMBDA+M+JX,    M-JX);
				AddTimes(X_side,     X+JX,      LAMBDA+3*M,     M-JX);

				AddTimes(X_side+JX+1,X,     LAMBDA+M+JX+1,M-JX-1);
				AddTimes(X_side+JX,  X+1,   LAMBDA+M+JX+1,M-JX-1);
				AddTimes(X_side,     X+JX+1,LAMBDA+3*M,   M-JX-1);
				AddTimes(X_side+1,   X+JX,  LAMBDA+3*M,   M-JX-1);

				AddTimes(X_side+2*JX,X,     LAMBDA+2*JX,    M-2*JX);
				AddTimes(X_side,     X+2*JX,LAMBDA+4*M,     M-2*JX);

				AddTimes(X_side+2*JX,  X+1,     LAMBDA+2*JX,  M-2*JX-1);
				AddTimes(X_side+2*JX+1,X,       LAMBDA+2*JX+1,M-2*JX-1);
				AddTimes(X_side+1,     X+2*JX,  LAMBDA+4*M+1, M-2*JX-1);
				AddTimes(X_side,       X+2*JX+1,LAMBDA+4*M,   M-2*JX-1);

				Norm(X_side,2.0,M);

				AddTimes(X_side+2,   X,     LAMBDA+2*M+2,   M-2);
				AddTimes(X_side,     X+2,   LAMBDA+2*M,     M-2);

				AddTimes(X_side+JX+2, X,     LAMBDA+M+JX+2,  M-JX-2);
				AddTimes(X_side+JX,   X+2,   LAMBDA+M+JX,    M-JX-2);
				AddTimes(X_side,      X+JX+2,LAMBDA+3*M,     M-JX-2);
				AddTimes(X_side+2,    X+JX,  LAMBDA+3*M+2,   M-JX-2);


				AddTimes(X_side+2*JX+2, X,       LAMBDA+2*JX+2,M-2*JX-2);
				AddTimes(X_side,        X+2*JX+2,LAMBDA+4*M,   M-2*JX-2);
				AddTimes(X_side+2*JX,   X+2,     LAMBDA+2*JX,  M-2*JX-2);
				AddTimes(X_side+2,      X+2*JX,  LAMBDA+4*M+2, M-2*JX-2);
				Norm(X_side,1.0/8.0,M);
				break;
			case 3:

				AddTimes(X_side,        X,        LAMBDA+3*M,         M);

				AddTimes(X_side+1,      X,        LAMBDA+3*M+1,       M-1);
				AddTimes(X_side,        X+1,      LAMBDA+3*M,         M-1);
				AddTimes(X_side+JX,     X,        LAMBDA+2*M+JX,      M-JX);
				AddTimes(X_side,        X+JX,     LAMBDA+4*M,         M-JX);

				AddTimes(X_side+JX+1,   X,        LAMBDA+2*M+JX+1,    M-JX-1);
				AddTimes(X_side,        X+JX+1,   LAMBDA+4*M,         M-JX-1);
				AddTimes(X_side+1,      X+JX,     LAMBDA+4*M+1,       M-JX-1);
				AddTimes(X_side+JX,     X+1,      LAMBDA+2*M+JX,      M-JX-1);

				AddTimes(X_side+2*JX,   X,        LAMBDA+M+2*JX,      M-2*JX);
				AddTimes(X_side,        X+2*JX,   LAMBDA+5*M,         M-2*JX);
				AddTimes(X_side+2,      X,        LAMBDA+3*M+2,       M-2);
				AddTimes(X_side,        X+2,      LAMBDA+3*M,         M-2);

				AddTimes(X_side+2*JX+1, X,        LAMBDA+M+2*JX+1,    M-2*JX-1);
				AddTimes(X_side,        X+2*JX+1, LAMBDA+5*M,         M-2*JX-1);
				AddTimes(X_side+1,      X+2*JX,   LAMBDA+5*M+1,       M-2*JX-1);
				AddTimes(X_side+2*JX,   X+1,      LAMBDA+M+2*JX,      M-2*JX-1);

				AddTimes(X_side+JX+2,   X,        LAMBDA+2*M+JX+2,    M-JX-2);
				AddTimes(X_side,        X+JX+2,   LAMBDA+4*M,         M-JX-2);
				AddTimes(X_side+2,      X+JX,     LAMBDA+4*M+2,       M-JX-2);
				AddTimes(X_side+JX,     X+2,      LAMBDA+2*M+JX,      M-JX-2);

				AddTimes(X_side+2*JX+2, X,        LAMBDA+1*M+2*JX+2,  M-2*JX-2);
				AddTimes(X_side,        X+2*JX+2, LAMBDA+5*M,         M-2*JX-2);
				AddTimes(X_side+2,      X+2*JX,   LAMBDA+5*M+2,       M-2*JX-2);
				AddTimes(X_side+2*JX,   X+2,      LAMBDA+1*M+2*JX,    M-2*JX-2);

				AddTimes(X_side+3*JX,   X,        LAMBDA+3*JX,      M-3*JX);
				AddTimes(X_side,        X+3*JX,   LAMBDA+6*M,       M-3*JX);

				AddTimes(X_side+3*JX,   X+1,      LAMBDA+3*JX,      M-3*JX-1);
				AddTimes(X_side+3*JX+1, X,        LAMBDA+3*JX+1,    M-3*JX-1);
				AddTimes(X_side+1,      X+3*JX,   LAMBDA+6*M+1,     M-3*JX-1);
				AddTimes(X_side,        X+3*JX+1, LAMBDA+6*M,       M-3*JX-1);

				AddTimes(X_side+3*JX,   X+2,      LAMBDA+3*JX,      M-3*JX-2);
				AddTimes(X_side+3*JX+2, X,        LAMBDA+3*JX+2,    M-3*JX-2);
				AddTimes(X_side+2,      X+3*JX,   LAMBDA+6*M+2,     M-3*JX-2);
				AddTimes(X_side,        X+3*JX+2, LAMBDA+6*M,       M-3*JX-2);

				Norm(X_side,2.0,M);

				AddTimes(X_side+3,      X,        LAMBDA+3*M+3,     M-3);
				AddTimes(X_side,        X+3,      LAMBDA+3*M,       M-3);

				AddTimes(X_side+JX+3,   X,        LAMBDA+2*M+JX+3,  M-JX-3);
				AddTimes(X_side+JX,     X+3,      LAMBDA+2*M+JX,    M-JX-3);
				AddTimes(X_side,        X+JX+3,   LAMBDA+4*M,       M-JX-3);
				AddTimes(X_side+3,      X+JX,     LAMBDA+4*M+3,     M-JX-3);

				AddTimes(X_side+2*JX+3, X,        LAMBDA+1*M+2*JX+3,M-2*JX-3);
				AddTimes(X_side+2*JX,   X+3,      LAMBDA+1*M+2*JX,  M-2*JX-3);
				AddTimes(X_side,        X+2*JX+3, LAMBDA+5*M,       M-2*JX-3);
				AddTimes(X_side+3,      X+2*JX,   LAMBDA+5*M+3,     M-2*JX-3);

				AddTimes(X_side+3*JX+3, X,        LAMBDA+3*JX+3,    M-3*JX-3);
				AddTimes(X_side,        X+3*JX+3, LAMBDA+6*M,       M-3*JX-3);
				AddTimes(X_side+3*JX,   X+3,      LAMBDA+3*JX,      M-3*JX-3);
				AddTimes(X_side+3,      X+3*JX,   LAMBDA+6*M+3,     M-3*JX-3);

				Norm(X_side,1.0/10.0,M);
				break;
			default:
				break;
		}
	}
}

void LGrad2::LReflect(Real *H, Real *P, Real *Q) {
	Times(H,l_1+JX,P,M-JX);
	AddTimes(H,l_11+JX,Q+JX,M-JX);
}

void LGrad2::UReflect(Real *H, Real *P, Real *Q) {
	Times (H+JX,l1,P+JX,M-JX);
	AddTimes(H+JX,l11,Q,M-JX);
}


void LGrad2::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	if (lattice_type == hexagonal) {
		Real *gs=G+M*7*s_to;
		Real *gs_1=G+M*7*s_from;

		Real *gz0=gs_1;
		Real *gz1=gs_1+M;
		Real *gz2=gs_1+2*M;
		Real *gz3=gs_1+3*M;
		Real *gz4=gs_1+4*M;
		Real *gz5=gs_1+5*M;
		Real *gz6=gs_1+6*M;

		Real *gx0=gs;
		Real *gx1=gs+M;
		Real *gx2=gs+2*M;
		Real *gx3=gs+3*M;
		Real *gx4=gs+4*M;
		Real *gx5=gs+5*M;
		Real *gx6=gs+6*M;
		Real *g=G1;


		Zero(gs,7*M);
		remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);
		set_bounds_x(gz0,gz6,0);set_bounds_x(gz1,gz5,0);set_bounds_x(gz2,gz4,0);set_bounds_x(gz3,0);

		LReflect(H,gz0,gz6); YplusisCtimesX(gx0+JX,H,2*P[0],M-JX);
		LReflect(H,gz1,gz5); YplusisCtimesX(gx0+JX,H,  P[0],M-JX);
		LReflect(H,gz2,gz4); YplusisCtimesX(gx0+JX,H,2*P[1],M-JX);
		LReflect(H,gz3,gz3); YplusisCtimesX(gx0+JX,H,2*P[1],M-JX);
		LReflect(H,gz4,gz2); YplusisCtimesX(gx0+JX,H,2*P[1],M-JX);

		UReflect(H,gz2,gz4); YplusisCtimesX(gx6,H+JX,2*P[1],M-JX);
		UReflect(H,gz3,gz3); YplusisCtimesX(gx6,H+JX,2*P[1],M-JX);
		UReflect(H,gz4,gz2); YplusisCtimesX(gx6,H+JX,2*P[1],M-JX);
		UReflect(H,gz5,gz1); YplusisCtimesX(gx6,H+JX,  P[0],M-JX);
		UReflect(H,gz6,gz0); YplusisCtimesX(gx6,H+JX,2*P[0],M-JX);

		remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);
		set_bounds_y(gz0,gz6,0);set_bounds_y(gz1,gz5,0);set_bounds_y(gz2,gz4,0);set_bounds_y(gz3,0);

		YplusisCtimesX(gx2+JY,gz0,2*P[1],M-JY);
		YplusisCtimesX(gx2+JY,gz1,  P[1],M-JY);
		YplusisCtimesX(gx2+JY,gz2,2*P[0],M-JY);
		YplusisCtimesX(gx2+JY,gz3,  P[0],M-JY);
		YplusisCtimesX(gx2+JY,gz5,  P[1],M-JY);
		YplusisCtimesX(gx2+JY,gz6,2*P[1],M-JY);

		YplusisCtimesX(gx4,gz0+JY,2*P[1],M-JY);
		YplusisCtimesX(gx4,gz1+JY,  P[1],M-JY);
		YplusisCtimesX(gx4,gz3+JY,  P[0],M-JY);
		YplusisCtimesX(gx4,gz4+JY,2*P[0],M-JY);
		YplusisCtimesX(gx4,gz5+JY,  P[1],M-JY);
		YplusisCtimesX(gx4,gz6+JY,2*P[1],M-JY);

		remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);//remove_bounds(gz5);remove_bounds(gz6);

		YplusisCtimesX(gx3,gz0,2*P[1],M);
		YplusisCtimesX(gx3,gz1,P[1],M);
		YplusisCtimesX(gx3,gz2,P[0],M);
		YplusisCtimesX(gx3,gz3,P[0],M);
		YplusisCtimesX(gx3,gz4,P[0],M);
		YplusisCtimesX(gx3,gz5,P[1],M);
		YplusisCtimesX(gx3,gz6,2*P[1],M);

		set_bounds_y(gz0,-1); set_bounds_y(gz1,-1);set_bounds_y(gz2,-1);set_bounds_y(gz4,-1);set_bounds_y(gz3,-1);
		LReflect(H,gz0,gz6); YplusisCtimesX(gx1+JX,H+JY,2*P[0],M-JX-JY);
		LReflect(H,gz1,gz5); YplusisCtimesX(gx1+JX,H+JY,  P[0],M-JX-JY);
		LReflect(H,gz2,gz4); YplusisCtimesX(gx1+JX,H+JY,2*P[1],M-JX-JY);
		LReflect(H,gz3,gz3); YplusisCtimesX(gx1+JX,H+JY,2*P[1],M-JX-JY);
		LReflect(H,gz4,gz2); YplusisCtimesX(gx1+JX,H+JY,2*P[1],M-JX-JY);

		remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);
		set_bounds_x(gz6,-1);set_bounds_x(gz5,-1);set_bounds_x(gz2,-1);set_bounds_x(gz4,-1);set_bounds_x(gz3,-1);
		UReflect(H,gz2,gz4); YplusisCtimesX(gx5+JY,H+JX,2*P[1],M-JX-JY);
		UReflect(H,gz3,gz3); YplusisCtimesX(gx5+JY,H+JX,2*P[1],M-JX-JY);
		UReflect(H,gz4,gz2); YplusisCtimesX(gx5+JY,H+JX,2*P[1],M-JX-JY);
		UReflect(H,gz5,gz2); YplusisCtimesX(gx5+JY,H+JX,  P[0],M-JX-JY);
		UReflect(H,gz6,gz0); YplusisCtimesX(gx5+JY,H+JX,2*P[0],M-JX-JY);

		for (int k=0; k<7; k++) Times(gs+k*M,gs+k*M,g,M);

	} else {
		Real *gs=G+M*5*s_to;
		Real *gs_1=G+M*5*s_from;
		Real *gz0=gs_1;
		Real *gz1=gs_1+M;
		Real *gz2=gs_1+2*M;
		Real *gz3=gs_1+3*M;
		Real *gz4=gs_1+4*M;
		set_bounds_x(gz0,gz4,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0); set_bounds_x(gz3,0);
		set_bounds_y(gz1,gz3,0); set_bounds_y(gz0,0); set_bounds_y(gz2,0); set_bounds_y(gz4,0);
		Real *gx0=gs;
		Real *gx1=gs+M;
		Real *gx2=gs+2*M;
		Real *gx3=gs+3*M;
		Real *gx4=gs+4*M;
		Real *g=G1;

		Zero(gs,5*M);
		LReflect(H,gz0,gz4); YplusisCtimesX(gx0+JX,H,P[0],M-JX);
		LReflect(H,gz1,gz3); YplusisCtimesX(gx0+JX,H,P[1],M-JX);
		LReflect(H,gz2,gz2); YplusisCtimesX(gx0+JX,H,2*P[1],M-JX);
		LReflect(H,gz3,gz1); YplusisCtimesX(gx0+JX,H,P[1],M-JX);

		YplusisCtimesX(gx1+JY,gz0,P[1],M-JY);
		YplusisCtimesX(gx1+JY,gz1,P[0],M-JY);
		YplusisCtimesX(gx1+JY,gz2,2*P[1],M-JY);
		YplusisCtimesX(gx1+JY,gz4,P[1],M-JY);

		YplusisCtimesX(gx2,gz0,P[1],M);
		YplusisCtimesX(gx2,gz1,P[1],M);
		YplusisCtimesX(gx2,gz2,P[0],M);
		YplusisCtimesX(gx2,gz3,P[1],M);
		YplusisCtimesX(gx2,gz4,P[1],M);

		YplusisCtimesX(gx3,gz0+JY,P[1],M-JY);
		YplusisCtimesX(gx3,gz2+JY,2*P[1],M-JY);
		YplusisCtimesX(gx3,gz3+JY,P[0],M-JY);
		YplusisCtimesX(gx3,gz4+JY,P[1],M-JY);

		UReflect(H,gz1,gz3); YplusisCtimesX(gx4,H+JX,P[1],M-JX);
		UReflect(H,gz2,gz2); YplusisCtimesX(gx4,H+JX,2*P[1],M-JX);
		UReflect(H,gz3,gz1); YplusisCtimesX(gx4,H+JX,P[1],M-JX);
		UReflect(H,gz4,gz0); YplusisCtimesX(gx4,H+JX,P[0],M-JX);

		for (int k=0; k<5; k++) Times(gs+k*M,gs+k*M,g,M);
	}
}

void LGrad2::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	if (lattice_type==hexagonal) {
		Real *gs=G+M*7*s_to;
		Real *gs_1=G+M*7*s_from;

		Real *gz0=gs_1;
		Real *gz1=gs_1+M;
		Real *gz2=gs_1+2*M;
		Real *gz3=gs_1+3*M;
		Real *gz4=gs_1+4*M;
		Real *gz5=gs_1+5*M;
		Real *gz6=gs_1+6*M;

		Real *gx0=gs;
		Real *gx1=gs+M;
		Real *gx2=gs+2*M;
		Real *gx3=gs+3*M;
		Real *gx4=gs+4*M;
		Real *gx5=gs+5*M;
		Real *gx6=gs+6*M;
		Real *g=G1;

		Zero(gs,7*M);
		remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);

		set_bounds_x(gz0,gz6,0);

		LReflect(H,gz6,gz0);
		YplusisCtimesX(gx2+JX,H,2*P[1],M-JX);
		YplusisCtimesX(gx3+JX,H,2*P[1],M-JX);
		YplusisCtimesX(gx4+JX,H,2*P[1],M-JX);
		YplusisCtimesX(gx5+JX,H,2*P[0],M-JX);
		YplusisCtimesX(gx6+JX,H,2*P[0],M-JX);

		UReflect(H,gz0,gz6);
		YplusisCtimesX(gx0,H+JX,2*P[0],M-JX);
		YplusisCtimesX(gx1,H+JX,2*P[0],M-JX);
		YplusisCtimesX(gx2,H+JX,2*P[1],M-JX);
		YplusisCtimesX(gx3,H+JX,2*P[1],M-JX);
		YplusisCtimesX(gx4,H+JX,2*P[1],M-JX);

		set_bounds_y(gz2,gz4,0);

		YplusisCtimesX(gx0+JY,gz4,2*P[1],M-JY);
		YplusisCtimesX(gx1+JY,gz4,2*P[1],M-JY);
		YplusisCtimesX(gx3+JY,gz4,P[0],M-JY);
		YplusisCtimesX(gx4+JY,gz4,2*P[0],M-JY);
		YplusisCtimesX(gx5+JY,gz4,2*P[1],M-JY);
		YplusisCtimesX(gx6+JY,gz4,2*P[1],M-JY);

		YplusisCtimesX(gx0,gz2+JY,2*P[1],M-JY);
		YplusisCtimesX(gx1,gz2+JY,2*P[1],M-JY);
		YplusisCtimesX(gx2,gz2+JY,2*P[0],M-JY);
		YplusisCtimesX(gx3,gz2+JY,P[0],M-JY);
		YplusisCtimesX(gx5,gz2+JY,2*P[1],M-JY);
		YplusisCtimesX(gx6,gz2+JY,2*P[1],M-JY);


		YplusisCtimesX(gx0,gz3,2*P[1],M);
		YplusisCtimesX(gx1,gz3,2*P[1],M);
		YplusisCtimesX(gx2,gz3,P[0],M);
		YplusisCtimesX(gx3,gz3,P[0],M);
		YplusisCtimesX(gx4,gz3,P[0],M);
		YplusisCtimesX(gx5,gz3,2*P[1],M);
		YplusisCtimesX(gx6,gz3,2*P[1],M);

		set_bounds_y(gz1,gz5,-1);
		//set_bounds_y(gz1,-1);set_bounds_y(gz5,-1);

		LReflect(H,gz5,gz1);
		YplusisCtimesX(gx2+JX,H+JY,P[1],M-JX-JY);
		YplusisCtimesX(gx3+JX,H+JY,P[1],M-JX-JY);
		YplusisCtimesX(gx4+JX,H+JY,P[1],M-JX-JY);
		YplusisCtimesX(gx5+JX,H+JY,P[0],M-JX-JY);
		YplusisCtimesX(gx6+JX,H+JY,P[0],M-JX-JY);

		remove_bounds(gz1);remove_bounds(gz5);
		set_bounds_x(gz1,gz5,-1);
		//set_bounds_x(gz1,-1);set_bounds_x(gz5,-1);

		UReflect(H,gz1,gz5);
		YplusisCtimesX(gx0+JY,H+JX,P[0],M-JY-JX);
		YplusisCtimesX(gx1+JY,H+JX,P[0],M-JY-JX);
		YplusisCtimesX(gx2+JY,H+JX,P[1],M-JY-JX);
		YplusisCtimesX(gx3+JY,H+JX,P[1],M-JY-JX);
		YplusisCtimesX(gx4+JY,H+JX,P[1],M-JY-JX);

		for (int k=0; k<7; k++) Times(gs+k*M,gs+k*M,g,M);

	} else {
		Real *gs=G+M*5*s_to;
		Real *gs_1=G+M*5*s_from;
		Real *gz0=gs_1;
		Real *gz1=gs_1+M;
		Real *gz2=gs_1+2*M;
		Real *gz3=gs_1+3*M;
		Real *gz4=gs_1+4*M;
		set_bounds_x(gz0,gz4,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0); set_bounds_x(gz3,0);
		set_bounds_y(gz1,gz3,0); set_bounds_y(gz0,0); set_bounds_y(gz2,0); set_bounds_y(gz4,0);
		Real *gx0=gs;
		Real *gx1=gs+M;
		Real *gx2=gs+2*M;
		Real *gx3=gs+3*M;
		Real *gx4=gs+4*M;
		Real *g=G1;

		Zero(gs,5*M);

		LReflect(H,gz4,gz0);
		YplusisCtimesX(gx1+JX,H,P[1],M-JX);
		YplusisCtimesX(gx2+JX,H,P[1],M-JX);
		YplusisCtimesX(gx3+JX,H,P[1],M-JX);
		YplusisCtimesX(gx4+JX,H,P[0],M-JX);

		YplusisCtimesX(gx0+JY,gz3,P[1],M-JY);
		YplusisCtimesX(gx2+JY,gz3,P[1],M-JY);
		YplusisCtimesX(gx3+JY,gz3,P[0],M-JY);
		YplusisCtimesX(gx4+JY,gz3,P[1],M-JY);

		YplusisCtimesX(gx0,gz2,2*P[1],M);
		YplusisCtimesX(gx1,gz2,2*P[1],M);
		YplusisCtimesX(gx2,gz2,P[0],M);
		YplusisCtimesX(gx3,gz2,2*P[1],M);
		YplusisCtimesX(gx4,gz2,2*P[1],M);

		YplusisCtimesX(gx0,gz1+JY,P[1],M-JY);
		YplusisCtimesX(gx1,gz1+JY,P[0],M-JY);
		YplusisCtimesX(gx2,gz1+JY,P[1],M-JY);
		YplusisCtimesX(gx4,gz1+JY,P[1],M-JY);

		UReflect(H,gz0,gz4);
		YplusisCtimesX(gx0,H+JX,P[0],M-JX);
		YplusisCtimesX(gx1,H+JX,P[1],M-JX);
		YplusisCtimesX(gx2,H+JX,P[1],M-JX);
		YplusisCtimesX(gx3,H+JX,P[1],M-JX);

		for (int k=0; k<5; k++) Times(gs+k*M,gs+k*M,g,M);
	}
}

void LGrad2::propagate(Real *G, Real *G1, int s_from, int s_to,int M) {
if (debug) cout <<" propagate in LGrad2 " << endl;
	Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	Zero(gs,M); set_bounds(gs_1);
	switch (fjc) {
		case 1:
			if (lattice_type==simple_cubic) {
				YplusisCtimesX(gs,gs_1,4.0/6.0,M);
				AddTimes(gs+JX,gs_1,lambda_1+JX,M-JX);
				AddTimes(gs,gs_1+JX,lambda1,    M-JX);
				YplusisCtimesX(gs+1,gs_1,1.0/6.0,M-1);
				YplusisCtimesX(gs,gs_1+1,1.0/6.0,M-1);
				Norm(gs,4.0,M);
				AddTimes(gs+JX+1,gs_1,lambda_1+JX+1,M-JX-1);
				AddTimes(gs+JX,gs_1+1,lambda_1+JX,  M-JX-1);
				AddTimes(gs+1,gs_1+JX,lambda1+1,M-JX-1);
				AddTimes(gs,gs_1+JX+1,lambda1,  M-JX-1);
				Norm(gs,1.0/6.0,M);
				Times(gs,gs,G1,M);
			} else { //9 point stencil; hexagonal

				YplusisCtimesX(gs,gs_1,0.5,M);
				AddTimes(gs+JX,gs_1,lambda_1+JX,M-JX);
				AddTimes(gs,gs_1+JX,lambda1,    M-JX);
				YplusisCtimesX(gs+1,gs_1,0.25,M-1);
				YplusisCtimesX(gs,gs_1+1,0.25,M-1);
				Norm(gs,2.0,M);
				AddTimes(gs+JX+1,gs_1,lambda_1+JX+1,M-JX-1);
				AddTimes(gs+JX,gs_1+1,lambda_1+JX,	M-JX-1);
				AddTimes(gs+1,gs_1+JX,lambda1+1,M-JX-1);
				AddTimes(gs,gs_1+JX+1,lambda1,  M-JX-1);
				Norm(gs,0.25,M);
				Times(gs,gs,G1,M);
			}
			break;
		case 2:  //25 point stencil
			AddTimes(gs,     gs_1,         LAMBDA+2*M,     M);

			AddTimes(gs+1,   gs_1,         LAMBDA+2*M+1,   M-1);
			AddTimes(gs,     gs_1+1,       LAMBDA+2*M,     M-1);
			AddTimes(gs+JX,  gs_1,         LAMBDA+M+JX,    M-JX);
			AddTimes(gs,     gs_1+JX,      LAMBDA+3*M,     M-JX);

			AddTimes(gs+JX+1,gs_1,         LAMBDA+M+JX+1,  M-JX-1);
			AddTimes(gs+JX,  gs_1+1,       LAMBDA+M+JX,    M-JX-1);
			AddTimes(gs,     gs_1+JX+1,    LAMBDA+3*M,     M-JX-1);
			AddTimes(gs+1,   gs_1+JX,      LAMBDA+3*M+1,   M-JX-1);

			AddTimes(gs+2*JX,  gs_1,       LAMBDA+  2*JX,  M-2*JX);
			AddTimes(gs,       gs_1+2*JX,  LAMBDA+4*M,     M-2*JX);

			AddTimes(gs+2*JX,  gs_1+1,     LAMBDA+  2*JX,  M-2*JX-1);
			AddTimes(gs+2*JX+1,gs_1,       LAMBDA+  2*JX+1,M-2*JX-1);
			AddTimes(gs+1,     gs_1+2*JX,  LAMBDA+4*M+1,   M-2*JX-1);
			AddTimes(gs,       gs_1+2*JX+1,LAMBDA+4*M,     M-2*JX-1);

			Norm(gs,2.0,M);
			AddTimes(gs+2,     gs_1,       LAMBDA+2*M+2,   M-2);
			AddTimes(gs,       gs_1+2,     LAMBDA+2*M,     M-2);

			AddTimes(gs+JX+2,  gs_1,       LAMBDA+1*M+JX+2,M-JX-2);
			AddTimes(gs+JX,    gs_1+2,     LAMBDA+1*M+JX,  M-JX-2);
			AddTimes(gs,       gs_1+JX+2,  LAMBDA+3*M,     M-JX-2);
			AddTimes(gs+2,     gs_1+JX,    LAMBDA+3*M+2,   M-JX-2);

			AddTimes(gs+2*JX+2,gs_1,       LAMBDA+   2*JX+2,M-2*JX-2);
			AddTimes(gs,       gs_1+2*JX+2,LAMBDA+4*M,     M-2*JX-2);
			AddTimes(gs+2*JX,  gs_1+2,     LAMBDA+   2*JX, M-2*JX-2);
			AddTimes(gs+2,     gs_1+2*JX,  LAMBDA+4*M+2,   M-2*JX-2);
			Norm(gs,1.0/8.0,M);
			Times(gs,gs,G1,M);
			break;
		case 3:

			AddTimes(gs,        gs_1,        LAMBDA+3*M,         M);

			AddTimes(gs+1,      gs_1,        LAMBDA+3*M+1,       M-1);
			AddTimes(gs,        gs_1+1,      LAMBDA+3*M,         M-1);
			AddTimes(gs+JX,     gs_1,        LAMBDA+2*M+JX,      M-JX);
			AddTimes(gs,        gs_1+JX,     LAMBDA+4*M,         M-JX);

			AddTimes(gs+JX+1,   gs_1,        LAMBDA+2*M+JX+1,    M-JX-1);
			AddTimes(gs,        gs_1+JX+1,   LAMBDA+4*M,         M-JX-1);
			AddTimes(gs+1,      gs_1+JX,     LAMBDA+4*M+1,       M-JX-1);
			AddTimes(gs+JX,     gs_1+1,      LAMBDA+2*M+JX,      M-JX-1);

			AddTimes(gs+2*JX,   gs_1,        LAMBDA+M+2*JX,      M-2*JX);
			AddTimes(gs,        gs_1+2*JX,   LAMBDA+5*M,         M-2*JX);
			AddTimes(gs+2,      gs_1,        LAMBDA+3*M+2,       M-2);
			AddTimes(gs,        gs_1+2,      LAMBDA+3*M,         M-2);

			AddTimes(gs+2*JX+1, gs_1,        LAMBDA+M+2*JX+1,    M-2*JX-1);
			AddTimes(gs,        gs_1+2*JX+1, LAMBDA+5*M,         M-2*JX-1);
			AddTimes(gs+1,      gs_1+2*JX,   LAMBDA+5*M+1,       M-2*JX-1);
			AddTimes(gs+2*JX,   gs_1+1,      LAMBDA+M+2*JX,      M-2*JX-1);

			AddTimes(gs+JX+2,   gs_1,        LAMBDA+2*M+JX+2,    M-JX-2);
			AddTimes(gs,        gs_1+JX+2,   LAMBDA+4*M,         M-JX-2);
			AddTimes(gs+2,      gs_1+JX,     LAMBDA+4*M+2,       M-JX-2);
			AddTimes(gs+JX,     gs_1+2,      LAMBDA+2*M+JX,      M-JX-2);

			AddTimes(gs+2*JX+2, gs_1,        LAMBDA+1*M+2*JX+2,  M-2*JX-2);
			AddTimes(gs,        gs_1+2*JX+2, LAMBDA+5*M,         M-2*JX-2);
			AddTimes(gs+2,      gs_1+2*JX,   LAMBDA+5*M+2,       M-2*JX-2);
			AddTimes(gs+2*JX,   gs_1+2,      LAMBDA+1*M+2*JX,    M-2*JX-2);

			AddTimes(gs+3*JX,   gs_1,        LAMBDA+   3*JX,   M-3*JX);
			AddTimes(gs,        gs_1+3*JX,   LAMBDA+6*M,       M-3*JX);

			AddTimes(gs+3*JX,   gs_1+1,      LAMBDA+   3*JX,   M-3*JX-1);
			AddTimes(gs+3*JX+1, gs_1,        LAMBDA+   3*JX+1, M-3*JX-1);
			AddTimes(gs+1,      gs_1+3*JX,   LAMBDA+6*M+1,     M-3*JX-1);
			AddTimes(gs,        gs_1+3*JX+1, LAMBDA+6*M,       M-3*JX-1);

			AddTimes(gs+3*JX,   gs_1+2,      LAMBDA+   3*JX,   M-3*JX-2);
			AddTimes(gs+3*JX+2, gs_1,        LAMBDA+   3*JX+2, M-3*JX-2);
			AddTimes(gs+2,      gs_1+3*JX,   LAMBDA+6*M+2,     M-3*JX-2);
			AddTimes(gs,        gs_1+3*JX+2, LAMBDA+6*M,       M-3*JX-2);

			Norm(gs,2.0,M);

			AddTimes(gs+3,      gs_1,        LAMBDA+3*M+3,     M-3);
			AddTimes(gs,        gs_1+3,      LAMBDA+3*M,       M-3);

			AddTimes(gs+JX+3,   gs_1,        LAMBDA+2*M+JX+3,  M-JX-3);
			AddTimes(gs+JX,     gs_1+3,      LAMBDA+2*M+JX,    M-JX-3);
			AddTimes(gs,        gs_1+JX+3,   LAMBDA+4*M,       M-JX-3);
			AddTimes(gs+3,      gs_1+JX,     LAMBDA+4*M+3,     M-JX-3);

			AddTimes(gs+2*JX+3, gs_1,        LAMBDA+1*M+2*JX+3,M-2*JX-3);
			AddTimes(gs+2*JX,   gs_1+3,      LAMBDA+1*M+2*JX,  M-2*JX-3);
			AddTimes(gs,        gs_1+2*JX+3, LAMBDA+5*M,       M-2*JX-3);
			AddTimes(gs+3,      gs_1+2*JX,   LAMBDA+5*M+3,     M-2*JX-3);

			AddTimes(gs+3*JX+3, gs_1,        LAMBDA+   3*JX+3, M-3*JX-3);
			AddTimes(gs,        gs_1+3*JX+3, LAMBDA+6*M,       M-3*JX-3);
			AddTimes(gs+3*JX,   gs_1+3,      LAMBDA+   3*JX,   M-3*JX-3);
			AddTimes(gs+3,      gs_1+3*JX,   LAMBDA+6*M+3,     M-3*JX-3);

			Norm(gs,1.0/10.0,M);
			Times(gs,gs,G1,M);
			break;
		default:
			break;
	}
}


bool LGrad2::ReadRange(int* r, int* H_p, int &n_pos, bool &block, string range, int var_pos, string seg_name, string range_type) {
if (debug) cout <<"ReadRange in LGrad2 " << endl;
	bool success=true;
	vector<string>set;
	vector<string>coor;
	vector<string>xyz;
	In[0]->split(range,';',set);
	coor.clear();
	block=true; In[0]->split(set[0],',',coor);
	if (coor.size()!=2) {cout << "In mon " + 	seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes do not come in set of two: 'x,y'" << endl; success=false;}
	else {
		r[0]=In[0]->Get_int(coor[0],0);
		r[1]=In[0]->Get_int(coor[1],0);
	}
	coor.clear(); In[0]->split(set[1],',',coor);

	if (coor.size()!=2) {cout << "In mon " + seg_name+ ", for 'pos 2', in '" + range_type + "', the coordinates do not come in set of two: 'x,y'" << endl; success=false;}
	else {
		r[3]=In[0]->Get_int(coor[0],0);
		r[4]=In[0]->Get_int(coor[1],0);
	}
	if (r[0] > r[3]) {cout << "In mon " + seg_name+ ", for 'pos 1', the x-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
	if (r[1] > r[4]) {cout << "In mon " + seg_name+ ", for 'pos 1', the y-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}

	return success;
}

bool LGrad2::ReadRangeFile(string filename,int* H_p, int &n_pos, string seg_name, string range_type) {
if (debug) cout <<"ReadRangeFile in LGrad2 " << endl;
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
	int px,py,p_i,x,y;
	int i=0;
	if (!In[0]->ReadFile(sub[0].append(".").append(filename),content)) {
		success=false;
		return success;
	}

	In[0]->split(content,'#',lines);
	length = lines.size();


	if (length == MX*MY) { //expect to read 'mask file';
		i=0;
		if (n_pos==0) {
			while (i<length){
				if (In[0]->Get_int(lines[i],0)==1) n_pos++;
				i++;
			};
			if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any unities." << endl;}
		} else {
			i=0; p_i=0;
			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++)  {
				if (In[0]->Get_int(lines[i],0)==1) {H_p[p_i]=P(x,y); p_i++;}
				i++;
			}
		}
	} else { //expect to read x,y
		px=0,py=0; i=0;
		if (n_pos==0) n_pos=length;
		else {
			while (i<length) {
				xyz.clear();
				In[0]->split(lines[i],',',xyz);
				length_xyz=xyz.size();
				if (length_xyz!=2) {
					cout << "In mon " + seg_name + " " +range_type+"_filename  the expected 'pair of coordinates' 'x,y' was not found. " << endl;  success = false;
				} else {
					px=In[0]->Get_int(xyz[0],0);
					if (px < 1 || px > MX) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the x-coordinate in "+range_type+"_filename out of bounds: 1.." << MX << endl; success =false;}
					py=In[0]->Get_int(xyz[1],0);
					if (py < 1 || py > MY) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the y-coordinate in "+range_type+"_filename out of bounds: 1.." << MY << endl; success =false;}
				}
				cout <<"reading px " << px << " and py " << py << endl;
				H_p[i]=P(px,py);
				i++;
			}
		}
	}
	return success;
}

bool LGrad2::FillMask(int* Mask, vector<int>px, vector<int>py, vector<int>pz, string filename) {
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
		if (MX*MY!=length) {success=false; cout <<"inputfile for filling delta_range has not expected length in x,y-directions" << endl;
		} else {
			for (int x=1; x<MX+1; x++)
			for (int y=1; y<MY+1; y++) Mask[x*JX + y]=In[0]->Get_int(lines[x*JX+y],-1);
		}
	} else  {
		for (int i=0; i<length_px; i++) {
			p=px[i]; if (p<1 || p>MX) {success=false; cout<<" x-value in delta_range out of bounds; " << endl; }
			p=py[i]; if (p<1 || p>MY) {success=false; cout<<" y-value in delta_range out of bounds; " << endl; }
			if (success) Mask[P(px[i],py[i])]=1; //Mask[px[i]*JX + fjc-1+ py[i]]=1;
		}
	}
	for (int i=0; i<M; i++) if (!(Mask[i]==0 || Mask[i]==1)) {success =false; cout <<"Delta_range does not contain '0' or '1' values. Check delta_inputfile values"<<endl; }
	return success;
}

bool LGrad2::CreateMASK(int* H_MASK, int* r, int* H_P, int n_pos, bool block) {
if (debug) cout <<"CreateMask for LGrad2 " + name << endl;
	bool success=true;
	H_Zero(H_MASK,M);
	if (block) {
		for (int x=r[0]; x<r[3]+1; x++)
		for (int y=r[1]; y<r[4]+1; y++)
			H_MASK[P(x,y)]=1;
	} else {
		for (int i = 0; i<n_pos; i++) H_MASK[H_P[i]]=1;
	}
	return success;
}


Real LGrad2::ComputeTheta(Real* phi) {
	Real result=0; remove_bounds(phi);
	if (geometry !="planar") Dot(result,phi,L,M);
	else {if (fjc==1) Sum(result,phi,M); else  Dot(result,phi,L,M);}
	return result;
}

void LGrad2::UpdateEE(Real* EE, Real* psi, Real* E) {
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_M_bounds(psi);
	Zero(EE,M);
	Real Exmin,Explus,Eymin,Eyplus;
	int x,y,z;
	int r;

	pf = pf/2;
	r=offset_first_layer*fjc;
	for (x=fjc; x<MX+fjc; x++) {
		r++;
		for (y=fjc; y<MY+fjc; y++) {
			z=x*JX+y;
			Exmin=psi[z]-psi[z-JX];
			Exmin*=(r-0.5)*Exmin*2*PIE;
			Explus=psi[z]-psi[z+JX];
			Explus*=(r+0.5)*Explus*2*PIE;
			Eymin=(psi[z]-psi[z-1])*fjc;
			Eymin*=Eymin;
			Eyplus=(psi[z]-psi[z+1])*fjc;
			Eyplus*=Eyplus;
			EE[z]=pf*((Exmin+Explus)/L[z]+Eymin+Eyplus);
		}
	}
}


void LGrad2::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, int* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	int x,y,i;

	Real r;
	Real epsXplus, epsXmin, epsYplus,epsYmin;
	//set_M_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);

	if (!fixedPsi0) {
		C=C/2/fjc/fjc;
		r=offset_first_layer*fjc;
		for (x=fjc; x<MX+fjc; x++) {
			r++;
			for (y=fjc; y<MY+fjc; y++) {
				i=x*JX+y;
				epsXmin=2*PIE*(r-1)*(eps[i]+eps[i-JX])/L[i]*fjc*fjc;
				epsXplus=2*PIE*r*(eps[i]+eps[i+JX])/L[i]*fjc*fjc;
				epsYmin=eps[i]+eps[i-1];
				epsYplus=eps[i]+eps[i+1];
				X[i]= (C*q[i]+epsXmin*psi[i-JX]+epsXplus*psi[i+JX]+epsYmin*psi[i-1]+epsYplus*psi[i+1])/
				(epsXmin+epsXplus+epsYmin+epsYplus);
			}
		}
		//Cp(psi,X,M);
		YisAminB(g,g,X,M);
  	 } else { //fixedPsi0 is true
		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++) {
				if (Mask[x*JX+y] == 0)
				X[x*JX+y]=0.25*(psi[(x-1)*JX+y]+psi[(x+1)*JX+y]
			        +psi[x*JX+y-1]  +psi[x*JX+y+1])
				 +0.5*q[x*JX+y]*C/eps[x*JX+y];
			}
		}


		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++){
				if (Mask[x*JX+y] == 0)
				X[x*JX+y]+=(psi[(x+1)*JX+y]-psi[(x-1)*JX+y])/(2.0*(offset_first_layer+x-fjc+0.5))*fjc;
			}
		}
		if (grad_epsilon) {
			for (x=fjc; x<MX+fjc; x++) {
				for (y=fjc; y<MY+fjc; y++) {
					if (Mask[x*JX+y] == 0) {
						X[x*JX+y]+=0.25*(eps[(x+1)*JX+y]-eps[(x-1)*JX+y])*(psi[(x+1)*JX+y]-psi[(x-1)*JX+y]+
                                              eps[x*JX+y+1]  -eps[x*JX+y-1])  *(psi[x*JX+y+1]  -psi[x*JX+y-1])/
					           eps[x*JX+y]*fjc*fjc;
					}
				}
			}
		}
		for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)
		if (Mask[x*JX+y] == 0) {
			psi[x*JX+y]=X[x*JX+y];
			g[x*JX+y]-=psi[x*JX+y];
		}
   	}
}


void LGrad2::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, int* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
	int x,y;

	Real C = -e*e/(eps0*k_BT*bond_length);
	for (x=fjc; x<MX+fjc; x++) {
		for (y=fjc; y<MY+fjc; y++){ //for all geometries
			if (Mask[x*JX+y] == 1)
			q[x*JX+y]=-0.5*(psi[(x-1)*JX+y]+psi[(x+1)*JX+y]
				        +psi[x*JX+y-1]  +psi[x*JX+y+1]
					 -4*psi[x*JX+y])*fjc*fjc*eps[x*JX+y]/C;
		}
	}

	for (x=fjc; x<MX+fjc; x++) {
		for (y=fjc; y<MY+fjc; y++){
			if (Mask[x*JX+y] == 1)
			q[x*JX+y]-=(psi[(x+1)*JX+y]-psi[(x-1)*JX+y])/(2.0*(offset_first_layer+x-fjc+0.5))*fjc*eps[x]/C;
		}
	}
	if (grad_epsilon) {
		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++) {
				if (Mask[x*JX+y] == 1)
				q[x*JX+y]-=0.25*(eps[(x+1)*JX+y]-eps[(x-1)*JX+y])*(psi[(x+1)*JX+y]-psi[(x-1)*JX+y]+
                                  		   eps[x*JX+y+1]  -eps[x*JX+y-1])  *(psi[x*JX+y+1]  -psi[x*JX+y-1])*fjc*fjc/C;
				}
		}
	}
	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)
	if (Mask[x*JX+y] == 1) {
		g[x*JX+y]=-q[x*JX+y];
	}

}

void LGrad2::remove_bounds(Real *X){
if (debug) cout <<"remove_bounds in LGrad2 " << endl;
	int x,y;
	int k=0;
	if (fjc==1) {
		for (x=1; x<MX+1; x++) {
			X[x*JX+0] = 0;
			X[x*JX+MY+1]=0;
		}
		for (y=1; y<MY+1; y++) {
			X[0+y] = 0;
			X[(MX+1)*JX+y]=0;
		}
		//corners
		for (x=0; x<1; x++) {
			X[x*JX+0] = 0;
			X[x*JX+MY+1]=0;
		}
		for (x=MX+1; x<MX+2; x++) {
			X[x*JX+0] = 0;
			X[x*JX+MY+1]=0;
		}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+k]=0;
				X[x*JX+MY+fjc+k]=0;
			}
		}
		for (y=0; y<MY+2*fjc; y++) { //this will also remove the corners...
			for (k=0; k<fjc; k++) {
				X[k*JX+y]=0;
				X[(MX+fjc+k)*JX+y]=0;
			}
		}
	}
}

void LGrad2::set_bounds_x(Real* X,Real*Y, int shifty){
if (debug) cout <<"set_bounds_x XY in LGrad2 " << endl;
	int y;
	int k=0;
	if (BX1>BXM)  {
		set_bounds_x(X,0); set_bounds_x(Y,0);
	} else  {
		if (fjc==1) {
			for (y=1; y<MY+1; y++) {
				X[0        +y]= Y[BX1*JX+(y+shifty)];
				X[(MX+1)*JX+y]= Y[BXM*JX+(y-shifty)];
				Y[0        +y]= X[BX1*JX+(y+shifty)];
				Y[(MX+1)*JX+y]= X[BXM*JX+(y-shifty)];

			}
		} else {
			cout <<"set_bounds_x error" << endl;
			for (y=0; y<MY+2*fjc; y++) { //this will also set the corners...fingers crossed ; this might go wrong when reflecting and periodic b.c. are mixed...
				for (k=0; k<fjc; k++) {
					X[k*JX+y]=Y[B_X1[k]*JX+(y+shifty)];
					X[(MX+fjc+k)*JX+y]=Y[B_XM[k]*JX+(y-shifty)];
					Y[k*JX+y]=X[B_X1[k]*JX+(y+shifty)];
					Y[(MX+fjc+k)*JX+y]=X[B_XM[k]*JX+(y-shifty)];
				}
			}
		}
	}
}

void LGrad2::set_bounds_y(Real* X,Real*Y, int shiftx){
if (debug) cout <<"set_bounds_y XY in LGrad2 " << endl;
	int x;
	int k=0;
	if (BY1>BYM) {
		set_bounds_y(X,0); set_bounds_y(Y,0);

	} else {
		if (fjc==1) {
			for (x=1; x<MX+1; x++) {
				X[x*JX+0   ] =Y[(x+shiftx)*JX+BY1];
				X[x*JX+MY+1] =Y[(x-shiftx)*JX+BYM];
				Y[x*JX+0   ] =X[(x+shiftx)*JX+BY1];
				Y[x*JX+MY+1] =X[(x-shiftx)*JX+BYM];

			}
		} else {
			for (x=fjc; x<MX+fjc; x++) {
				for (k=0; k<fjc; k++) {
					X[x*JX+k]=Y[(x+shiftx)*JX+B_Y1[k]];
					X[x*JX+MY+fjc+k]=Y[(x-shiftx)*JX+B_YM[k]];
					Y[x*JX+k]=X[(x+shiftx)*JX+B_Y1[k]];
					Y[x*JX+MY+fjc+k]=X[(x-shiftx)*JX+B_YM[k]];
				}
			}
		}
	}
}

void LGrad2::set_bounds_x(Real* X,int shifty){
if (debug) cout <<"set_bounds_x X in LGrad2 " << endl;
	int y;
	int k=0;
	//if (BX1>BXM) shifty=0; //periodic

	if (fjc==1) {
		for (y=1; y<MY+1; y++) { //not yet 0 and MY+1....
			X[0        +y] = X[BX1*JX+(y+shifty)];
			X[(MX+1)*JX+y] = X[BXM*JX+(y-shifty)];
		}
	} else {
		for (y=0; y<MY+2*fjc; y++) { //this will also set the corners...fingers crossed ; this might go wrong when reflecting and periodic b.c. are mixed...
			for (k=0; k<fjc; k++) {
				X[k*JX+y]=X[B_X1[k]*JX+(y+shifty)];
				X[(MX+fjc+k)*JX+y]=X[B_XM[k]*JX+(y-shifty)];
			}
		}
	}
}

void LGrad2::set_bounds_y(Real* X,int shiftx){
if (debug) cout <<"set_bounds_y X in LGrad2 " << endl;
	int x;
	int k=0;
	//if (BY1>BYM) shiftx=0;

	if (fjc==1) {
		for (x=1; x<MX+1; x++) {
			X[x*JX+0   ]= X[(x+shiftx)*JX+BY1];
			X[x*JX+MY+1]= X[(x-shiftx)*JX+BYM];
		}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+k]=X[(x+shiftx)*JX+B_Y1[k]];
				X[x*JX+MY+fjc+k]=X[(x-shiftx)*JX+B_YM[k]];
			}
		}
	}
}

void LGrad2::set_bounds(Real* X){
if (debug) cout <<"set_bounds in LGrad2 " << endl;
	int x,y;
	int k=0;
	if (fjc==1) {
		for (x=1; x<MX+1; x++) {
			X[x*JX+0] = X[x*JX+BY1];
			X[x*JX+MY+1]=X[x*JX+BYM];
		}
		for (y=1; y<MY+1; y++) {
			X[0+y] = X[BX1*JX+y];
			X[(MX+1)*JX+y]=X[BXM*JX+y];
		}
		//corners....
		//for (x=0; x<1; x++) {
			X[0] = X[BY1];
			X[MY+1]=X[BYM];
		//}
		//for (x=MX+1; x<MX+2; x++) {
			X[(MX+1)*JX+0] = X[(MX+1)*JX+BY1];
			X[(MX+1)*JX+MY+1]=X[(MX+1)*JX+BYM];
		//}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+k]=X[x*JX+B_Y1[k]];
				X[x*JX+MY+fjc+k]=X[x*JX+B_YM[k]];
			}
		}
		for (y=fjc; y<MY+fjc; y++) {
			for (k=0; k<fjc; k++) {
				X[k*JX+y]=X[B_X1[k]*JX+y];
				X[(MX+fjc+k)*JX+y]=X[B_XM[k]*JX+y];
			}
		}
		for (int k=0; k<fjc; k++) for (int m=0; m<fjc; m++)  {
			X[k*JX         +m      ]=X[k*JX+B_Y1[m]];
			X[(MX+fjc+k)*JX+m      ]=X[(MX+fjc+k)*JX+B_Y1[m]];
			X[k*JX         +(MY+fjc+m)]=X[k*JX+B_YM[m]];
			X[(MX+fjc+k)*JX+(MY+fjc+m)]=X[(MX+fjc+k)*JX+B_YM[m]];
			//X[k*JX         +m      ]=X[B_X1[k]*JX+B_Y1[m]];
			//X[(MX+fjc+k)*JX+m      ]=X[B_XM[k]*JX+B_Y1[m]];
			//X[k*JX         +(MY+fjc+m)]=X[B_X1[k]*JX+B_YM[m]];
			//X[(MX+fjc+k)*JX+(MY+fjc+m)]=X[B_XM[k]*JX+B_YM[m]];
		}
	}
}

void LGrad2::set_M_bounds(Real* X){
if (debug) cout <<"set_bounds in LGrad2 " << endl;
	int x,y;
	int k=0;
	if (fjc==1) {
		for (x=1; x<MX+1; x++) {
			X[x*JX+0] = X[x*JX+1];
			X[x*JX+MY+1]=X[x*JX+MY];
		}
		for (y=1; y<MY+1; y++) {
			X[0+y] = X[1*JX+y];
			X[(MX+1)*JX+y]=X[MX*JX+y];
		}
		//corners
		for (x=0; x<1; x++) {
			X[x*JX+0] = X[x*JX+1];
			X[x*JX+MY+1]=X[x*JX+MY];
		}
		for (x=MX+1; x<MX+2; x++) {
			X[x*JX+0] = X[x*JX+1];
			X[x*JX+MY+1]=X[x*JX+MY];
		}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+k]=X[x*JX+2*fjc-1-k];
				X[x*JX+MY+fjc+k]=X[x*JX+MY+fjc-k-1];
			}
		}
		for (y=0; y<MY+2*fjc; y++) { //this will also set the corners...fingers crossed ; this might go wrong when reflecting and periodic b.c. are mixed...
			for (k=0; k<fjc; k++) {
				X[k*JX+y]=X[(2*fjc-1-k)*JX+y];
				X[(MX+fjc+k)*JX+y]=X[(MX+fjc-k-1)*JX+y];
			}
		}
	}
}


void LGrad2::remove_bounds(int *X){
if (debug) cout <<"remove_bounds in LGrad2 " << endl;
	int x,y;
	int k;
	if (fjc==1) {
		for (x=0; x<MX+2; x++) {
			X[P(x,0)] = 0;  //needs testing if this is okay to put to zero in case of surface...
			X[P(x,MY+1)]=0;
		}
		for (y=0; y<MY+2; y++) {
			X[P(0,y)] = 0;
			X[P(MX+1,y)]=0;
		}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+k]=0;
				X[x*JX+MY+fjc+k]=0;
			}
		}
		for (y=0; y<MY+2*fjc; y++) { //this will also set the corners...fingers crossed...
			for (k=0; k<fjc; k++) {
				X[k*JX+y]=0;
				X[(MX+fjc+k)*JX+y]=0;
			}
		}
	}
}

void LGrad2::set_bounds(int* X){
if (debug) cout <<"set_bounds in LGrad2 " << endl;
	int x,y;
	int k=0;
	if (fjc==1) {
		for (x=1; x<MX+1; x++) {
			X[x*JX+0] = X[x*JX+BY1];
			X[x*JX+MY+1]=X[x*JX+BYM];
		}
		for (y=1; y<MY+1; y++) {
			X[0+y] = X[BX1*JX+y];
			X[(MX+1)*JX+y]=X[BXM*JX+y];
		}
		//corners
		//for (x=0; x<1; x++) {
			X[0] = X[BY1];
			X[MY+1]=X[BYM];
		//}
		//for (x=MX+1; x<MX+2; x++) {
			X[(MX+1)*JX+0] = X[(MX+1)*JX+BY1];
			X[(MX+1)*JX+MY+1]=X[(MX+1)*JX+BYM];
		//}
	} else {
		for (x=fjc; x<MX+fjc; x++) {
			for (k=0; k<fjc; k++) {
				X[x*JX+k]=X[x*JX+B_Y1[k]];
				X[x*JX+MY+fjc+k]=X[x*JX+B_YM[k]];
			}
		}
		for (y=fjc; y<MY+fjc; y++) { //this will also set the corners...fingers crossed ; this might go wrong when reflecting and periodic b.c. are mixed...
			for (k=0; k<fjc; k++) {
				X[k*JX+y]=X[B_X1[k]*JX+y];
				X[(MX+fjc+k)*JX+y]=X[B_XM[k]*JX+y];
			}
		}
		for (int k=0; k<fjc; k++) for (int m=0; m<fjc; m++)  {
			X[k*JX         +m      ]   =X[B_X1[k]*JX+B_Y1[m]];
			X[(MX+fjc+k)*JX+m      ]   =X[B_XM[k]*JX+B_Y1[m]];
			X[k*JX         +(MY+fjc+m)]=X[B_X1[k]*JX+B_YM[m]];
			X[(MX+fjc+k)*JX+(MY+fjc+m)]=X[B_XM[k]*JX+B_YM[m]];
		}
	}
}

Real LGrad2::ComputeGN(Real* G,int Markov, int M){
	Real GN=0;
	if (Markov==2) {
		if (lattice_type == hexagonal) {
			for (int k=0; k<7; k++) {
				if (k==1 || k==5)
					GN += WeightedSum(G+k*M);
				else
					GN += 2*WeightedSum(G+k*M);
			}
			GN /=12.0;

		} else {
			for (int k=0; k<5; k++) {
				if (k== 2)
					GN += 2*WeightedSum(G+k*M);
				else
					GN += WeightedSum(G+k*M);
			}
			GN /=6.0;
		}
	} else GN = WeightedSum(G);
	return GN;
}
void LGrad2::AddPhiS(Real* phi,Real* Gf,Real* Gb,int Markov, int M){
	Real one=1.0;
	if (Markov==2) {
		if (lattice_type == hexagonal) {
			for (int k=0; k<7; k++) {
				if (k==1 || k==5)
					YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,1.0/12.0*one,M);
				else
					YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,1.0/6.0*one,M);
			}

		} else {
			for (int k=0; k<5; k++) {
				if (k==2)
					YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,1.0/3.0*one,M);
				else
					YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,1.0/6.0*one,M);
			}
		}
	} else AddTimes(phi,Gf,Gb,M);
}
void LGrad2::AddPhiS(Real* phi,Real* Gf,Real* Gb,Real degeneracy, int Markov, int M){
	if (Markov==2) {
		if (lattice_type == hexagonal) {
			for (int k=0; k<7; k++) {
				if (k==1 || k==5)
					YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,degeneracy/12.0,M);
				else
					YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,degeneracy/6.0,M);
			}

		} else {
			for (int k=0; k<5; k++) {
				if (k==2)
					YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,degeneracy/3.0,M);
				else
					YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,degeneracy/6.0,M);
			}
		}
	} else YplusisCtimesAtimesB(phi,Gf,Gb,degeneracy,M);
}


void LGrad2::AddPhiS(Real* phi,Real* Gf,Real* Gb, Real* G1, Real norm, int Markov, int M){
	cout << "composition not (yet) implemented for alias in LGrad2 " << endl;
}

void LGrad2::Initiate(Real* G,Real* Gz,int Markov, int M){
	if (Markov==2) {
		if (lattice_type == hexagonal) {
			for (int k=0; k<7; k++) Cp(G+k*M,Gz,M);
		} else {
			for (int k=0; k<5; k++) Cp(G+k*M,Gz,M);
		}
	} else Cp(G,Gz,M);
}

void LGrad2::Terminate(Real* Gz,Real* G,int Markov, int M){
	Real one=1.0;
	if (Markov==2) {
		Zero(Gz,M);
		if (lattice_type == hexagonal) {
			Add(Gz,G,M); Add(Gz,G+2*M,M); Add(Gz,G+3*M,M); Add(Gz,G+4*M,M); Add(Gz,G+6*M,M);
			Norm(Gz,2*one,M);
			Add(Gz,G+M,M); Add(Gz,G+5*M,M);
			Norm(Gz,1.0/12.0*one,M);
		} else {
			for (int k=0; k<5; k++) Add(Gz,G+k*M,M);
			Norm(Gz,1.0/6.0*one,M);
		}
	} else Cp(Gz,G,M);
}

bool LGrad2:: PutMask(int* MASK,vector<int>px,vector<int>py,vector<int>pz,int R){
if (debug) cout <<"PutMask in LGrad2 " << endl;
	bool success=true;
	int length =px.size();
	int X,Y;
	if (length > 1) {
		cout <<"In two gradient system, we can have just one particle: we found " <<length <<"particles. " << endl;
		return false;
	}
	for (int i =0; i<length; i++) {
		int xx,yy;
		xx=px[i]; yy=py[i];
		if (xx !=0) {
			cout <<"In two gradients system, we expect the particle at the central axis" << endl;
			return false;
		}
		for (int x=1; x<xx+R+1; x++)
		for (int y=yy-R; y<yy+R+1; y++){
			if ((xx-x)*(xx-x)+(yy-y)*(yy-y) <=R*R) {
				X=x; Y=y;
				if (y<1) {cout << "particle too close to y=0 boundary " << endl;
					return false;
				}
				if (x>MX) {
					cout <<"in two gradient system, particle should be smaller than size of system in radial direction" << endl;
					return false;
				}
				if (y>MY) { cout <<"particle too close to y upperbound " << endl;
					return false;
				}
				MASK[P(X,Y)]++;
			}
		}
	}

	return success;
}

Real LGrad2::DphiDt(Real* g, Real* B_phitot, Real* phiA, Real* phiB, Real* alphaA, Real* alphaB, Real B_A, Real B_B) {
	cout <<"LGrad2: DphiDt not implemented yet " << endl;
	return 0;
}
