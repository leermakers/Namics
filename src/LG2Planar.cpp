#include <iostream>
#include <string>
#include "lattice.h"
#include "LG2Planar.h"

LG2Planar::LG2Planar(vector<Input*> In_,string name_): LGrad2(In_,name_) {
	JY=1;
}

LG2Planar::~LG2Planar() {
if (debug) cout <<"LG2Planar destructor " << endl;
}


void LG2Planar:: ComputeLambdas() {

	if (fjc==1) {
		for (int i=0; i<M; i++) L[i]=1;
	}

	if (fjc==2) {
		for (int i=0; i<M; i++) L[i]=1.0/fjc;
	}
}


void LG2Planar::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
if (debug) cout <<" Side in LG2Planar " << endl;
	if (ignore_sites) {
		Cp(X_side,X,M); return;
	}
	Zero(X_side,M);set_bounds(X);

	if (fcc_sites) {
		Add(X_side,X,     M);
		Add(X_side+1,X,   M-1);
		Add(X_side,X+1,   M-1);
		Add(X_side+JX,X,  M-JX);
		Add(X_side,X+JX,  M-JX);
		Add(X_side+JX+1,X,M-JX-1);
		Add(X_side+JX,X+1,M-JX);
		Add(X_side+1,X+JX,M-JX);
		Add(X_side,X+JX+1,M-JX-1);
		Norm(X_side,1.0/9.0,M);
	} else {
		if (fjc==1) {
			if (!stencil_full) {
				if (lattice_type==simple_cubic) { //6 point stencil (voorheen 9-punts
					Add(X_side+JX,X,M-JX);
					Add(X_side,X+JX,M-JX);
					Add(X_side+1,X,M-1);
					Add(X_side,X+1,M-1);
					Norm(X_side,1.0/2.0,M);
					Add(X_side,X,M);
					Norm(X_side,1.0/3.0,M);
				} else { //Johan's method
					Add(X_side+JX,X,   M-JX);
					Add(X_side,   X+JX,M-JX);
					Add(X_side+JY,X   ,M-JY);
					Add(X_side,   X+JY,M-JY);
					Add(X_side,   X   ,M);
					Norm(X_side,2.0,M);

					remove_bounds(X);
					set_bounds_x(X,-1);
					Add(X_side+JX,X+JY,M-JX-JY);
					Add(X_side+JY,X+JX,M-JX-JY);

					Norm(X_side,1.0/12.0,M);

						//YplusisCtimesX(X_side,X,    12.0/48.0*one,M);
						//YplusisCtimesX(X_side+1,X,   6.0/48.0*one,M-1);
						//YplusisCtimesX(X_side,X+1,   6.0/48.0*one,M-1);
						//YplusisCtimesX(X_side+JX,X,  6.0/48.0*one,M-JX);
						//YplusisCtimesX(X_side,X+JX,  6.0/48.0*one,M-JX);
						//YplusisCtimesX(X_side+JX+1,X,3.0/48.0*one,M-JX-1);
						//YplusisCtimesX(X_side+JX,X+1,3.0/48.0*one,M-JX);
						//YplusisCtimesX(X_side+1,X+JX,3.0/48.0*one,M-JX);
						//YplusisCtimesX(X_side,X+JX+1,3.0/48.0*one,M-JX-1);
					}
			} else {
				if (lattice_type==simple_cubic) {  //9punts stencil, voorheen 6pnts
					//cout <<"not fully tested" << endl;
					Real C1=16.0/36.0;
					Real C2=4.0/36.0;
					Real C3=1.0/36.0;
					YplusisCtimesX(X_side,X,    C1,M);
					YplusisCtimesX(X_side+1,X,   C2,M-1);
					YplusisCtimesX(X_side,X+1,   C2,M-1);
					YplusisCtimesX(X_side+JX,X,  C2,M-JX);
					YplusisCtimesX(X_side,X+JX,  C2,M-JX);
					YplusisCtimesX(X_side+JX+1,X,C3,M-JX-1);
					YplusisCtimesX(X_side+JX,X+1,C3,M-JX);
					YplusisCtimesX(X_side+1,X+JX,C3,M-JX);
					YplusisCtimesX(X_side,X+JX+1,C3,M-JX-1);

				} else {
						//hexagonal //9 point stencil
						//Add(gs,gs_1,M);
					Real Two=2.0;
					Real C=1.0/16.0;
					YplusisCtimesX(X_side,X,Two,M);
					Add(X_side+JX,X,   M-JX);
					Add(X_side,   X+JX,M-JX);
					Add(X_side+JY,X,   M-JY);
					Add(X_side,   X+JY,M-JY);
					Norm(X_side,Two,M);
					Add(X_side+JX+JY,X,      M-JX-JY);
					Add(X_side,      X+JX+JY,M-JX-JY);
					Add(X_side+JX,   X+JY,   M-JX-JY);
					Add(X_side+JY,   X+JX,   M-JX-JY);

					Norm(X_side,C,M);
				}
			}

		} else {  //fjc>1
			for (int block=0; block<3; block++){
				int bk;
				int a,b;
				for (int x=-fjc; x<fjc+1; x++) for (int y=-fjc; y<fjc+1; y++) {
					bk=a=b=0;
					if (x==-fjc || x==fjc) bk++;
					if (y==-fjc || y==fjc) bk++;
					if (bk==block) {
						if (x<0) a =-x*JX; else b=x*JX;
						if (y<0) a -=y*JY; else b+=y*JY;
						Add(X_side+a,X+b,M-a-b);
					}
				}
				if (block !=2) Norm(X_side,2.0,M); else Norm(X_side,1.0/(4.0*(FJC-2)*FJC+1),M);
			}
		}
	}
}

void LG2Planar::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	if (!stencil_full) {
		if (lattice_type==hexagonal) {
			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;

			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			Real *gx0=gs,   *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M, *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			Zero(gs,12*M);
			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			//set_bounds_x(gz0,gz11,0,0); set_bounds_x(gz1,gz10,0,0);set_bounds_x(gz2,gz9,0,0); set_bounds_x(gz3,gz8,0,0); set_bounds_x(gz4,gz7,0,0); set_bounds_x(gz5,gz6,0,0);
			set_bounds_x(gz0,gz11,0); set_bounds_x(gz1,gz10,0);set_bounds_x(gz2,gz9,0); set_bounds_x(gz3,gz8,0); set_bounds_x(gz4,gz7,0); set_bounds_x(gz5,gz6,0);

			YplusisCtimesX(gx0+JX,gz0,P[0],M-JX);
			YplusisCtimesX(gx0+JX,gz1,P[0],M-JX);
			YplusisCtimesX(gx0+JX,gz2,P[0],M-JX);
			YplusisCtimesX(gx0+JX,gz3,P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz4,P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz5,P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz6,P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz7,P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz8,P[1],M-JX);

			YplusisCtimesX(gx11,gz3+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz4+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz5+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz6+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz7+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz8+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz9+JX,P[0],M-JX);
			YplusisCtimesX(gx11,gz10+JX,P[0],M-JX);
			YplusisCtimesX(gx11,gz11+JX,P[0],M-JX);

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			//set_bounds_y(gz2,gz9,0,0);  set_bounds_y(gz3,gz8,0,0); set_bounds_y(gz4,gz7,0,0); set_bounds_y(gz0,gz11,0,0); set_bounds_y(gz1,gz10,0,0); set_bounds_y(gz5,gz6,0,0);
			//set_bounds_z(gz2,gz9,0);  set_bounds_z(gz3,gz8,0); set_bounds_z(gz4,gz7,0); set_bounds_z(gz0,gz11,0); set_bounds_z(gz1,gz10,0); set_bounds_z(gz5,gz6,0,0);

			YplusisCtimesX(gx3,gz0,P[1],M);
			YplusisCtimesX(gx3,gz1,P[1],M);
			YplusisCtimesX(gx3,gz2,P[1],M);
			YplusisCtimesX(gx3,gz3,P[0],M);
			YplusisCtimesX(gx3,gz4,P[0],M);
			YplusisCtimesX(gx3,gz5,P[0],M);
			YplusisCtimesX(gx3,gz9,P[1],M);
			YplusisCtimesX(gx3,gz10,P[1],M);
			YplusisCtimesX(gx3,gz11,P[1],M);

			YplusisCtimesX(gx8,gz0,P[1],M);
			YplusisCtimesX(gx8,gz1,P[1],M);
			YplusisCtimesX(gx8,gz2,P[1],M);
			YplusisCtimesX(gx8,gz6,P[0],M);
			YplusisCtimesX(gx8,gz7,P[0],M);
			YplusisCtimesX(gx8,gz8,P[0],M);
			YplusisCtimesX(gx8,gz9,P[1],M);
			YplusisCtimesX(gx8,gz10,P[1],M);
			YplusisCtimesX(gx8,gz11,P[1],M);

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			//set_bounds_z(gz1,gz10,0,0); set_bounds_z(gz4,gz7,0,0); set_bounds_z(gz5,gz6,0,0); set_bounds_z(gz0,gz11,0,0); set_bounds_z(gz2,gz9,0,0); set_bounds_z(gz3,gz8,0,0);
			set_bounds_y(gz1,gz10,0); set_bounds_y(gz4,gz7,0); set_bounds_y(gz5,gz6,0); set_bounds_y(gz0,gz11,0); set_bounds_y(gz2,gz9,0); set_bounds_y(gz3,gz8,0);

			YplusisCtimesX(gx5+JY,gz0,P[1],M-JY);
			YplusisCtimesX(gx5+JY,gz1,P[1],M-JY);
			YplusisCtimesX(gx5+JY,gz2,P[1],M-JY);
			YplusisCtimesX(gx5+JY,gz3,P[0],M-JY);
			YplusisCtimesX(gx5+JY,gz4,P[0],M-JY);
			YplusisCtimesX(gx5+JY,gz5,P[0],M-JY);
			YplusisCtimesX(gx5+JY,gz9,P[1],M-JY);
			YplusisCtimesX(gx5+JY,gz10,P[1],M-JY);
			YplusisCtimesX(gx5+JY,gz11,P[1],M-JY);

			YplusisCtimesX(gx6,gz0+JY,P[1],M-JY);
			YplusisCtimesX(gx6,gz1+JY,P[1],M-JY);
			YplusisCtimesX(gx6,gz2+JY,P[1],M-JY);
			YplusisCtimesX(gx6,gz6+JY,P[0],M-JY);
			YplusisCtimesX(gx6,gz7+JY,P[0],M-JY);
			YplusisCtimesX(gx6,gz8+JY,P[0],M-JY);
			YplusisCtimesX(gx6,gz9+JY,P[1],M-JY);
			YplusisCtimesX(gx6,gz10+JY,P[1],M-JY);
			YplusisCtimesX(gx6,gz11+JY,P[1],M-JY);

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			//set_bounds_x(gz0,gz11,0,-1); set_bounds_x(gz1,gz10,0,-1);set_bounds_x(gz2,gz9,0,-1); set_bounds_x(gz3,gz8,0,-1); set_bounds_x(gz4,gz7,0,-1); set_bounds_x(gz5,gz6,0,-1);
			set_bounds_x(gz0,gz11,-1); set_bounds_x(gz1,gz10,-1);set_bounds_x(gz2,gz9,-1); set_bounds_x(gz3,gz8,-1); set_bounds_x(gz4,gz7,-1); set_bounds_x(gz5,gz6,-1);

			YplusisCtimesX(gx1+JX,gz0+JY,P[0],M-JX-JY);
			YplusisCtimesX(gx1+JX,gz1+JY,P[0],M-JX-JY);
			YplusisCtimesX(gx1+JX,gz2+JY,P[0],M-JX-JY);
			YplusisCtimesX(gx1+JX,gz3+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx1+JX,gz4+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx1+JX,gz5+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx1+JX,gz6+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx1+JX,gz7+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx1+JX,gz8+JY,P[1],M-JX-JY);

			YplusisCtimesX(gx10+JY,gz3+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx10+JY,gz4+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx10+JY,gz5+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx10+JY,gz6+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx10+JY,gz7+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx10+JY,gz8+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx10+JY,gz9+JX,P[0],M-JX-JY);
			YplusisCtimesX(gx10+JY,gz10+JX,P[0],M-JX-JY);
			YplusisCtimesX(gx10+JY,gz11+JX,P[0],M-JX-JY);

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			//set_bounds_y(gz2,gz9,-1,0);  set_bounds_y(gz3,gz8,-1,0); set_bounds_y(gz4,gz7,-1,0); set_bounds_y(gz0,gz11,-1,0); set_bounds_y(gz1,gz10,-1,0); set_bounds_y(gz5,gz6,-1,0);
			//set_bounds_z(gz2,gz9,-1);  set_bounds_z(gz3,gz8,-1); set_bounds_z(gz4,gz7,-1); set_bounds_z(gz0,gz11,-1); set_bounds_z(gz1,gz10,-1); set_bounds_z(gz5,gz6,-1);
			set_bounds_x(gz0,gz11,0); set_bounds_x(gz1,gz10,0);set_bounds_x(gz2,gz9,0); set_bounds_x(gz3,gz8,0); set_bounds_x(gz4,gz7,0); set_bounds_x(gz5,gz6,0);

			YplusisCtimesX(gx2+JX,gz0,P[0],M-JX);
			YplusisCtimesX(gx2+JX,gz1,P[0],M-JX);
			YplusisCtimesX(gx2+JX,gz2,P[0],M-JX);
			YplusisCtimesX(gx2+JX,gz3,P[1],M-JX);
			YplusisCtimesX(gx2+JX,gz4,P[1],M-JX);
			YplusisCtimesX(gx2+JX,gz5,P[1],M-JX);
			YplusisCtimesX(gx2+JX,gz6,P[1],M-JX);
			YplusisCtimesX(gx2+JX,gz7,P[1],M-JX);
			YplusisCtimesX(gx2+JX,gz8,P[1],M-JX);

			YplusisCtimesX(gx9,gz3+JX,P[1],M-JX);
			YplusisCtimesX(gx9,gz4+JX,P[1],M-JX);
			YplusisCtimesX(gx9,gz5+JX,P[1],M-JX);
			YplusisCtimesX(gx9,gz6+JX,P[1],M-JX);
			YplusisCtimesX(gx9,gz7+JX,P[1],M-JX);
			YplusisCtimesX(gx9,gz8+JX,P[1],M-JX);
			YplusisCtimesX(gx9,gz9+JX,P[0],M-JX);
			YplusisCtimesX(gx9,gz10+JX,P[0],M-JX);
			YplusisCtimesX(gx9,gz11+JX,P[0],M-JX);

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			//set_bounds_z(gz1,gz10,0,-1); set_bounds_z(gz4,gz7,0,-1); set_bounds_z(gz5,gz6,0,-1); set_bounds_z(gz0,gz11,0,-1); set_bounds_z(gz2,gz9,0,-1); set_bounds_z(gz3,gz8,0,-1);
			set_bounds_y(gz1,gz10,0); set_bounds_y(gz4,gz7,0); set_bounds_y(gz5,gz6,0); set_bounds_y(gz0,gz11,0); set_bounds_y(gz2,gz9,0); set_bounds_y(gz3,gz8,0);

			YplusisCtimesX(gx4,gz0+JY,P[1],M-JY);
			YplusisCtimesX(gx4,gz1+JY,P[1],M-JY);
			YplusisCtimesX(gx4,gz2+JY,P[1],M-JY);
			YplusisCtimesX(gx4,gz3+JY,P[0],M-JY);
			YplusisCtimesX(gx4,gz4+JY,P[0],M-JY);
			YplusisCtimesX(gx4,gz5+JY,P[0],M-JY);
			YplusisCtimesX(gx4,gz9+JY,P[1],M-JY);
			YplusisCtimesX(gx4,gz10+JY,P[1],M-JY);
			YplusisCtimesX(gx4,gz11+JY,P[1],M-JY);

			YplusisCtimesX(gx7+JY,gz0,P[1],M-JY);
			YplusisCtimesX(gx7+JY,gz1,P[1],M-JY);
			YplusisCtimesX(gx7+JY,gz2,P[1],M-JY);
			YplusisCtimesX(gx7+JY,gz6,P[0],M-JY);
			YplusisCtimesX(gx7+JY,gz7,P[0],M-JY);
			YplusisCtimesX(gx7+JY,gz8,P[0],M-JY);
			YplusisCtimesX(gx7+JY,gz9,P[1],M-JY);
			YplusisCtimesX(gx7+JY,gz10,P[1],M-JY);
			YplusisCtimesX(gx7+JY,gz11,P[1],M-JY);

			for (int k=0; k<12; k++) Times(gs+k*M,gs+k*M,g,M);

/*
			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;

			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			Real *gx0=gs,   *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M, *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			Zero(gs,12*M);
			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_x(gz0,gz11,0); set_bounds_x(gz1,gz10,0);set_bounds_x(gz2,gz9,0); set_bounds_x(gz3,gz8,0); set_bounds_x(gz4,gz7,0); set_bounds_x(gz5,gz6,0);

			YplusisCtimesX(gx0+JX,gz0,P[0],M-JX); //0 and 1 are equivalent
			YplusisCtimesX(gx0+JX,gz1,P[0],M-JX); //10 and 11 are equivalent
			YplusisCtimesX(gx0+JX,gz2,P[0],M-JX); //3 and 4 equivalent
			YplusisCtimesX(gx0+JX,gz3,P[1],M-JX); //7 and 8
			YplusisCtimesX(gx0+JX,gz4,P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz5,P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz6,P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz7,P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz8,P[1],M-JX);

			YplusisCtimesX(gx11,gz3+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz4+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz5+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz6+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz7+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz8+JX,P[1],M-JX);
			YplusisCtimesX(gx11,gz9+JX,P[0],M-JX);
			YplusisCtimesX(gx11,gz10+JX,P[0],M-JX);
			YplusisCtimesX(gx11,gz11+JX,P[0],M-JX);

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_y(gz2,gz9,0);  set_bounds_y(gz3,gz8,0); set_bounds_y(gz4,gz7,0); set_bounds_y(gz0,gz11,0); set_bounds_y(gz1,gz10,0); set_bounds_y(gz5,gz6,0);

			YplusisCtimesX(gx3+JY,gz0,P[1],M-JY);
			YplusisCtimesX(gx3+JY,gz1,P[1],M-JY);
			YplusisCtimesX(gx3+JY,gz2,P[1],M-JY);
			YplusisCtimesX(gx3+JY,gz3,P[0],M-JY);
			YplusisCtimesX(gx3+JY,gz4,P[0],M-JY);
			YplusisCtimesX(gx3+JY,gz5,P[0],M-JY);
			YplusisCtimesX(gx3+JY,gz9,P[1],M-JY);
			YplusisCtimesX(gx3+JY,gz10,P[1],M-JY);
			YplusisCtimesX(gx3+JY,gz11,P[1],M-JY);


			YplusisCtimesX(gx8,gz0+JY,P[1],M-JY);
			YplusisCtimesX(gx8,gz1+JY,P[1],M-JY);
			YplusisCtimesX(gx8,gz2+JY,P[1],M-JY);
			YplusisCtimesX(gx8,gz6+JY,P[0],M-JY);
			YplusisCtimesX(gx8,gz7+JY,P[0],M-JY);
			YplusisCtimesX(gx8,gz8+JY,P[0],M-JY);
			YplusisCtimesX(gx8,gz9+JY,P[1],M-JY);
			YplusisCtimesX(gx8,gz10+JY,P[1],M-JY);
			YplusisCtimesX(gx8,gz11+JY,P[1],M-JY);


			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);

			YplusisCtimesX(gx5,gz0,P[1],M);
			YplusisCtimesX(gx5,gz1,P[1],M);
			YplusisCtimesX(gx5,gz2,P[1],M);
			YplusisCtimesX(gx5,gz3,P[0],M);
			YplusisCtimesX(gx5,gz4,P[0],M);
			YplusisCtimesX(gx5,gz5,P[0],M);
			YplusisCtimesX(gx5,gz9,P[1],M);
			YplusisCtimesX(gx5,gz10,P[1],M);
			YplusisCtimesX(gx5,gz11,P[1],M);

			YplusisCtimesX(gx6,gz0,P[1],M);
			YplusisCtimesX(gx6,gz1,P[1],M);
			YplusisCtimesX(gx6,gz2,P[1],M);
			YplusisCtimesX(gx6,gz6,P[0],M);
			YplusisCtimesX(gx6,gz7,P[0],M);
			YplusisCtimesX(gx6,gz8,P[0],M);
			YplusisCtimesX(gx6,gz9,P[1],M);
			YplusisCtimesX(gx6,gz10,P[1],M);
			YplusisCtimesX(gx6,gz11,P[1],M);

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_x(gz0,gz11,0); set_bounds_x(gz1,gz10,0);set_bounds_x(gz2,gz9,0); set_bounds_x(gz3,gz8,0); set_bounds_x(gz4,gz7,0); set_bounds_x(gz5,gz6,0);

			YplusisCtimesX(gx1+JX,gz0,P[0],M-JX);
			YplusisCtimesX(gx1+JX,gz1,P[0],M-JX);
			YplusisCtimesX(gx1+JX,gz2,P[0],M-JX);
			YplusisCtimesX(gx1+JX,gz3,P[1],M-JX);
			YplusisCtimesX(gx1+JX,gz4,P[1],M-JX);
			YplusisCtimesX(gx1+JX,gz5,P[1],M-JX);
			YplusisCtimesX(gx1+JX,gz6,P[1],M-JX);
			YplusisCtimesX(gx1+JX,gz7,P[1],M-JX);
			YplusisCtimesX(gx1+JX,gz8,P[1],M-JX);

			YplusisCtimesX(gx10,gz3+JX,P[1],M-JX);
			YplusisCtimesX(gx10,gz4+JX,P[1],M-JX);
			YplusisCtimesX(gx10,gz5+JX,P[1],M-JX);
			YplusisCtimesX(gx10,gz6+JX,P[1],M-JX);
			YplusisCtimesX(gx10,gz7+JX,P[1],M-JX);
			YplusisCtimesX(gx10,gz8+JX,P[1],M-JX);
			YplusisCtimesX(gx10,gz9+JX,P[0],M-JX);
			YplusisCtimesX(gx10,gz10+JX,P[0],M-JX);
			YplusisCtimesX(gx10,gz11+JX,P[0],M-JX);

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_y(gz2,gz9,-1);  set_bounds_y(gz3,gz8,-1); set_bounds_y(gz4,gz7,-1); set_bounds_y(gz0,gz11,-1); set_bounds_y(gz1,gz10,-1); set_bounds_y(gz5,gz6,-1);

			YplusisCtimesX(gx2+JX,gz0+JY,P[0],M-JX-JY);
			YplusisCtimesX(gx2+JX,gz1+JY,P[0],M-JX-JY);
			YplusisCtimesX(gx2+JX,gz2+JY,P[0],M-JX-JY);
			YplusisCtimesX(gx2+JX,gz3+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx2+JX,gz4+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx2+JX,gz5+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx2+JX,gz6+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx2+JX,gz7+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx2+JX,gz8+JY,P[1],M-JX-JY);

			YplusisCtimesX(gx9+JY,gz3+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx9+JY,gz4+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx9+JY,gz5+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx9+JY,gz6+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx9+JY,gz7+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx9+JY,gz8+JX,P[1],M-JX-JY);
			YplusisCtimesX(gx9+JY,gz9+JX,P[0],M-JX-JY);
			YplusisCtimesX(gx9+JY,gz10+JX,P[0],M-JX-JY);
			YplusisCtimesX(gx9+JY,gz11+JX,P[0],M-JX-JY);

			remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);

			YplusisCtimesX(gx4+JY,gz0,P[1],M-JY);
			YplusisCtimesX(gx4+JY,gz1,P[1],M-JY);
			YplusisCtimesX(gx4+JY,gz2,P[1],M-JY);
			YplusisCtimesX(gx4+JY,gz3,P[0],M-JY);
			YplusisCtimesX(gx4+JY,gz4,P[0],M-JY);
			YplusisCtimesX(gx4+JY,gz5,P[0],M-JY);
			YplusisCtimesX(gx4+JY,gz9,P[1],M-JY);
			YplusisCtimesX(gx4+JY,gz10,P[1],M-JY);
			YplusisCtimesX(gx4+JY,gz11,P[1],M-JY);

			YplusisCtimesX(gx7,gz0+JY,P[1],M-JY);
			YplusisCtimesX(gx7,gz1+JY,P[1],M-JY);
			YplusisCtimesX(gx7,gz2+JY,P[1],M-JY);
			YplusisCtimesX(gx7,gz6+JY,P[0],M-JY);
			YplusisCtimesX(gx7,gz7+JY,P[0],M-JY);
			YplusisCtimesX(gx7,gz8+JY,P[0],M-JY);
			YplusisCtimesX(gx7,gz9+JY,P[1],M-JY);
			YplusisCtimesX(gx7,gz10+JY,P[1],M-JY);
			YplusisCtimesX(gx7,gz11+JY,P[1],M-JY);

			for (int k=0; k<12; k++) Times(gs+k*M,gs+k*M,g,M);
*/

		} else { //simple_cubic
			Real *gs=G+M*5*s_to;
			Real *gs_1=G+M*5*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M;
			set_bounds_x(gz0,gz4,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0); set_bounds_x(gz3,0);
			set_bounds_y(gz1,gz3,0); set_bounds_y(gz0,0); set_bounds_y(gz2,0); set_bounds_y(gz4,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M;
			Real *g=G1;

			Zero(gs,5*M);
			YplusisCtimesX(gx0+JX,gz0,  P[0],M-JX);
			YplusisCtimesX(gx0+JX,gz1,  P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz2,2*P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz3,  P[1],M-JX);

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

			YplusisCtimesX(gx4,gz1+JX,P[1],M-JX);
			YplusisCtimesX(gx4,gz2+JX,2*P[1],M-JX);
			YplusisCtimesX(gx4,gz3+JX,P[1],M-JX);
			YplusisCtimesX(gx4,gz4+JX,P[0],M-JX);

			for (int k=0; k<5; k++) Times(gs+k*M,gs+k*M,g,M);
		}
	} else {
		if (lattice_type == simple_cubic) { //size=2*FJC-1 here FJC=3
			if (fjc==1) {

			} else { //fjc==2

			/* geheugensteun
			0   x   x
			1   x   y
			2   x  -y
			3   x   0
			4   y   x
			5   y   y
			6   y   -x
			7   y  0
			8   0   y
			9   0   x
			10  0   0
			11  0  -x
			12  0  -y
			13 -y   0
			14 -y  x
			15 -y  -y
			16 -y -x
			17 -x   0
			18 -x  -y
			19 -x   y
			20 -x  -x
			*/
				int size = 2*FJC-1;
				Real *gs=G+M*size*s_to;
				Real *gs_1=G+M*size*s_from;
				Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M;
				Real            *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M;
				set_bounds_x(gz0,gz8,0);
				set_bounds_x(gz1,gz7,0);
				set_bounds_x(gz2,gz6,0);
				set_bounds_x(gz3,0);
				set_bounds_x(gz4,0);
				set_bounds_x(gz5,0);

				set_bounds_y(gz1,gz6,0);
				set_bounds_y(gz2,gz7,0);
				set_bounds_y(gz3,gz5,0);
				set_bounds_y(gz0,0);
				set_bounds_y(gz5,0);
				set_bounds_y(gz8,0);
				Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M;
				//Real *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M, *gx8=gs+8*M;
				Real *g=G1;

				Zero(gs,size*M);
				YplusisCtimesX(gx0+2*JX,gz0,  P[0],M-2*JX);
				YplusisCtimesX(gx0+2*JX,gz1,  P[1],M-2*JX);
				YplusisCtimesX(gx0+2*JX,gz2,  P[1],M-2*JX);
				YplusisCtimesX(gx0+2*JX,gz3,  P[2],M-2*JX);
				YplusisCtimesX(gx0+2*JX,gz4,  P[2],M-2*JX);
				YplusisCtimesX(gx0+2*JX,gz5,  P[2],M-2*JX);
				YplusisCtimesX(gx0+2*JX,gz6,  P[3],M-2*JX);
				YplusisCtimesX(gx0+2*JX,gz7,  P[3],M-2*JX);

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

				YplusisCtimesX(gx4,gz1+JX,P[1],M-JX);
				YplusisCtimesX(gx4,gz2+JX,2*P[1],M-JX);
				YplusisCtimesX(gx4,gz3+JX,P[1],M-JX);
				YplusisCtimesX(gx4,gz4+JX,P[0],M-JX);
				for (int k=0; k<size; k++) Times(gs+k*M,gs+k*M,g,M);
			}
		} else {//hexagonal // not correct because molecules can not 'diffuse' from one plane to the other.... ('behoud' wel segmenten....)
			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M;
			Real *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M;
			Real *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			set_bounds_x(gz0,gz3,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0);
			set_bounds_y(gz1,gz2,0); set_bounds_y(gz0,0); set_bounds_y(gz3,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M;
			Real *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M;
			Real *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			Zero(gs,12*M);
			YplusisCtimesX(gx0+JX,gz0,P[0],M-JX);
			YplusisCtimesX(gx0+JX,gz1,P[1],M-JX);
			YplusisCtimesX(gx0+JX,gz2,P[1],M-JX);

			YplusisCtimesX(gx1+JY,gz0,P[1],M-JY);
			YplusisCtimesX(gx1+JY,gz1,P[0],M-JY);
			YplusisCtimesX(gx1+JY,gz3,P[1],M-JY);

			YplusisCtimesX(gx2,gz0+JY,P[1],M-JY);
			YplusisCtimesX(gx2,gz2+JY,P[0],M-JY);
			YplusisCtimesX(gx2,gz3+JY,P[1],M-JY);

			YplusisCtimesX(gx3,gz1+JX,P[1],M-JX);
			YplusisCtimesX(gx3,gz2+JX,P[1],M-JX);
			YplusisCtimesX(gx3,gz3+JX,P[0],M-JX);

			YplusisCtimesX(gx4+JX,gz4,P[0],M-JX);
			YplusisCtimesX(gx4+JX,gz5,P[1],M-JX);
			YplusisCtimesX(gx4+JX,gz6,P[1],M-JX);

			YplusisCtimesX(gx5,gz4,P[1],M);
			YplusisCtimesX(gx5,gz6,P[0],M);
			YplusisCtimesX(gx5,gz7,P[1],M);

			YplusisCtimesX(gx6,gz4,P[1],M);
			YplusisCtimesX(gx6,gz5,P[0],M);
			YplusisCtimesX(gx6,gz7,P[1],M);

			YplusisCtimesX(gx7,gz5+JX,P[1],M-JX);
			YplusisCtimesX(gx7,gz6+JX,P[1],M-JX);
			YplusisCtimesX(gx7,gz7+JX,P[0],M-JX);


			YplusisCtimesX(gx8,gz8,P[0],M);
			YplusisCtimesX(gx8,gz9,P[1],M);
			YplusisCtimesX(gx8,gz10,P[1],M);

			YplusisCtimesX(gx9+JY,gz8,P[1],M-JY);
			YplusisCtimesX(gx9+JY,gz9,P[0],M-JY);
			YplusisCtimesX(gx9+JY,gz11,P[1],M-JY);

			YplusisCtimesX(gx10,gz8+JY,P[1],M-JY);
			YplusisCtimesX(gx10,gz10+JY,P[0],M-JY);
			YplusisCtimesX(gx10,gz11+JY,P[1],M-JY);

			YplusisCtimesX(gx11,gz9,P[1],M);
			YplusisCtimesX(gx11,gz10,P[1],M);
			YplusisCtimesX(gx11,gz11,P[0],M);

			for (int k=0; k<12; k++) Times(gs+k*M,gs+k*M,g,M);		}
	}
}
void LG2Planar::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	if (!stencil_full) {
		if (lattice_type==hexagonal) {
			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;

			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M, *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			Zero(gs,12*M);
			for (int k=0; k<12; k++) remove_bounds(gs_1+k*M);
			//remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);
			//remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
			set_bounds_x(gz0,gz11,0);

			YplusisCtimesX(gx3+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx4+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx5+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx6+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx7+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx8+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx9+JX,gz11,P[0],M-JX);
			YplusisCtimesX(gx10+JX,gz11,P[0],M-JX);
			YplusisCtimesX(gx11+JX,gz11,P[0],M-JX);

			YplusisCtimesX(gx0,gz0+JX,P[0],M-JX);
			YplusisCtimesX(gx1,gz0+JX,P[0],M-JX);
			YplusisCtimesX(gx2,gz0+JX,P[0],M-JX);
			YplusisCtimesX(gx3,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx4,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx5,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx6,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx7,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx8,gz0+JX,P[1],M-JX);

			remove_bounds(gz0);remove_bounds(gz11);
			//set_bounds_z(gz3,gz8,0);

			YplusisCtimesX(gx0,gz8,P[1],M);
			YplusisCtimesX(gx1,gz8,P[1],M);
			YplusisCtimesX(gx2,gz8,P[1],M);
			YplusisCtimesX(gx6,gz8,P[0],M);
			YplusisCtimesX(gx7,gz8,P[0],M);
			YplusisCtimesX(gx8,gz8,P[0],M);
			YplusisCtimesX(gx9,gz8,P[1],M);
			YplusisCtimesX(gx10,gz8,P[1],M);
			YplusisCtimesX(gx11,gz8,P[1],M);

			YplusisCtimesX(gx0,gz3,P[1],M);
			YplusisCtimesX(gx1,gz3,P[1],M);
			YplusisCtimesX(gx2,gz3,P[1],M);
			YplusisCtimesX(gx3,gz3,P[0],M);
			YplusisCtimesX(gx4,gz3,P[0],M);
			YplusisCtimesX(gx5,gz3,P[0],M);
			YplusisCtimesX(gx9,gz3,P[1],M);
			YplusisCtimesX(gx10,gz3,P[1],M);
			YplusisCtimesX(gx11,gz3,P[1],M);

			remove_bounds(gz3);remove_bounds(gz8);
			set_bounds_y(gz5,gz6,0);

			YplusisCtimesX(gx0+JY,gz6,P[1],M-JY);
			YplusisCtimesX(gx1+JY,gz6,P[1],M-JY);
			YplusisCtimesX(gx2+JY,gz6,P[1],M-JY);
			YplusisCtimesX(gx6+JY,gz6,P[0],M-JY);
			YplusisCtimesX(gx7+JY,gz6,P[0],M-JY);
			YplusisCtimesX(gx8+JY,gz6,P[0],M-JY);
			YplusisCtimesX(gx9+JY,gz6,P[1],M-JY);
			YplusisCtimesX(gx10+JY,gz6,P[1],M-JY);
			YplusisCtimesX(gx11+JY,gz6,P[1],M-JY);

			YplusisCtimesX(gx0,gz5+JY,P[1],M-JY);
			YplusisCtimesX(gx1,gz5+JY,P[1],M-JY);
			YplusisCtimesX(gx2,gz5+JY,P[1],M-JY);
			YplusisCtimesX(gx3,gz5+JY,P[0],M-JY);
			YplusisCtimesX(gx4,gz5+JY,P[0],M-JY);
			YplusisCtimesX(gx5,gz5+JY,P[0],M-JY);
			YplusisCtimesX(gx9,gz5+JY,P[1],M-JY);
			YplusisCtimesX(gx10,gz5+JY,P[1],M-JY);
			YplusisCtimesX(gx11,gz5+JY,P[1],M-JY);

			remove_bounds(gz5); remove_bounds(gz6);
			set_bounds_x(gz1,gz10,-1);

			YplusisCtimesX(gx3+JX,gz10+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx4+JX,gz10+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx5+JX,gz10+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx6+JX,gz10+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx7+JX,gz10+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx8+JX,gz10+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx9+JX,gz10+JY,P[0],M-JX-JY);
			YplusisCtimesX(gx10+JX,gz10+JY,P[0],M-JX-JY);
			YplusisCtimesX(gx11+JX,gz10+JY,P[0],M-JX-JY);

			YplusisCtimesX(gx0+JY,gz1+JX,P[0],M-JY-JX);
			YplusisCtimesX(gx1+JY,gz1+JX,P[0],M-JY-JX);
			YplusisCtimesX(gx2+JY,gz1+JX,P[0],M-JY-JX);
			YplusisCtimesX(gx3+JY,gz1+JX,P[1],M-JY-JX);
			YplusisCtimesX(gx4+JY,gz1+JX,P[1],M-JY-JX);
			YplusisCtimesX(gx5+JY,gz1+JX,P[1],M-JY-JX);
			YplusisCtimesX(gx6+JY,gz1+JX,P[1],M-JY-JX);
			YplusisCtimesX(gx7+JY,gz1+JX,P[1],M-JY-JX);
			YplusisCtimesX(gx8+JY,gz1+JX,P[1],M-JY-JX);

			remove_bounds(gz1);remove_bounds(gz10);
			set_bounds_x(gz2,gz9,0);

			YplusisCtimesX(gx3+JX,gz9,P[1],M-JX);
			YplusisCtimesX(gx4+JX,gz9,P[1],M-JX);
			YplusisCtimesX(gx5+JX,gz9,P[1],M-JX);
			YplusisCtimesX(gx6+JX,gz9,P[1],M-JX);
			YplusisCtimesX(gx7+JX,gz9,P[1],M-JX);
			YplusisCtimesX(gx8+JX,gz9,P[1],M-JX);
			YplusisCtimesX(gx9+JX,gz9,P[0],M-JX);
			YplusisCtimesX(gx10+JX,gz9,P[0],M-JX);
			YplusisCtimesX(gx11+JX,gz9,P[0],M-JX);

			YplusisCtimesX(gx0,gz2+JX,P[0],M-JX);
			YplusisCtimesX(gx1,gz2+JX,P[0],M-JX);
			YplusisCtimesX(gx2,gz2+JX,P[0],M-JX);
			YplusisCtimesX(gx3,gz2+JX,P[1],M-JX);
			YplusisCtimesX(gx4,gz2+JX,P[1],M-JX);
			YplusisCtimesX(gx5,gz2+JX,P[1],M-JX);
			YplusisCtimesX(gx6,gz2+JX,P[1],M-JX);
			YplusisCtimesX(gx7,gz2+JX,P[1],M-JX);
			YplusisCtimesX(gx8,gz2+JX,P[1],M-JX);

			remove_bounds(gz2);remove_bounds(gz9);
			set_bounds_y(gz4,gz7,0);

			YplusisCtimesX(gx0,gz7+JY,P[1],M-JY);
			YplusisCtimesX(gx1,gz7+JY,P[1],M-JY);
			YplusisCtimesX(gx2,gz7+JY,P[1],M-JY);
			YplusisCtimesX(gx6,gz7+JY,P[0],M-JY);
			YplusisCtimesX(gx7,gz7+JY,P[0],M-JY);
			YplusisCtimesX(gx8,gz7+JY,P[0],M-JY);
			YplusisCtimesX(gx9,gz7+JY,P[1],M-JY);
			YplusisCtimesX(gx10,gz7+JY,P[1],M-JY);
			YplusisCtimesX(gx11,gz7+JY,P[1],M-JY);

			YplusisCtimesX(gx0+JY,gz4,P[1],M-JY);
			YplusisCtimesX(gx1+JY,gz4,P[1],M-JY);
			YplusisCtimesX(gx2+JY,gz4,P[1],M-JY);
			YplusisCtimesX(gx3+JY,gz4,P[0],M-JY);
			YplusisCtimesX(gx4+JY,gz4,P[0],M-JY);
			YplusisCtimesX(gx5+JY,gz4,P[0],M-JY);
			YplusisCtimesX(gx9+JY,gz4,P[1],M-JY);
			YplusisCtimesX(gx10+JY,gz4,P[1],M-JY);
			YplusisCtimesX(gx11+JY,gz4,P[1],M-JY);

			for (int k=0; k<12; k++) Times(gs+k*M,gs+k*M,g,M);

/*
			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;

			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M, *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M, *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			Zero(gs,12*M);
			for (int k=0; k<12; k++) remove_bounds(gs_1+k*M);

			set_bounds_x(gz0,gz11,0);

			YplusisCtimesX(gx3+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx4+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx5+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx6+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx7+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx8+JX,gz11,P[1],M-JX);
			YplusisCtimesX(gx9+JX,gz11,P[0],M-JX);
			YplusisCtimesX(gx10+JX,gz11,P[0],M-JX);
			YplusisCtimesX(gx11+JX,gz11,P[0],M-JX);

			YplusisCtimesX(gx0,gz0+JX,P[0],M-JX);
			YplusisCtimesX(gx1,gz0+JX,P[0],M-JX);
			YplusisCtimesX(gx2,gz0+JX,P[0],M-JX);
			YplusisCtimesX(gx3,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx4,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx5,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx6,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx7,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx8,gz0+JX,P[1],M-JX);

			remove_bounds(gz0);remove_bounds(gz11);
			set_bounds_y(gz3,gz8,0);

			YplusisCtimesX(gx0+JY,gz8,P[1],M-JY);
			YplusisCtimesX(gx1+JY,gz8,P[1],M-JY);
			YplusisCtimesX(gx2+JY,gz8,P[1],M-JY);
			YplusisCtimesX(gx6+JY,gz8,P[0],M-JY);
			YplusisCtimesX(gx7+JY,gz8,P[0],M-JY);
			YplusisCtimesX(gx8+JY,gz8,P[0],M-JY);
			YplusisCtimesX(gx9+JY,gz8,P[1],M-JY);
			YplusisCtimesX(gx10+JY,gz8,P[1],M-JY);
			YplusisCtimesX(gx11+JY,gz8,P[1],M-JY);

			YplusisCtimesX(gx0,gz3+JY,P[1],M-JY);
			YplusisCtimesX(gx1,gz3+JY,P[1],M-JY);
			YplusisCtimesX(gx2,gz3+JY,P[1],M-JY);
			YplusisCtimesX(gx3,gz3+JY,P[0],M-JY);
			YplusisCtimesX(gx4,gz3+JY,P[0],M-JY);
			YplusisCtimesX(gx5,gz3+JY,P[0],M-JY);
			YplusisCtimesX(gx9,gz3+JY,P[1],M-JY);
			YplusisCtimesX(gx10,gz3+JY,P[1],M-JY);
			YplusisCtimesX(gx11,gz3+JY,P[1],M-JY);

			remove_bounds(gz3);remove_bounds(gz8);

			YplusisCtimesX(gx0,gz6,P[1],M);
			YplusisCtimesX(gx1,gz6,P[1],M);
			YplusisCtimesX(gx2,gz6,P[1],M);
			YplusisCtimesX(gx6,gz6,P[0],M);
			YplusisCtimesX(gx7,gz6,P[0],M);
			YplusisCtimesX(gx8,gz6,P[0],M);
			YplusisCtimesX(gx9,gz6,P[1],M);
			YplusisCtimesX(gx10,gz6,P[1],M);
			YplusisCtimesX(gx11,gz6,P[1],M);

			YplusisCtimesX(gx0,gz5,P[1],M);
			YplusisCtimesX(gx1,gz5,P[1],M);
			YplusisCtimesX(gx2,gz5,P[1],M);
			YplusisCtimesX(gx3,gz5,P[0],M);
			YplusisCtimesX(gx4,gz5,P[0],M);
			YplusisCtimesX(gx5,gz5,P[0],M);
			YplusisCtimesX(gx9,gz5,P[1],M);
			YplusisCtimesX(gx10,gz5,P[1],M);
			YplusisCtimesX(gx11,gz5,P[1],M);

			remove_bounds(gz5); remove_bounds(gz6);
			set_bounds_x(gz1,gz10,0);

			YplusisCtimesX(gx3+JX,gz10,P[1],M-JX);
			YplusisCtimesX(gx4+JX,gz10,P[1],M-JX);
			YplusisCtimesX(gx5+JX,gz10,P[1],M-JX);
			YplusisCtimesX(gx6+JX,gz10,P[1],M-JX);
			YplusisCtimesX(gx7+JX,gz10,P[1],M-JX);
			YplusisCtimesX(gx8+JX,gz10,P[1],M-JX);
			YplusisCtimesX(gx9+JX,gz10,P[0],M-JX);
			YplusisCtimesX(gx10+JX,gz10,P[0],M-JX);
			YplusisCtimesX(gx11+JX,gz10,P[0],M-JX);

			YplusisCtimesX(gx0,gz1+JX,P[0],M-JX);
			YplusisCtimesX(gx1,gz1+JX,P[0],M-JX);
			YplusisCtimesX(gx2,gz1+JX,P[0],M-JX);
			YplusisCtimesX(gx3,gz1+JX,P[1],M-JX);
			YplusisCtimesX(gx4,gz1+JX,P[1],M-JX);
			YplusisCtimesX(gx5,gz1+JX,P[1],M-JX);
			YplusisCtimesX(gx6,gz1+JX,P[1],M-JX);
			YplusisCtimesX(gx7,gz1+JX,P[1],M-JX);
			YplusisCtimesX(gx8,gz1+JX,P[1],M-JX);

			remove_bounds(gz1);remove_bounds(gz10);
			set_bounds_y(gz2,gz9,-1);

			YplusisCtimesX(gx3+JX,gz9+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx4+JX,gz9+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx5+JX,gz9+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx6+JX,gz9+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx7+JX,gz9+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx8+JX,gz9+JY,P[1],M-JX-JY);
			YplusisCtimesX(gx9+JX,gz9+JY,P[0],M-JX-JY);
			YplusisCtimesX(gx10+JX,gz9+JY,P[0],M-JX-JY);
			YplusisCtimesX(gx11+JX,gz9+JY,P[0],M-JX-JY);

			YplusisCtimesX(gx0+JY,gz2+JX,P[0],M-JY-JX);
			YplusisCtimesX(gx1+JY,gz2+JX,P[0],M-JY-JX);
			YplusisCtimesX(gx2+JY,gz2+JX,P[0],M-JY-JX);
			YplusisCtimesX(gx3+JY,gz2+JX,P[1],M-JY-JX);
			YplusisCtimesX(gx4+JY,gz2+JX,P[1],M-JY-JX);
			YplusisCtimesX(gx5+JY,gz2+JX,P[1],M-JY-JX);
			YplusisCtimesX(gx6+JY,gz2+JX,P[1],M-JY-JX);
			YplusisCtimesX(gx7+JY,gz2+JX,P[1],M-JY-JX);
			YplusisCtimesX(gx8+JY,gz2+JX,P[1],M-JY-JX);

			remove_bounds(gz2);remove_bounds(gz9);

			YplusisCtimesX(gx0+JY,gz7,P[1],M-JY);
			YplusisCtimesX(gx1+JY,gz7,P[1],M-JY);
			YplusisCtimesX(gx2+JY,gz7,P[1],M-JY);
			YplusisCtimesX(gx6+JY,gz7,P[0],M-JY);
			YplusisCtimesX(gx7+JY,gz7,P[0],M-JY);
			YplusisCtimesX(gx8+JY,gz7,P[0],M-JY);
			YplusisCtimesX(gx9+JY,gz7,P[1],M-JY);
			YplusisCtimesX(gx10+JY,gz7,P[1],M-JY);
			YplusisCtimesX(gx11+JY,gz7,P[1],M-JY);

			YplusisCtimesX(gx0,gz4+JY,P[1],M-JY);
			YplusisCtimesX(gx1,gz4+JY,P[1],M-JY);
			YplusisCtimesX(gx2,gz4+JY,P[1],M-JY);
			YplusisCtimesX(gx3,gz4+JY,P[0],M-JY);
			YplusisCtimesX(gx4,gz4+JY,P[0],M-JY);
			YplusisCtimesX(gx5,gz4+JY,P[0],M-JY);
			YplusisCtimesX(gx9,gz4+JY,P[1],M-JY);
			YplusisCtimesX(gx10,gz4+JY,P[1],M-JY);
			YplusisCtimesX(gx11,gz4+JY,P[1],M-JY);

			for (int k=0; k<12; k++) Times(gs+k*M,gs+k*M,g,M);
*/
		} else {
			Real *gs=G+M*5*s_to;
			Real *gs_1=G+M*5*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M,*gz2=gs_1+2*M, *gz3=gs_1+3*M, *gz4=gs_1+4*M;
			set_bounds_x(gz0,gz4,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0); set_bounds_x(gz3,0);
			set_bounds_y(gz1,gz3,0); set_bounds_y(gz0,0); set_bounds_y(gz2,0); set_bounds_y(gz4,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M, *gx4=gs+4*M;
			Real *g=G1;

			Zero(gs,5*M);
			YplusisCtimesX(gx1+JX,gz4,P[1],M-JX);
			YplusisCtimesX(gx2+JX,gz4,P[1],M-JX);
			YplusisCtimesX(gx3+JX,gz4,P[1],M-JX);
			YplusisCtimesX(gx4+JX,gz4,P[0],M-JX);

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

			YplusisCtimesX(gx0,gz0+JX,P[0],M-JX);
			YplusisCtimesX(gx1,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx2,gz0+JX,P[1],M-JX);
			YplusisCtimesX(gx3,gz0+JX,P[1],M-JX);

			for (int k=0; k<5; k++) Times(gs+k*M,gs+k*M,g,M);		}

	} else {
		if (lattice_type==simple_cubic) {
			cout <<"simple_cubic, markov 2, planar, stencil_full, not implemented " << endl;
		} else { //hexagonal //not finished and also not physically acceptable....there are three sub-lattices and there is no cross-over from one to the other.
			Real *gs=G+M*12*s_to;
			Real *gs_1=G+M*12*s_from;
			Real *gz0=gs_1, *gz1=gs_1+M, *gz2=gs_1+2*M, *gz3=gs_1+3*M;
			Real *gz4=gs_1+4*M, *gz5=gs_1+5*M, *gz6=gs_1+6*M, *gz7=gs_1+7*M;
			Real *gz8=gs_1+8*M, *gz9=gs_1+9*M, *gz10=gs_1+10*M, *gz11=gs_1+11*M;
			set_bounds_x(gz0,gz3,0); set_bounds_x(gz1,0); set_bounds_x(gz2,0);
			set_bounds_y(gz1,gz2,0); set_bounds_y(gz0,0); set_bounds_y(gz3,0);
			Real *gx0=gs, *gx1=gs+M, *gx2=gs+2*M, *gx3=gs+3*M;
			Real *gx4=gs+4*M, *gx5=gs+5*M, *gx6=gs+6*M, *gx7=gs+7*M;
			Real *gx8=gs+8*M, *gx9=gs+9*M, *gx10=gs+10*M, *gx11=gs+11*M;
			Real *g=G1;

			Zero(gs,12*M);
			YplusisCtimesX(gx0+JY,gz2,   P[1],  M-JY);
			YplusisCtimesX(gx0,   gz1+JY,P[1],  M-JY);
			YplusisCtimesX(gx0,   gz0+JX,P[0],  M-JX);

			YplusisCtimesX(gx1+JX,gz3,   P[1],  M-JX);
			YplusisCtimesX(gx1,   gz1+JY,P[0],  M-JY);
			YplusisCtimesX(gx1,   gz0+JX,P[1],  M-JX);

			YplusisCtimesX(gx2+JX,gz3,   P[1],  M-JX);
			YplusisCtimesX(gx2+JY,gz2,   P[0],  M-JY);
			YplusisCtimesX(gx2,   gz0+JX,P[1],  M-JX);

			YplusisCtimesX(gx3+JX, gz3,   P[0],  M-JX);
			YplusisCtimesX(gx3+JY, gz2,   P[1],  M-JY);
			YplusisCtimesX(gx3,    gz1+JY,P[1],  M-JY);


			YplusisCtimesX(gx4,   gz6,   P[1],  M);
			YplusisCtimesX(gx4,   gz5,   P[1],  M);
			YplusisCtimesX(gx4,   gz4+JX,P[0],  M-JX);

			YplusisCtimesX(gx5+JX,gz7,   P[1],  M-JX);
			YplusisCtimesX(gx5,   gz5,   P[0],  M);
			YplusisCtimesX(gx5,   gz4+JX,P[1],  M-JX);

			YplusisCtimesX(gx6+JX,gz7,   P[1],  M-JX);
			YplusisCtimesX(gx6,   gz6,   P[0],  M);
			YplusisCtimesX(gx6,   gz4+JX,P[1],  M-JX);

			YplusisCtimesX(gx7+JX, gz7,   P[0],  M-JX);
			YplusisCtimesX(gx7,    gz6,   P[1],  M);
			YplusisCtimesX(gx7,    gz5,   P[1],  M);


			YplusisCtimesX(gx8+JY,gz10,   P[1],  M-JY);
			YplusisCtimesX(gx8,   gz9+JY, P[1],  M-JY);
			YplusisCtimesX(gx8,   gz8,    P[0],  M);

			YplusisCtimesX(gx9,   gz11,   P[1],  M);
			YplusisCtimesX(gx9,   gz9+JY, P[0],  M-JY);
			YplusisCtimesX(gx9,   gz8,    P[1],  M);

			YplusisCtimesX(gx10,   gz11,   P[1],  M);
			YplusisCtimesX(gx10+JY,gz10,   P[0],  M-JY);
			YplusisCtimesX(gx10,   gz8,    P[1],  M);

			YplusisCtimesX(gx11,    gz11,  P[0],  M);
			YplusisCtimesX(gx11+JY, gz10,  P[1],  M-JY);
			YplusisCtimesX(gx11,    gz9+JY,P[1],  M-JY);
			for (int k=0; k<12; k++) Times(gs+k*M,gs+k*M,g,M);		}
	}
}


void LG2Planar::propagate(Real *G, Real *G1, int s_from, int s_to,int M) { //this procedure should function on simple cubic lattice.
if (debug) cout <<" propagate in LGrad2 " << endl;
	Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	Zero(gs,M); set_bounds(gs_1);
	if (fjc==1) {
		if (!stencil_full) {
			if (lattice_type==simple_cubic) { //9 point stencil..
				Add(gs+JX,gs_1,M-JX);
				Add(gs,gs_1+JX,M-JX);
				Add(gs+JY,gs_1,M-JY);
				Add(gs,gs_1+JY,M-JY);
				Norm(gs,1.0/2.0,M);
				Add(gs,gs_1,M);
				Norm(gs,1.0/3.0,M);
				Times(gs,gs,G1,M);


			} else { //hexagonal Johan's method //kept for nostalgic reasons

				Add(gs+JX,gs_1,   M-JX);
				Add(gs,   gs_1+JX,M-JX);
				Add(gs+JY,gs_1   ,M-JY);
				Add(gs,   gs_1+JY,M-JY);
				Add(gs,   gs_1   ,M);
				Norm(gs,2.0,M);

				remove_bounds(gs_1);
				set_bounds_x(gs_1,-1);
				Add(gs+JX,gs_1+JY,M-JX-JY);
				Add(gs+JY,gs_1+JX,M-JX-JY);

				Norm(gs,1.0/12.0,M);
				Times(gs,gs,G1,M);
				}
		} else {
			if (lattice_type==simple_cubic) {
				Real C1=16.0/36.0;
				Real C2=4.0/36.0;
				Real C3=1.0/36.0;
				YplusisCtimesX(gs,gs_1,    C1,M);
				YplusisCtimesX(gs+1,gs_1,   C2,M-1);
				YplusisCtimesX(gs,gs_1+1,   C2,M-1);
				YplusisCtimesX(gs+JX,gs_1,  C2,M-JX);
				YplusisCtimesX(gs,gs_1+JX,  C2,M-JX);
				YplusisCtimesX(gs+JX+1,gs_1,C3,M-JX-1);
				YplusisCtimesX(gs+JX,gs_1+1,C3,M-JX);
				YplusisCtimesX(gs+1,gs_1+JX,C3,M-JX);
				YplusisCtimesX(gs,gs_1+JX+1,C3,M-JX-1);
				Times(gs,gs,G1,M);
			} else { //hexagonal //9 point stencil

				//Add(gs,gs_1,M);
				Real Two=2.0;
				Real C=1.0/16.0;
				YplusisCtimesX(gs,gs_1,Two,M);
				Add(gs+JX,gs_1,   M-JX);
				Add(gs,   gs_1+JX,M-JX);
				Add(gs+JY,gs_1,   M-JY);
				Add(gs,   gs_1+JY,M-JY);
				Norm(gs,Two,M);
				Add(gs+JX+JY,gs_1,      M-JX-JY);
				Add(gs,      gs_1+JX+JY,M-JX-JY);
				Add(gs+JX,   gs_1+JY,   M-JX-JY);
				Add(gs+JY,   gs_1+JX,   M-JX-JY);

				Norm(gs,C,M);
				Times(gs,gs,G1,M);
			}
		}
	} else {
		for (int block=0; block<3; block++){
			Real Two=2.0;
			int bk;
			int a,b;
			for (int x=-fjc; x<fjc+1; x++) for (int y=-fjc; y<fjc+1; y++) {
				bk=a=b=0;
				if (x==-fjc || x==fjc) bk++;
				if (y==-fjc || y==fjc) bk++;
				if (bk==block) {
					if (x<0) a =-x*JX; else b=x*JX;
					if (y<0) a -=y*JY; else b+=y*JY;
					Add(gs+a,gs_1+b,M-a-b);
				}
			}
			if (block !=2) Norm(gs,Two,M); else Norm(gs,1.0/(4.0*(FJC-2)*FJC+1),M);
		}
		Times(gs,gs,G1,M);
	}
}


void LG2Planar::UpdateEE(Real* EE, Real* psi, Real* E) {
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_M_bounds(psi);
	Zero(EE,M);
	Real Exmin,Explus,Eymin,Eyplus;
	int x,y,z;

	pf = pf/2.0;
	for (x=fjc; x<MX+fjc; x++) {
		for (y=fjc; y<MY+fjc; y++) {
			z=x*JX+y;
			Exmin=psi[z]-psi[z-JX];
			Exmin*=Exmin;
			Explus=psi[z]-psi[z+JX];
			Explus*=Explus;
			Eymin=psi[z]-psi[z-1];
			Eymin*=Eymin;
			Eyplus=psi[z]-psi[z+1];
			Eyplus*=Eyplus;
			EE[x*JX+y]=pf*(Exmin+Explus+Eymin+Eyplus);
		}
	}
}


void LG2Planar::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, int* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	int x,y,i;

	Real epsXplus, epsXmin, epsYplus,epsYmin;
	//set_M_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);
	Real ax,ay;

	if (!fixedPsi0) {
		C=C*2.0/fjc/fjc;
		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++) {
				i=x*JX+y;
				epsXmin=eps[i]+eps[i-JX];
				epsXplus=eps[i]+eps[i+JX];
				epsYmin=eps[i]+eps[i-1];
				epsYplus=eps[i]+eps[i+1];
				if (x==fjc) ax=psi[i-JX]; else ax=X[i-JX]; //upwind
				if (y==fjc) ay=psi[i-1]; else ay=X[i-1]; //upwind
				//X[i]= (C*q[i]+epsXmin*psi[i-JX]+epsXplus*psi[i+JX]+epsYmin*psi[i-1]+epsYplus*psi[i+1])/(epsXmin+epsXplus+epsYmin+epsYplus);
				X[i]= (C*q[i]+epsXmin*ax+epsXplus*psi[i+JX]+epsYmin*ay+epsYplus*psi[i+1])/(epsXmin+epsXplus+epsYmin+epsYplus);
			}
		}
		//Cp(psi,X,M);
		YisAminB(g,g,X,M);
  	} else { //fixedPsi0 is true
		for (x=fjc; x<MX+fjc; x++) {
			for (y=fjc; y<MY+fjc; y++){
				if (Mask[x*JX+y] == 0)
				X[x*JX+y]=0.25*(psi[(x-1)*JX+y]+psi[(x+1)*JX+y]
			        +psi[x*JX+y-1]  +psi[x*JX+y+1])
				 +0.5*q[x*JX+y]*C/eps[x*JX+y];
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


void LG2Planar::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, int* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
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

bool LG2Planar:: PutMask(int* MASK,vector<int>px,vector<int>py,vector<int>pz,int R){
	bool success=false;
	cout <<"PutMask does not make sence in planar 2 gradient system " << endl;
	return success;
}

Real LG2Planar::DphiDt(Real* g, Real* B_phitot, Real* phiA, Real* phiB, Real* alphaA, Real* alphaB, Real B_A, Real B_B) {
	cout <<"LG2Planar : DphiDt not implemented yet " << endl;
	return 0;
}
