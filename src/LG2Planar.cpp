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
	Real one=1.0;
	if (ignore_sites) {
		Cp(X_side,X,M); return;
	}
	Zero(X_side,M);set_bounds(X);

	if (fcc_sites) {
		YplusisCtimesX(X_side,X,     1.0/9.0*one,M);
		YplusisCtimesX(X_side+1,X,   1.0/9.0*one,M-1);
		YplusisCtimesX(X_side,X+1,   1.0/9.0*one,M-1);
		YplusisCtimesX(X_side+JX,X,  1.0/9.0*one,M-JX);
		YplusisCtimesX(X_side,X+JX,  1.0/9.0*one,M-JX);
		YplusisCtimesX(X_side+JX+1,X,1.0/9.0*one,M-JX-1);
		YplusisCtimesX(X_side+JX,X+1,1.0/9.0*one,M-JX);
		YplusisCtimesX(X_side+1,X+JX,1.0/9.0*one,M-JX);
		YplusisCtimesX(X_side,X+JX+1,1.0/9.0*one,M-JX-1);
	} else {
		if (fjc==1) {
			if (stencil_full) {
				if (lattice_type==simple_cubic) { //9 point stencil
					YplusisCtimesX(X_side,X,    16.0/36.0*one,M);
					YplusisCtimesX(X_side+1,X,   4.0/36.0*one,M-1);
					YplusisCtimesX(X_side,X+1,   4.0/36.0*one,M-1);
					YplusisCtimesX(X_side+JX,X,  4.0/36.0*one,M-JX);
					YplusisCtimesX(X_side,X+JX,  4.0/36.0*one,M-JX);
					YplusisCtimesX(X_side+JX+1,X,1.0/36.0*one,M-JX-1);
					YplusisCtimesX(X_side+JX,X+1,1.0/36.0*one,M-JX);
					YplusisCtimesX(X_side+1,X+JX,1.0/36.0*one,M-JX);
					YplusisCtimesX(X_side,X+JX+1,1.0/36.0*one,M-JX-1);
				} else {
					YplusisCtimesX(X_side,X,    12.0/48.0*one,M);
					YplusisCtimesX(X_side+1,X,   6.0/48.0*one,M-1);
					YplusisCtimesX(X_side,X+1,   6.0/48.0*one,M-1);
					YplusisCtimesX(X_side+JX,X,  6.0/48.0*one,M-JX);
					YplusisCtimesX(X_side,X+JX,  6.0/48.0*one,M-JX);
					YplusisCtimesX(X_side+JX+1,X,3.0/48.0*one,M-JX-1);
					YplusisCtimesX(X_side+JX,X+1,3.0/48.0*one,M-JX);
					YplusisCtimesX(X_side+1,X+JX,3.0/48.0*one,M-JX);
					YplusisCtimesX(X_side,X+JX+1,3.0/48.0*one,M-JX-1);
				}
			} else {
				if (lattice_type==simple_cubic) {//classical
					Add(X_side+JX,X,M-JX);
					Add(X_side,X+JX,M-JX);
					Add(X_side+1,X,M-1);
					Add(X_side,X+1,M-1);
					Norm(X_side,1.0/2.0*one,M);
					Add(X_side,X,M);
					Norm(X_side,1.0/3.0*one,M);
				} else { //not fully tested...
					cout <<" in side fractions the bc are not fully tested yet..." << endl;
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

					Norm(X_side,1.0/12.0*one,M);
				}
			}
		}
		if (fjc==2) {
			Add(X_side,X,M);
			Add(X_side+JX,X,M-JX);
			Add(X_side,X+JX,M-JX);
			Add(X_side+1,X,M-1);
			Add(X_side,X+1,M-1);
			Add(X_side+JX+1,X,M-JX-1);
			Add(X_side,X+JX+1,M-JX-1);
			Add(X_side+1,X+JX,M-JX);
			Add(X_side+JX,X+1,M-JX);
			Norm(X_side,2.0,M);
			Add(X_side+2*JX,X,M-2*JX);
			Add(X_side,X+2*JX,M-2*JX);
			Add(X_side+2*JX,X+1,M-2*JX);
			Add(X_side+2*JX+1,X,M-2*JX-1);
			Add(X_side+1,X+2*JX,M-2*JX);
			Add(X_side,X+2*JX+1,M-2*JX-1);
			Add(X_side+JX+2,X,M-JX-2);
			Add(X_side+JX,X+2,M-JX);
			Add(X_side,X+JX+2,M-JX-2);
			Add(X_side+2,X+JX,M-JX);
			Add(X_side+2,X,M-2);
			Add(X_side,X+2,M-2);
			Norm(X_side,2.0,M);
			Add(X_side+2*JX+2,X,M-2*JX-2);
			Add(X_side,X+2*JX+2,M-2*JX-2);
			Add(X_side+2*JX,X+2,M-2*JX);
			Add(X_side+2,X+2*JX,M-2*JX);
			Norm(X_side,1.0/64.0,M);
		}
	}
}


void LG2Planar::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
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
		set_bounds_x(gz0,gz6,0);set_bounds_x(gz1,gz5,0);set_bounds_x(gz2,gz4,0);set_bounds_x(gz3,0);

		YplusisCtimesX(gx0+JX,gz0,2*P[0],M-JX);
		YplusisCtimesX(gx0+JX,gz1,  P[0],M-JX);
		YplusisCtimesX(gx0+JX,gz2,2*P[1],M-JX);
		YplusisCtimesX(gx0+JX,gz3,2*P[1],M-JX);
		YplusisCtimesX(gx0+JX,gz4,2*P[1],M-JX);

		YplusisCtimesX(gx6,gz2+JX,2*P[1],M-JX);
		YplusisCtimesX(gx6,gz3+JX,2*P[1],M-JX);
		YplusisCtimesX(gx6,gz4+JX,2*P[1],M-JX);
		YplusisCtimesX(gx6,gz5+JX,  P[0],M-JX);
		YplusisCtimesX(gx6,gz6+JX,2*P[0],M-JX);

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
		YplusisCtimesX(gx1+JX,gz0+JY,2*P[0],M-JX-JY);
		YplusisCtimesX(gx1+JX,gz1+JY,  P[0],M-JX-JY);
		YplusisCtimesX(gx1+JX,gz2+JY,2*P[1],M-JX-JY);
		YplusisCtimesX(gx1+JX,gz3+JY,2*P[1],M-JX-JY);
		YplusisCtimesX(gx1+JX,gz4+JY,2*P[1],M-JX-JY);

		remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);
		set_bounds_x(gz6,-1);set_bounds_x(gz5,-1);set_bounds_x(gz2,-1);set_bounds_x(gz4,-1);set_bounds_x(gz3,-1);
		YplusisCtimesX(gx5+JY,gz2+JX,2*P[1],M-JX-JY);
		YplusisCtimesX(gx5+JY,gz3+JX,2*P[1],M-JX-JY);
		YplusisCtimesX(gx5+JY,gz4+JX,2*P[1],M-JX-JY);
		YplusisCtimesX(gx5+JY,gz5+JX,  P[0],M-JX-JY);
		YplusisCtimesX(gx5+JY,gz6+JX,2*P[0],M-JX-JY);

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
}
void LG2Planar::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
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

		YplusisCtimesX(gx2+JX,gz6,2*P[1],M-JX);
		YplusisCtimesX(gx3+JX,gz6,2*P[1],M-JX);
		YplusisCtimesX(gx4+JX,gz6,2*P[1],M-JX);
		YplusisCtimesX(gx5+JX,gz6,2*P[0],M-JX);
		YplusisCtimesX(gx6+JX,gz6,2*P[0],M-JX);

		YplusisCtimesX(gx0,gz0+JX,2*P[0],M-JX);
		YplusisCtimesX(gx1,gz0+JX,2*P[0],M-JX);
		YplusisCtimesX(gx2,gz0+JX,2*P[1],M-JX);
		YplusisCtimesX(gx3,gz0+JX,2*P[1],M-JX);
		YplusisCtimesX(gx4,gz0+JX,2*P[1],M-JX);

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

		set_bounds_y(gz1,gz5,-1); //waarom dit werkt is niet helemaal duidelijk. (zie asymmetry met 'forward')
		//set_bounds_y(gz1,-1);set_bounds_y(gz5,-1);

		YplusisCtimesX(gx2+JX,gz5+JY,P[1],M-JX-JY);
		YplusisCtimesX(gx3+JX,gz5+JY,P[1],M-JX-JY);
		YplusisCtimesX(gx4+JX,gz5+JY,P[1],M-JX-JY);
		YplusisCtimesX(gx5+JX,gz5+JY,P[0],M-JX-JY);
		YplusisCtimesX(gx6+JX,gz5+JY,P[0],M-JX-JY);

		remove_bounds(gz1);remove_bounds(gz5);
		set_bounds_x(gz1,gz5,-1);
		//set_bounds_x(gz1,-1);set_bounds_x(gz5,-1);

		YplusisCtimesX(gx0+JY,gz1+JX,P[0],M-JY-JX);
		YplusisCtimesX(gx1+JY,gz1+JX,P[0],M-JY-JX);
		YplusisCtimesX(gx2+JY,gz1+JX,P[1],M-JY-JX);
		YplusisCtimesX(gx3+JY,gz1+JX,P[1],M-JY-JX);
		YplusisCtimesX(gx4+JY,gz1+JX,P[1],M-JY-JX);

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


		for (int k=0; k<5; k++) Times(gs+k*M,gs+k*M,g,M);
	}
}


void LG2Planar::propagate(Real *G, Real *G1, int s_from, int s_to,int M) { //this procedure should function on simple cubic lattice.
if (debug) cout <<" propagate in LGrad2 " << endl;
	Real one=1.0;
	Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);

	Zero(gs,M); set_bounds(gs_1);
	if (fjc==1) {
		if (stencil_full) {
			if (lattice_type==simple_cubic) { //9 point stencil
				YplusisCtimesX(gs,gs_1,    16.0/36.0*one,M);
				YplusisCtimesX(gs+1,gs_1,   4.0/36.0*one,M-1);
				YplusisCtimesX(gs,gs_1+1,   4.0/36.0*one,M-1);
				YplusisCtimesX(gs+JX,gs_1,  4.0/36.0*one,M-JX);
				YplusisCtimesX(gs,gs_1+JX,  4.0/36.0*one,M-JX);
				YplusisCtimesX(gs+JX+1,gs_1,1.0/36.0*one,M-JX-1);
				YplusisCtimesX(gs+JX,gs_1+1,1.0/36.0*one,M-JX);
				YplusisCtimesX(gs+1,gs_1+JX,1.0/36.0*one,M-JX);
				YplusisCtimesX(gs,gs_1+JX+1,1.0/36.0*one,M-JX-1);
				Times(gs,gs,G1,M);
			} else { //hexagonal //9 point stencil
				YplusisCtimesX(gs,gs_1,    12.0/48.0*one,M);
				YplusisCtimesX(gs+1,gs_1,   6.0/48.0*one,M-1);
				YplusisCtimesX(gs,gs_1+1,   6.0/48.0*one,M-1);
				YplusisCtimesX(gs+JX,gs_1,  6.0/48.0*one,M-JX);
				YplusisCtimesX(gs,gs_1+JX,  6.0/48.0*one,M-JX);
				YplusisCtimesX(gs+JX+1,gs_1,3.0/48.0*one,M-JX-1);
				YplusisCtimesX(gs+JX,gs_1+1,3.0/48.0*one,M-JX);
				YplusisCtimesX(gs+1,gs_1+JX,3.0/48.0*one,M-JX);
				YplusisCtimesX(gs,gs_1+JX+1,3.0/48.0*one,M-JX-1);
				Times(gs,gs,G1,M);
			}
		} else { // classical!
			if (lattice_type==simple_cubic) {
				Add(gs+JX,gs_1,M-JX);
				Add(gs,gs_1+JX,M-JX);
				Add(gs+JY,gs_1,M-1);
				Add(gs,gs_1+JY,M-1);
				Norm(gs,1.0/2.0*one,M);
				Add(gs,gs_1,M);
				Norm(gs,1.0/3.0*one,M);
				Times(gs,gs,G1,M);
			} else { //hexagonal Johan's method
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
		}
	}
	if (fjc==2) { //25 point stencil only fjc==2 implemented....
		 //lattice_type = hexagonal
		Add(gs,gs_1,M);
		Add(gs+JX,gs_1,M-JX);
		Add(gs,gs_1+JX,M-JX);
		Add(gs+1,gs_1,M-1);
		Add(gs,gs_1+1,M-1);
		Add(gs+JX+1,gs_1,M-JX-1);
		Add(gs,gs_1+JX+1,M-JX-1);
		Add(gs+1,gs_1+JX,M-JX);
		Add(gs+JX,gs_1+1,M-JX);
		Norm(gs,2.0,M);
		Add(gs+2*JX,gs_1,M-2*JX);
		Add(gs,gs_1+2*JX,M-2*JX);
		Add(gs+2*JX,gs_1+1,M-2*JX);
		Add(gs+2*JX+1,gs_1,M-2*JX-1);
		Add(gs+1,gs_1+2*JX,M-2*JX);
		Add(gs,gs_1+2*JX+1,M-2*JX-1);
		Add(gs+JX+2,gs_1,M-JX-2);
		Add(gs+JX,gs_1+2,M-JX);
		Add(gs,gs_1+JX+2,M-JX-2);
		Add(gs+2,gs_1+JX,M-JX);
		Add(gs+2,gs_1,M-2);
		Add(gs,gs_1+2,M-2);
		Norm(gs,2.0,M);
		Add(gs+2*JX+2,gs_1,M-2*JX-2);
		Add(gs,gs_1+2*JX+2,M-2*JX-2);
		Add(gs+2*JX,gs_1+2,M-2*JX);
		Add(gs+2,gs_1+2*JX,M-2*JX);
		Norm(gs,1.0/64.0,M);
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
