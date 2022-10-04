#include <iostream>
#include <string>
#include "lattice.h"
#include "LGrad3.h"

LGrad3::LGrad3(vector<Input*> In_,string name_): Lattice(In_,name_) {}

LGrad3::~LGrad3() {
if (debug) cout <<"LGrad3 destructor " << endl;
}

void LGrad3:: ComputeLambdas() {
}

bool LGrad3::PutM() {
if (debug) cout << "PutM in LGrad3 " << endl;
	bool success=true;
	volume = MX*MY*MZ;
	JX=(MZ+2*fjc)*(MY+2*fjc); JY=MZ+2*fjc; JZ=1; M = (MX+2*fjc)*(MY+2*fjc)*(MZ+2*fjc);

	Accesible_volume=volume;
	return success;
}

void LGrad3::TimesL(Real* X){
if (debug) cout << "TimesL in LGrad3 " << endl;
}

void LGrad3::DivL(Real* X){
if (debug) cout << "DivL in LGrad3 " << endl;
}

Real LGrad3:: Moment(Real* X,Real Xb, int n) {
if (debug) cout << "Moment in LGrad3 " << endl;
	Real Result=0;
	//cout <<"Moment analysis not implemented in LGrad3 " << endl;
	return Result/fjc;
}

Real LGrad3::WeightedSum(Real* X){
if (debug) cout << "weighted sum in LGrad3 " << endl;
	Real sum{0};
	remove_bounds(X);
	Sum(sum,X,M);
	return sum;
}

void LGrad3::vtk(string filename, Real* X, string id,bool writebounds) {
if (debug) cout << "vtk in LGrad3 " << endl;
	FILE *fp;
	int i;
	fp = fopen(filename.c_str(),"w+");
	fprintf(fp,"# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\n",MZ,MY,MX);

	if (writebounds) {
		fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",(MX+2*fjc)*(MY+2*fjc)*(MZ+2*fjc));
	} else {
		fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",MX*MY*MZ);
	}
	fprintf(fp,"SCALARS %s double\nLOOKUP_TABLE default\n",id.c_str());
#ifdef LongReal
	if (writebounds) for(i=0; i<M; i++) fprintf(fp,"%Lf\n",X[i]);
	else {
		for (int x=fjc; x<MX+fjc; x++)
		for (int y=fjc; y<MY+fjc; y++)
		for (int z=fjc; z<MZ+fjc; z++)
		fprintf(fp,"%Le\n",X[P(x,y,z)]);
	}
#else
	if (writebounds) for(i=0; i<M; i++) fprintf(fp,"%f\n",X[i]);
	else {
		for (int x=fjc; x<MX+fjc; x++)
		for (int y=fjc; y<MY+fjc; y++)
		for (int z=fjc; z<MZ+fjc; z++)
		fprintf(fp,"%e\n",X[P(x,y,z)]);
	}
#endif
	fclose(fp);
}

void LGrad3::PutProfiles(FILE* pf,vector<Real*> X,bool writebounds){
if (debug) cout <<"PutProfiles in LGrad3 " << endl;
	Real one=1.0;
	int x,y,z,i;
	int length=X.size();
	int a;
	if (writebounds) a=0; else a = fjc;
	for (x=a; x<MX+2*fjc-a; x++)
	for (y=a; y<MY+2*fjc-a; y++)
	for (z=a; z<MZ+2*fjc-a; z++) {
#ifdef LongReal
		fprintf(pf,"%Le\t%Le\t%Le\t",one*(x-fjc+1)/fjc-0.5/fjc,one*(y-fjc+1)/fjc-0.5/fjc,one*(z-fjc+1)/fjc-0.5/fjc);
		for (i=0; i<length; i++) fprintf(pf,"%.20Le\t",X[i][P(x,y,z)]);
#else
		fprintf(pf,"%e\t%e\t%e\t",one*(x-fjc+1)/fjc-0.5/fjc,one*(y-fjc+1)/fjc-0.5/fjc,one*(z-fjc+1)/fjc-0.5/fjc);
		for (i=0; i<length; i++) fprintf(pf,"%.20e\t",X[i][P(x,y,z)]);

#endif
		fprintf(pf,"\n");
	}
}

void LGrad3::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
if (debug) cout <<" Side in LGrad3 " << endl;
	if (ignore_sites) {
		Cp(X_side,X,M); return;
	}
	Zero(X_side,M);//set_bounds(X);

	if (stencil_full) {
		Add(X_side+JX,X   ,M-JX);
		Add(X_side   ,X+JX,M-JX);
		Add(X_side+JY,X   ,M-JY);
		Add(X_side   ,X+JY,M-JY);
		Add(X_side+1 ,X   ,M-1);
		Add(X_side   ,X+1 ,M-1);

		if (lattice_type == simple_cubic) {
			Norm(X_side,4.0,M);
		} else {
			Norm(X_side,2.0,M);
		}
		Add(X_side+JX+JY, X,       M-JX-JY);
		Add(X_side,       X+JX+JY, M-JX-JY);
		Add(X_side+JY,    X+JX,    M-JY-JX);
		Add(X_side+JX,    X+JY,    M-JY-JX);
		Add(X_side+JX+1,  X,       M-JX-1);
		Add(X_side,       X+JX+1,  M-JX-1);
		Add(X_side+JX,    X+1,     M-JX);
		Add(X_side+1,     X+JX,    M-JX);
		Add(X_side+JY+1,  X,       M-JY-1);
		Add(X_side,       X+JY+1,  M-JX-1);
		Add(X_side+JY,    X+1,     M-JY);
		Add(X_side+1,     X+JY,    M-JY);
		if (lattice_type == simple_cubic) {
			Norm(X_side,4.0,M);
		} else {
			Norm(X_side,2.0,M);
		}
		Add(X_side+JX+JY+1,  X,	    M-JX-JY-1);
		Add(X_side,          X+JX+JY+1, M-JX-JY-1);
		Add(X_side+JX+JY,    X+1,       M-JX-JY-1);
		Add(X_side+1,        X+JX+JY,   M-JX-JY-1);
		Add(X_side+JX+1,     X+JY,      M-JX-JY-1);
		Add(X_side+JY,       X+JX+1,    M-JX-JY-1);
		Add(X_side+JY+1,     X+JX,      M-JX-JY-1);
		Add(X_side+JX,       X+JY+1,    M-JX-JY-1);
		if (lattice_type == simple_cubic) {
			Norm(X_side,1.0/152.0,M);
		} else {
			Norm(X_side,1.0/56.0,M);
		}
	} else {
		if (lattice_type==simple_cubic) {
			Add(X_side+JX,X   ,M-JX);
			Add(X_side   ,X+JX,M-JX);
			Add(X_side+JY,X   ,M-JY);
			Add(X_side   ,X+JY,M-JY);
			Add(X_side+JZ,X   ,M-JZ);
			Add(X_side   ,X+JZ,M-JZ);
	 		Norm(X_side,1.0/6.0,M);
		} else { //hexagonal
			if (fjc==1) {
				Add(X_side+JX ,X,    M-JX);
				Add(X_side    ,X+JX ,M-JX);
				Add(X_side+JX ,X+JY ,M-JX-JY);
				Add(X_side    ,X+JY ,M-JY);
				Add(X_side+JY ,X    ,M-JY);
				Add(X_side+JY ,X+JX ,M-JX-JY);
				Add(X_side+JZ ,X+JY ,M-JY-JZ);
				Add(X_side+JZ ,X    ,M-JZ);
				Add(X_side+JZ ,X+JX ,M-JX-JZ);
				Add(X_side+JX ,X+JZ ,M-JX-JZ);
				Add(X_side    ,X+JZ ,M-JZ);
				Add(X_side+JY ,X+JZ ,M-JY-JZ);
				Norm(X_side,1.0/12.0,M);
			} else { //fjc==2
				Add(X_side,        X,        M);   //1
				Add(X_side+JX,     X,        M-JX);//2
				Add(X_side,        X+JX,     M-JX);//3
				Add(X_side+JY,     X,        M-JY);//4
				Add(X_side,        X+JY,     M-JY);//5
				Add(X_side+1,      X,        M-1);//6
				Add(X_side,        X+1,      M-1);//7
				Add(X_side+JX+1,   X,        M-JX-1);//8
				Add(X_side,        X+JX+1,   M-JX-1);//9
				Add(X_side+1,      X+JX,     M-JX-1);//10
				Add(X_side+JX,     X+1,      M-JX-1);//11
				Add(X_side+JY+1,   X,        M-JY-1);//12
				Add(X_side,        X+JY+1,   M-JY-1);//13
				Add(X_side+1,      X+JY,     M-JY-1);//14
				Add(X_side+JY,     X+1,      M-JY-1);//15
				Add(X_side+JX+JY,  X,        M-JX-JY);//16
				Add(X_side,        X+JX+JY,  M-JX-JY);//17
				Add(X_side+JY,     X+JX,     M-JX-JY);//18
				Add(X_side+JX,     X+JY,     M-JX-JY);//19
				Add(X_side+JX+JY+1,X,        M-JX-JY-1);//20
				Add(X_side,        X+JX+JY+1,M-JX-JY-1);//21
				Add(X_side+JX+JY,  X+1,      M-JX-JY-1);//22
				Add(X_side+1,      X+JX+JY,  M-JX-JY-1);//23
				Add(X_side+JX+1,   X+JY,     M-JX-JY-1);//24
				Add(X_side+JY,     X+JX+1,   M-JX-JY-1);//25
				Add(X_side+JY+1,   X+JX,     M-JX-JY-1);//26
				Add(X_side+JX,     X+JY+1,   M-JX-JY-1);//27
				Norm(X_side,2.0,M);

				Add(X_side+2*JX,  X,          M-2*JX);//28
				Add(X_side,       X+2*JX,     M-2*JX);//29
				Add(X_side+2*JY,  X,          M-2*JY);//30
				Add(X_side,       X+2*JY,     M-2*JY);//31
				Add(X_side+2,     X,          M-2);//32
				Add(X_side,       X+2,        M-2);//33
				Add(X_side+2*JX+1, X,     	 M-2*JX-1);//34
				Add(X_side+2*JX,   X+1,       M-2*JX-1);//35
				Add(X_side+1,      X+2*JX,    M-2*JX-1);//36
				Add(X_side,        X+2*JX+1,  M-2*JX-1);//37
				Add(X_side+2*JY,   X+1,       M-2*JY-1);//38
				Add(X_side+2*JY+1, X,         M-2*JY-1);//39
				Add(X_side+1,      X+2*JY,    M-2*JY-1);//40
				Add(X_side,        X+2*JY+1,  M-2*JY-1);//41
				Add(X_side+JX+2,   X,         M-JX-2);//42
				Add(X_side+JX,     X+2,       M-JX-2);//43
				Add(X_side,        X+JX+2,    M-JX-2);//44
				Add(X_side+2,      X+JX,      M-JX-2);//45
				Add(X_side+JY+2,   X,         M-JY-2);//46
				Add(X_side+JY,     X+2,       M-JY-2);//47
				Add(X_side,        X+JY+2,    M-JY-2);//48
				Add(X_side+2,      X+JY,      M-JY-2);//49
				Add(X_side+2*JX,   X+JY,      M-2*JX-JY);//50
				Add(X_side+2*JX+JY,X,         M-2*JX-JY);//51
				Add(X_side+JY,     X+2*JX,    M-2*JX-JY);//52
				Add(X_side,        X+2*JX+JY, M-2*JX-JY);//53
				Add(X_side+2*JY,     X+JX,        M-2*JY-JX);//54
				Add(X_side+2*JY+JX,  X,           M-2*JY-JX);//55
				Add(X_side+JX,       X+2*JY,      M-2*JY-JX);//56
				Add(X_side,          X+2*JY+JX,   M-2*JY-JX);//57
				Add(X_side+2*JX+JY+1,X,           M-2*JX-JY-1);//58
				Add(X_side+2*JX+1,   X+JY,        M-2*JX-JY-1);//59
				Add(X_side+JY+1,     X+2*JX,      M-2*JX-JY-1);//60
				Add(X_side+1,        X+2*JX+JY,   M-2*JX-JY-1);//61
				Add(X_side+2*JX+JY,  X+1,    	 M-2*JX-JY-1);//62
				Add(X_side+2*JX,     X+JY+1,      M-2*JX-JY-1);//63
				Add(X_side+JY,       X+2*JX+1,    M-2*JX-JY-1);//64
				Add(X_side,          X+2*JX+JY+1, M-2*JX-JY-1);//65
				Add(X_side+JX+2*JY+1,X,           M-JX-2*JY-1);//66
				Add(X_side+JX+1,     X+2*JY,      M-JX-2*JY-1);//67
				Add(X_side+2*JY+1,   X+JX,        M-JX-2*JY-1);//68
				Add(X_side+1,        X+JX+2*JY,   M-JX-2*JY-1);//69
				Add(X_side+JX+2*JY,  X+1,         M-JX-2*JY-1);//70
				Add(X_side+JX,       X+2*JY+1,    M-JX-2*JY-1);//71
				Add(X_side+2*JY,     X+JX+1,      M-JX-2*JY-1);//72
				Add(X_side,          X+JX+2*JY+1, M-JX-2*JY-1);//73
				Add(X_side+JX+JY+2,  X,           M-JX-JY-2);//74
				Add(X_side+JX+2,     X+JY,        M-JX-JY-2);//75
				Add(X_side+JY+2,     X+JX,        M-JX-JY-2);//76
				Add(X_side+2,        X+JX+JY,     M-JX-JY-2);//77
				Add(X_side+JX+JY,    X+2,         M-JX-JY-2);//78
				Add(X_side+JX,       X+JY+2,      M-JX-JY-2);//79
				Add(X_side+JY,       X+JX+2,      M-JX-JY-2);//80
				Add(X_side,          X+JX+JY+2,   M-JX-JY-2);//81
				Norm(X_side,2.0,M);

				Add(X_side+2*JX+2*JY,X,           M-2*JX-2*JY);//82
				Add(X_side,          X+2*JX+2*JY, M-2*JX-2*JY);//83
				Add(X_side+2*JX,     X+2*JY,      M-2*JX-2*JY);//84
				Add(X_side+2*JY,     X+2*JX,      M-2*JX-2*JY);//85
				Add(X_side+2*JX+2*JY+1,X,           M-2*JX-2*JY-1);//86
				Add(X_side+1,          X+2*JX+2*JY, M-2*JX-2*JY-1);//87
				Add(X_side+2*JX+1,     X+2*JY,      M-2*JX-2*JY-1);//88
				Add(X_side+2*JY+1,     X+2*JX,      M-2*JX-2*JY-1);//89
				Add(X_side+2*JX+2*JY,X+1,           M-2*JX-2*JY-1);//90
				Add(X_side,          X+2*JX+2*JY+1, M-2*JX-2*JY-1);//91
				Add(X_side+2*JX,     X+2*JY+1,      M-2*JX-2*JY-1);//92
				Add(X_side+2*JY,     X+2*JX+1,      M-2*JX-2*JY-1);//93
				Add(X_side+2*JX+2,     X,           M-2*JX-2);//94
				Add(X_side+2*JX+JY+2,  X,           M-2*JX-JY-2);//95
				Add(X_side+2*JX+2,     X+JY,        M-2*JX-JY-2);//96
				Add(X_side+2,          X+2*JX,      M-2*JX-2);//97
				Add(X_side+JY+2,       X+2*JX,      M-2*JX-JY-2);//98
				Add(X_side+2,          X+2*JX+JY,   M-2*JX-JY-2);//99
				Add(X_side+2*JY+2,     X,           M-2*JY-2);//100
				Add(X_side+2*JY+JX+2,  X,           M-2*JY-JX-2);//101
				Add(X_side+2*JY+2,     X+JX,        M-2*JY-JX-2);//102
				Add(X_side+2,          X+2*JY,      M-2*JY-2);//103
				Add(X_side+JX+2,       X+2*JY,      M-2*JY-JX-2);//104
				Add(X_side+2,          X+2*JY+JX,   M-2*JY-JX-2);//105
				Add(X_side+2*JX,     X+2,           M-2*JX-2);//106
				Add(X_side+2*JX+JY,  X+2,           M-2*JX-JY-2);//107
				Add(X_side+2*JX,     X+JY+2,        M-2*JX-JY-2);//108
				Add(X_side,          X+2*JX+2,      M-2*JX-2);//109
				Add(X_side+JY,       X+2*JX+2,      M-2*JX-JY-2);//110
				Add(X_side,          X+2*JX+JY+2,   M-2*JX-JY-2);//111
				Add(X_side+2*JY,     X+2,           M-2*JY-2);//112
				Add(X_side+2*JY+JX,  X+2,           M-2*JY-JX-2);//113
				Add(X_side+2*JY,     X+JX+2,        M-2*JY-JX-2);//114
				Add(X_side,          X+2*JY+2,      M-2*JY-2);//115
				Add(X_side+JX,       X+2*JY+2,      M-2*JY-JX-2);//116
				Add(X_side,          X+2*JY+JX+2,   M-2*JY-JX-2);//117
				Norm(X_side,2.0,M);

				Add(X_side+2*JX+2*JY+2,X,           M-2*JX-2*JY-2);//118
				Add(X_side+2,          X+2*JX+2*JY, M-2*JX-2*JY-2);//119
				Add(X_side+2*JX+2,     X+2*JY,      M-2*JX-2*JY-2);//120
				Add(X_side+2*JY+2,     X+2*JX,      M-2*JX-2*JY-2);//121
				Add(X_side+2*JX+2*JY,X+2,           M-2*JX-2*JY-2);//122
				Add(X_side,          X+2*JX+2*JY+2, M-2*JX-2*JY-2);//123
				Add(X_side+2*JX,     X+2*JY+2,      M-2*JX-2*JY-2);//124
				Add(X_side+2*JY,     X+2*JX+2,      M-2*JX-2*JY-2);//125
				Norm(X_side,1.0/512.0,M);
			}
		}
	}
}

void LGrad3::propagateF(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	if (lattice_type == hexagonal ) {
		Real *gs=G+M*12*s_to;
		Real *gs_1=G+M*12*s_from;

		Real *gz0=gs_1;
		Real *gz1=gs_1+M;
		Real *gz2=gs_1+2*M;
		Real *gz3=gs_1+3*M;
		Real *gz4=gs_1+4*M;
		Real *gz5=gs_1+5*M;
		Real *gz6=gs_1+6*M;
		Real *gz7=gs_1+7*M;
		Real *gz8=gs_1+8*M;
		Real *gz9=gs_1+9*M;
		Real *gz10=gs_1+10*M;
		Real *gz11=gs_1+11*M;

		Real *gx0=gs;
		Real *gx1=gs+M;
		Real *gx2=gs+2*M;
		Real *gx3=gs+3*M;
		Real *gx4=gs+4*M;
		Real *gx5=gs+5*M;
		Real *gx6=gs+6*M;
		Real *gx7=gs+7*M;
		Real *gx8=gs+8*M;
		Real *gx9=gs+9*M;
		Real *gx10=gs+10*M;
		Real *gx11=gs+11*M;
		Real *g=G1;

		Zero(gs,12*M);
		remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
		set_bounds_x(gz0,gz11,0,0); set_bounds_x(gz1,gz10,0,0);set_bounds_x(gz2,gz9,0,0); set_bounds_x(gz3,gz8,0,0); set_bounds_x(gz4,gz7,0,0); set_bounds_x(gz5,gz6,0,0);

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
		set_bounds_y(gz2,gz9,0,0);  set_bounds_y(gz3,gz8,0,0); set_bounds_y(gz4,gz7,0,0); set_bounds_y(gz0,gz11,0,0); set_bounds_y(gz1,gz10,0,0); set_bounds_y(gz5,gz6,0,0);

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
		set_bounds_z(gz1,gz10,0,0); set_bounds_z(gz4,gz7,0,0); set_bounds_z(gz5,gz6,0,0); set_bounds_z(gz0,gz11,0,0); set_bounds_z(gz2,gz9,0,0); set_bounds_z(gz3,gz8,0,0);

		YplusisCtimesX(gx5+JZ,gz0,P[1],M-JZ);
		YplusisCtimesX(gx5+JZ,gz1,P[1],M-JZ);
		YplusisCtimesX(gx5+JZ,gz2,P[1],M-JZ);
		YplusisCtimesX(gx5+JZ,gz3,P[0],M-JZ);
		YplusisCtimesX(gx5+JZ,gz4,P[0],M-JZ);
		YplusisCtimesX(gx5+JZ,gz5,P[0],M-JZ);
		YplusisCtimesX(gx5+JZ,gz9,P[1],M-JZ);
		YplusisCtimesX(gx5+JZ,gz10,P[1],M-JZ);
		YplusisCtimesX(gx5+JZ,gz11,P[1],M-JZ);

		YplusisCtimesX(gx6,gz0+JZ,P[1],M-JZ);
		YplusisCtimesX(gx6,gz1+JZ,P[1],M-JZ);
		YplusisCtimesX(gx6,gz2+JZ,P[1],M-JZ);
		YplusisCtimesX(gx6,gz6+JZ,P[0],M-JZ);
		YplusisCtimesX(gx6,gz7+JZ,P[0],M-JZ);
		YplusisCtimesX(gx6,gz8+JZ,P[0],M-JZ);
		YplusisCtimesX(gx6,gz9+JZ,P[1],M-JZ);
		YplusisCtimesX(gx6,gz10+JZ,P[1],M-JZ);
		YplusisCtimesX(gx6,gz11+JZ,P[1],M-JZ);

		remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
		set_bounds_x(gz0,gz11,0,-1); set_bounds_x(gz1,gz10,0,-1);set_bounds_x(gz2,gz9,0,-1); set_bounds_x(gz3,gz8,0,-1); set_bounds_x(gz4,gz7,0,-1); set_bounds_x(gz5,gz6,0,-1);

		YplusisCtimesX(gx1+JX,gz0+JZ,P[0],M-JX-JZ);
		YplusisCtimesX(gx1+JX,gz1+JZ,P[0],M-JX-JZ);
		YplusisCtimesX(gx1+JX,gz2+JZ,P[0],M-JX-JZ);
		YplusisCtimesX(gx1+JX,gz3+JZ,P[1],M-JX-JZ);
		YplusisCtimesX(gx1+JX,gz4+JZ,P[1],M-JX-JZ);
		YplusisCtimesX(gx1+JX,gz5+JZ,P[1],M-JX-JZ);
		YplusisCtimesX(gx1+JX,gz6+JZ,P[1],M-JX-JZ);
		YplusisCtimesX(gx1+JX,gz7+JZ,P[1],M-JX-JZ);
		YplusisCtimesX(gx1+JX,gz8+JZ,P[1],M-JX-JZ);

		YplusisCtimesX(gx10+JZ,gz3+JX,P[1],M-JX-JZ);
		YplusisCtimesX(gx10+JZ,gz4+JX,P[1],M-JX-JZ);
		YplusisCtimesX(gx10+JZ,gz5+JX,P[1],M-JX-JZ);
		YplusisCtimesX(gx10+JZ,gz6+JX,P[1],M-JX-JZ);
		YplusisCtimesX(gx10+JZ,gz7+JX,P[1],M-JX-JZ);
		YplusisCtimesX(gx10+JZ,gz8+JX,P[1],M-JX-JZ);
		YplusisCtimesX(gx10+JZ,gz9+JX,P[0],M-JX-JZ);
		YplusisCtimesX(gx10+JZ,gz10+JX,P[0],M-JX-JZ);
		YplusisCtimesX(gx10+JZ,gz11+JX,P[0],M-JX-JZ);

		remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
		set_bounds_y(gz2,gz9,-1,0);  set_bounds_y(gz3,gz8,-1,0); set_bounds_y(gz4,gz7,-1,0); set_bounds_y(gz0,gz11,-1,0); set_bounds_y(gz1,gz10,-1,0); set_bounds_y(gz5,gz6,-1,0);

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
		set_bounds_z(gz1,gz10,0,-1); set_bounds_z(gz4,gz7,0,-1); set_bounds_z(gz5,gz6,0,-1); set_bounds_z(gz0,gz11,0,-1); set_bounds_z(gz2,gz9,0,-1); set_bounds_z(gz3,gz8,0,-1);


		YplusisCtimesX(gx4+JY,gz0+JZ,P[1],M-JY-JZ);
		YplusisCtimesX(gx4+JY,gz1+JZ,P[1],M-JY-JZ);
		YplusisCtimesX(gx4+JY,gz2+JZ,P[1],M-JY-JZ);
		YplusisCtimesX(gx4+JY,gz3+JZ,P[0],M-JY-JZ);
		YplusisCtimesX(gx4+JY,gz4+JZ,P[0],M-JY-JZ);
		YplusisCtimesX(gx4+JY,gz5+JZ,P[0],M-JY-JZ);
		YplusisCtimesX(gx4+JY,gz9+JZ,P[1],M-JY-JZ);
		YplusisCtimesX(gx4+JY,gz10+JZ,P[1],M-JY-JZ);
		YplusisCtimesX(gx4+JY,gz11+JZ,P[1],M-JY-JZ);

		YplusisCtimesX(gx7+JZ,gz0+JY,P[1],M-JY-JZ);
		YplusisCtimesX(gx7+JZ,gz1+JY,P[1],M-JY-JZ);
		YplusisCtimesX(gx7+JZ,gz2+JY,P[1],M-JY-JZ);
		YplusisCtimesX(gx7+JZ,gz6+JY,P[0],M-JY-JZ);
		YplusisCtimesX(gx7+JZ,gz7+JY,P[0],M-JY-JZ);
		YplusisCtimesX(gx7+JZ,gz8+JY,P[0],M-JY-JZ);
		YplusisCtimesX(gx7+JZ,gz9+JY,P[1],M-JY-JZ);
		YplusisCtimesX(gx7+JZ,gz10+JY,P[1],M-JY-JZ);
		YplusisCtimesX(gx7+JZ,gz11+JY,P[1],M-JY-JZ);

		for (int k=0; k<12; k++) Times(gs+k*M,gs+k*M,g,M);

	} else {
		Real *gs=G+M*6*s_to;
		Real *gs_1=G+M*6*s_from;

		Real *gz0=gs_1;
		Real *gz1=gs_1+M;
		Real *gz2=gs_1+2*M;
		Real *gz3=gs_1+3*M;
		Real *gz4=gs_1+4*M;
		Real *gz5=gs_1+5*M;
		set_bounds_x(gz0,gz5,0,0); set_bounds_x(gz1,0,0); set_bounds_x(gz2,0,0); set_bounds_x(gz3,0,0); set_bounds_x(gz4,0,0);
		set_bounds_y(gz1,gz4,0,0); set_bounds_y(gz0,0,0); set_bounds_y(gz2,0,0); set_bounds_y(gz3,0,0); set_bounds_y(gz5,0,0);
		set_bounds_z(gz2,gz3,0,0); set_bounds_z(gz0,0,0); set_bounds_z(gz1,0,0); set_bounds_z(gz4,0,0); set_bounds_z(gz5,0,0);
		Real *gx0=gs;
		Real *gx1=gs+M;
		Real *gx2=gs+2*M;
		Real *gx3=gs+3*M;
		Real *gx4=gs+4*M;
		Real *gx5=gs+5*M;
		Real *g=G1;

		Zero(gs,6*M);
		YplusisCtimesX(gx0+JX,gz0,P[0],M-JX);
		YplusisCtimesX(gx0+JX,gz1,P[1],M-JX);
		YplusisCtimesX(gx0+JX,gz2,P[1],M-JX);
		YplusisCtimesX(gx0+JX,gz3,P[1],M-JX);
		YplusisCtimesX(gx0+JX,gz4,P[1],M-JX);

		YplusisCtimesX(gx1+JY,gz0,P[1],M-JY);
		YplusisCtimesX(gx1+JY,gz1,P[0],M-JY);
		YplusisCtimesX(gx1+JY,gz2,P[1],M-JY);
		YplusisCtimesX(gx1+JY,gz3,P[1],M-JY);
		YplusisCtimesX(gx1+JY,gz5,P[1],M-JY);

		YplusisCtimesX(gx2+JZ,gz0,P[1],M-JZ);
		YplusisCtimesX(gx2+JZ,gz1,P[1],M-JZ);
		YplusisCtimesX(gx2+JZ,gz2,P[0],M-JZ);
		YplusisCtimesX(gx2+JZ,gz4,P[1],M-JZ);
		YplusisCtimesX(gx2+JZ,gz5,P[1],M-JZ);

		YplusisCtimesX(gx3,gz0+JZ,P[1],M-JZ);
		YplusisCtimesX(gx3,gz1+JZ,P[1],M-JZ);
		YplusisCtimesX(gx3,gz3+JZ,P[0],M-JZ);
		YplusisCtimesX(gx3,gz4+JZ,P[1],M-JZ);
		YplusisCtimesX(gx3,gz5+JZ,P[1],M-JZ);

		YplusisCtimesX(gx4,gz0+JY,P[1],M-JY);
		YplusisCtimesX(gx4,gz2+JY,P[1],M-JY);
		YplusisCtimesX(gx4,gz3+JY,P[1],M-JY);
		YplusisCtimesX(gx4,gz4+JY,P[0],M-JY);
		YplusisCtimesX(gx4,gz5+JY,P[1],M-JY);

		YplusisCtimesX(gx5,gz1+JX,P[1],M-JX);
		YplusisCtimesX(gx5,gz2+JX,P[1],M-JX);
		YplusisCtimesX(gx5,gz3+JX,P[1],M-JX);
		YplusisCtimesX(gx5,gz4+JX,P[1],M-JX);
		YplusisCtimesX(gx5,gz5+JX,P[0],M-JX);

		for (int k=0; k<6; k++) Times(gs+k*M,gs+k*M,g,M);
	}
}
void LGrad3::propagateB(Real *G, Real *G1, Real* P, int s_from, int s_to,int M) {
	if (lattice_type == hexagonal) {
		Real *gs=G+M*12*s_to;
		Real *gs_1=G+M*12*s_from;

		Real *gz0=gs_1;
		Real *gz1=gs_1+M;
		Real *gz2=gs_1+2*M;
		Real *gz3=gs_1+3*M;
		Real *gz4=gs_1+4*M;
		Real *gz5=gs_1+5*M;
		Real *gz6=gs_1+6*M;
		Real *gz7=gs_1+7*M;
		Real *gz8=gs_1+8*M;
		Real *gz9=gs_1+9*M;
		Real *gz10=gs_1+10*M;
		Real *gz11=gs_1+11*M;

		Real *gx0=gs;
		Real *gx1=gs+M;
		Real *gx2=gs+2*M;
		Real *gx3=gs+3*M;
		Real *gx4=gs+4*M;
		Real *gx5=gs+5*M;
		Real *gx6=gs+6*M;
		Real *gx7=gs+7*M;
		Real *gx8=gs+8*M;
		Real *gx9=gs+9*M;
		Real *gx10=gs+10*M;
		Real *gx11=gs+11*M;
		Real *g=G1;

		Zero(gs,12*M);
		remove_bounds(gz0);remove_bounds(gz1);remove_bounds(gz2);remove_bounds(gz3);remove_bounds(gz4);remove_bounds(gz5);remove_bounds(gz6);remove_bounds(gz7);remove_bounds(gz8);remove_bounds(gz9);remove_bounds(gz10);remove_bounds(gz11);
		set_bounds_x(gz0,gz11,0,0);

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
		set_bounds_y(gz3,gz8,0,0);

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
		set_bounds_z(gz5,gz6,0,0);

		YplusisCtimesX(gx0+JZ,gz6,P[1],M-JZ);
		YplusisCtimesX(gx1+JZ,gz6,P[1],M-JZ);
		YplusisCtimesX(gx2+JZ,gz6,P[1],M-JZ);
		YplusisCtimesX(gx6+JZ,gz6,P[0],M-JZ);
		YplusisCtimesX(gx7+JZ,gz6,P[0],M-JZ);
		YplusisCtimesX(gx8+JZ,gz6,P[0],M-JZ);
		YplusisCtimesX(gx9+JZ,gz6,P[1],M-JZ);
		YplusisCtimesX(gx10+JZ,gz6,P[1],M-JZ);
		YplusisCtimesX(gx11+JZ,gz6,P[1],M-JZ);

		YplusisCtimesX(gx0,gz5+JZ,P[1],M-JZ);
		YplusisCtimesX(gx1,gz5+JZ,P[1],M-JZ);
		YplusisCtimesX(gx2,gz5+JZ,P[1],M-JZ);
		YplusisCtimesX(gx3,gz5+JZ,P[0],M-JZ);
		YplusisCtimesX(gx4,gz5+JZ,P[0],M-JZ);
		YplusisCtimesX(gx5,gz5+JZ,P[0],M-JZ);
		YplusisCtimesX(gx9,gz5+JZ,P[1],M-JZ);
		YplusisCtimesX(gx10,gz5+JZ,P[1],M-JZ);
		YplusisCtimesX(gx11,gz5+JZ,P[1],M-JZ);

		remove_bounds(gz5); remove_bounds(gz6);
		set_bounds_x(gz1,gz10,0,-1);

		YplusisCtimesX(gx3+JX,gz10+JZ,P[1],M-JX-JZ);
		YplusisCtimesX(gx4+JX,gz10+JZ,P[1],M-JX-JZ);
		YplusisCtimesX(gx5+JX,gz10+JZ,P[1],M-JX-JZ);
		YplusisCtimesX(gx6+JX,gz10+JZ,P[1],M-JX-JZ);
		YplusisCtimesX(gx7+JX,gz10+JZ,P[1],M-JX-JZ);
		YplusisCtimesX(gx8+JX,gz10+JZ,P[1],M-JX-JZ);
		YplusisCtimesX(gx9+JX,gz10+JZ,P[0],M-JX-JZ);
		YplusisCtimesX(gx10+JX,gz10+JZ,P[0],M-JX-JZ);
		YplusisCtimesX(gx11+JX,gz10+JZ,P[0],M-JX-JZ);

		YplusisCtimesX(gx0+JZ,gz1+JX,P[0],M-JZ-JX);
		YplusisCtimesX(gx1+JZ,gz1+JX,P[0],M-JZ-JX);
		YplusisCtimesX(gx2+JZ,gz1+JX,P[0],M-JZ-JX);
		YplusisCtimesX(gx3+JZ,gz1+JX,P[1],M-JZ-JX);
		YplusisCtimesX(gx4+JZ,gz1+JX,P[1],M-JZ-JX);
		YplusisCtimesX(gx5+JZ,gz1+JX,P[1],M-JZ-JX);
		YplusisCtimesX(gx6+JZ,gz1+JX,P[1],M-JZ-JX);
		YplusisCtimesX(gx7+JZ,gz1+JX,P[1],M-JZ-JX);
		YplusisCtimesX(gx8+JZ,gz1+JX,P[1],M-JZ-JX);

		remove_bounds(gz1);remove_bounds(gz10);
		set_bounds_y(gz2,gz9,-1,0);

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
		set_bounds_z(gz4,gz7,0,-1);

		YplusisCtimesX(gx0+JY,gz7+JZ,P[1],M-JY-JZ);
		YplusisCtimesX(gx1+JY,gz7+JZ,P[1],M-JY-JZ);
		YplusisCtimesX(gx2+JY,gz7+JZ,P[1],M-JY-JZ);
		YplusisCtimesX(gx6+JY,gz7+JZ,P[0],M-JY-JZ);
		YplusisCtimesX(gx7+JY,gz7+JZ,P[0],M-JY-JZ);
		YplusisCtimesX(gx8+JY,gz7+JZ,P[0],M-JY-JZ);
		YplusisCtimesX(gx9+JY,gz7+JZ,P[1],M-JY-JZ);
		YplusisCtimesX(gx10+JY,gz7+JZ,P[1],M-JY-JZ);
		YplusisCtimesX(gx11+JY,gz7+JZ,P[1],M-JY-JZ);

		YplusisCtimesX(gx0+JZ,gz4+JY,P[1],M-JY-JZ);
		YplusisCtimesX(gx1+JZ,gz4+JY,P[1],M-JY-JZ);
		YplusisCtimesX(gx2+JZ,gz4+JY,P[1],M-JY-JZ);
		YplusisCtimesX(gx3+JZ,gz4+JY,P[0],M-JY-JZ);
		YplusisCtimesX(gx4+JZ,gz4+JY,P[0],M-JY-JZ);
		YplusisCtimesX(gx5+JZ,gz4+JY,P[0],M-JY-JZ);
		YplusisCtimesX(gx9+JZ,gz4+JY,P[1],M-JY-JZ);
		YplusisCtimesX(gx10+JZ,gz4+JY,P[1],M-JY-JZ);
		YplusisCtimesX(gx11+JZ,gz4+JY,P[1],M-JY-JZ);

		for (int k=0; k<12; k++) Times(gs+k*M,gs+k*M,g,M);

	} else {
		Real *gs=G+M*6*s_to;
		Real *gs_1=G+M*6*s_from;
		Real *gz0=gs_1;
		Real *gz1=gs_1+M;
		Real *gz2=gs_1+2*M;
		Real *gz3=gs_1+3*M;
		Real *gz4=gs_1+4*M;
		Real *gz5=gs_1+5*M;

		Real *gx0=gs;
		Real *gx1=gs+M;
		Real *gx2=gs+2*M;
		Real *gx3=gs+3*M;
		Real *gx4=gs+4*M;
		Real *gx5=gs+5*M;
		Real *g=G1;

		Zero(gs,6*M);
		set_bounds_x(gz0,gz5,0,0);

		YplusisCtimesX(gx1+JX,gz5,P[1],M-JX);
		YplusisCtimesX(gx2+JX,gz5,P[1],M-JX);
		YplusisCtimesX(gx3+JX,gz5,P[1],M-JX);
		YplusisCtimesX(gx4+JX,gz5,P[1],M-JX);
		YplusisCtimesX(gx5+JX,gz5,P[0],M-JX);

		set_bounds_y(gz1,gz4,0,0);

		YplusisCtimesX(gx0+JY,gz4,P[1],M-JY);
		YplusisCtimesX(gx2+JY,gz4,P[1],M-JY);
		YplusisCtimesX(gx3+JY,gz4,P[1],M-JY);
		YplusisCtimesX(gx4+JY,gz4,P[0],M-JY);
		YplusisCtimesX(gx5+JY,gz4,P[1],M-JY);

		set_bounds_z(gz2,gz3,0,0);

		YplusisCtimesX(gx0+JZ,gz3,P[1],M-JZ);
		YplusisCtimesX(gx1+JZ,gz3,P[1],M-JZ);
		YplusisCtimesX(gx3+JZ,gz3,P[0],M-JZ);
		YplusisCtimesX(gx4+JZ,gz3,P[1],M-JZ);
		YplusisCtimesX(gx5+JZ,gz3,P[1],M-JZ);

		YplusisCtimesX(gx0,gz2+JZ,P[1],M-JZ);
		YplusisCtimesX(gx1,gz2+JZ,P[1],M-JZ);
		YplusisCtimesX(gx2,gz2+JZ,P[0],M-JZ);
		YplusisCtimesX(gx4,gz2+JZ,P[1],M-JZ);
		YplusisCtimesX(gx5,gz2+JZ,P[1],M-JZ);

		YplusisCtimesX(gx0,gz1+JY,P[1],M-JY);
		YplusisCtimesX(gx1,gz1+JY,P[0],M-JY);
		YplusisCtimesX(gx2,gz1+JY,P[1],M-JY);
		YplusisCtimesX(gx3,gz1+JY,P[1],M-JY);
		YplusisCtimesX(gx5,gz1+JY,P[1],M-JY);

		YplusisCtimesX(gx0,gz0+JX,P[0],M-JX);
		YplusisCtimesX(gx1,gz0+JX,P[1],M-JX);
		YplusisCtimesX(gx2,gz0+JX,P[1],M-JX);
		YplusisCtimesX(gx3,gz0+JX,P[1],M-JX);
		YplusisCtimesX(gx4,gz0+JX,P[1],M-JX);

		for (int k=0; k<6; k++) Times(gs+k*M,gs+k*M,g,M);
	}
}

void LGrad3::propagate(Real *G, Real *G1, int s_from, int s_to,int M) { //this procedure should function on simple cubic lattice.
if (debug) cout <<" propagate in LGrad3 " << endl;
	Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	int JX_=JX, JY_=JY;
	int k=sub_box_on;

	Zero(gs,M);
	set_bounds(gs_1);
	if (k>0) {
		JX_=jx[k];
		JY_=jy[k];
	}

	if (stencil_full) {
		Add(gs+JX_,gs_1,M-JX_);
		Add(gs,gs_1+JX_,M-JX_);
		Add(gs+JY_,gs_1,M-JY_);
		Add(gs,gs_1+JY_,M-JY_);
		Add(gs+1,gs_1,M-1);
		Add(gs,gs_1+1, M-1);
		if (lattice_type == simple_cubic) {
			Norm(gs,4.0,M);
		} else {
			Norm(gs,2.0,M);
		}
		Add(gs+JX_+JY_,gs_1,M-JX_-JY_);
		Add(gs,gs_1+JX_+JY_,M-JX_-JY_);
		Add(gs+JY_,gs_1+JX,M-JY_-JX_);
		Add(gs+JX,gs_1+JY_,M-JY_-JX_);
		Add(gs+JX_+1,gs_1,M-JX_-1);
		Add(gs,gs_1+JX_+1,M-JX_-1);
		Add(gs+JX_,gs_1+1,M-JX_);
		Add(gs+1,gs_1+JX_,M-JX_);
		Add(gs+JY_+1,gs_1,M-JY_-1);
		Add(gs,gs_1+JY_+1,M-JX_-1);
		Add(gs+JY_,gs_1+1,M-JY_);
		Add(gs+1,gs_1+JY_,M-JY_);
		if (lattice_type == simple_cubic) {
			Norm(gs,4.0,M);
		} else {
			Norm(gs,2.0,M);
		}
		Add(gs+JX_+JY_+1,gs_1,M-JX_-JY_-1);
		Add(gs,gs_1+JX_+JY_+1,M-JX_-JY_-1);
		Add(gs+JX_+JY_,gs_1+1,M-JX_-JY_-1);
		Add(gs+1,gs_1+JX_+JY_,M-JX_-JY_-1);
		Add(gs+JX_+1,gs_1+JY_,M-JX_-JY_-1);
		Add(gs+JY_,gs_1+JX_+1,M-JX_-JY_-1);
		Add(gs+JY_+1,gs_1+JX_,M-JX_-JY_-1);
		Add(gs+JX_,gs_1+JY_+1,M-JX_-JY_-1);
		if (lattice_type == simple_cubic) {
			Norm(gs,1.0/152.0,M);
		} else {
			Norm(gs,1.0/56.0,M);
		}
		Times(gs,gs,G1,M);
	} else {
		if (lattice_type==simple_cubic) {
#ifdef CUDA
			Propagate_gs_locality(gs, gs_1, G1, JX, JY, JZ, M);
#else
			Add(gs+JX_,gs_1    ,M-JX_);
			Add(gs    ,gs_1+JX_,M-JX_);
			Add(gs+JY_,gs_1    ,M-JY_);
			Add(gs    ,gs_1+JY_,M-JY_);
			Add(gs+1  ,gs_1    ,M-1);
			Add(gs    ,gs_1+1  ,M-1);
			Norm(gs,1.0/6.0,M);
			Times(gs,gs,G1,M);
#endif
		} else { //hexagonal
			if (fjc==1) {
				Add(gs+JX_,gs_1,    M-JX_);
				Add(gs,    gs_1+JX_,M-JX_);
				Add(gs+JY_,gs_1    ,M-JY_);
				Add(gs,    gs_1+JY_,M-JY_);
				Add(gs+JZ, gs_1    ,M-JZ);
				Add(gs,    gs_1+JZ ,M-JZ);

				remove_bounds(gs_1);
				//set_bounds_x(gs_1,-1,0);
				Add(gs+JX_,gs_1+JY_,M-JX_-JY_);
				Add(gs+JY_,gs_1+JX_,M-JX_-JY_);

				remove_bounds(gs_1);
				//set_bounds_x(gs_1,0,-1);
				Add(gs+JX_,gs_1+JZ ,M-JX_-JZ);
				Add(gs+JZ, gs_1+JX_,M-JX_-JZ);

				remove_bounds(gs_1);
				//set_bounds_y(gs_1,0,-1);
				Add(gs+JY_,gs_1+JZ ,M-JY_-JZ);
				Add(gs+JZ, gs_1+JY_,M-JY_-JZ);

				Norm(gs,1.0/12.0,M);
				Times(gs,gs,G1,M);
			} else { //hexagonal and fjc=2
				Add(gs,            gs_1,             M);   //1
				Add(gs+JX,         gs_1,             M-JX);//2
				Add(gs,            gs_1+JX,          M-JX);//3
				Add(gs+JY,         gs_1,             M-JY);//4
				Add(gs,            gs_1+JY,          M-JY);//5
				Add(gs+1,          gs_1,             M-1);//6
				Add(gs,            gs_1+1,           M-1);//7
				Add(gs+JX+1,       gs_1,             M-JX-1);//8
				Add(gs,            gs_1+JX+1,        M-JX-1);//9
				Add(gs+1,          gs_1+JX,          M-JX-1);//10
				Add(gs+JX,         gs_1+1,           M-JX-1);//11
				Add(gs+JY+1,       gs_1,             M-JY-1);//12
				Add(gs,            gs_1+JY+1,        M-JY-1);//13
				Add(gs+1,          gs_1+JY,          M-JY-1);//14
				Add(gs+JY,         gs_1+1,           M-JY-1);//15
				Add(gs+JX+JY,      gs_1,             M-JX-JY);//16
				Add(gs,            gs_1+JX+JY,       M-JX-JY);//17
				Add(gs+JY,         gs_1+JX,          M-JX-JY);//18
				Add(gs+JX,         gs_1+JY,          M-JX-JY);//19
				Add(gs+JX+JY+1,    gs_1,             M-JX-JY-1);//20
				Add(gs,            gs_1+JX+JY+1,     M-JX-JY-1);//21
				Add(gs+JX+JY,      gs_1+1,           M-JX-JY-1);//22
				Add(gs+1,          gs_1+JX+JY,       M-JX-JY-1);//23
				Add(gs+JX+1,       gs_1+JY,          M-JX-JY-1);//24
				Add(gs+JY,         gs_1+JX+1,        M-JX-JY-1);//25
				Add(gs+JY+1,       gs_1+JX,          M-JX-JY-1);//26
				Add(gs+JX,         gs_1+JY+1,        M-JX-JY-1);//27
				Norm(gs,2.0,M);

				Add(gs+2*JX,       gs_1,             M-2*JX);//28
				Add(gs,            gs_1+2*JX,        M-2*JX);//29
				Add(gs+2*JY,       gs_1,             M-2*JY);//30
				Add(gs,            gs_1+2*JY,        M-2*JY);//31
				Add(gs+2,          gs_1,             M-2);//32
				Add(gs,            gs_1+2,           M-2);//33
				Add(gs+2*JX+1,     gs_1,     	     M-2*JX-1);//34
				Add(gs+2*JX,       gs_1+1,           M-2*JX-1);//35
				Add(gs+1,          gs_1+2*JX,        M-2*JX-1);//36
				Add(gs,            gs_1+2*JX+1,      M-2*JX-1);//37
				Add(gs+2*JY,       gs_1+1,           M-2*JY-1);//38
				Add(gs+2*JY+1,     gs_1,             M-2*JY-1);//39
				Add(gs+1,          gs_1+2*JY,        M-2*JY-1);//40
				Add(gs,            gs_1+2*JY+1,      M-2*JY-1);//41
				Add(gs+JX+2,       gs_1,             M-JX-2);//42
				Add(gs+JX,         gs_1+2,           M-JX-2);//43
				Add(gs,            gs_1+JX+2,        M-JX-2);//44
				Add(gs+2,          gs_1+JX,          M-JX-2);//45
				Add(gs+JY+2,       gs_1,             M-JY-2);//46
				Add(gs+JY,         gs_1+2,           M-JY-2);//47
				Add(gs,            gs_1+JY+2,        M-JY-2);//48
				Add(gs+2,          gs_1+JY,          M-JY-2);//49
				Add(gs+2*JX,       gs_1+JY,          M-2*JX-JY);//50
				Add(gs+2*JX+JY,    gs_1,             M-2*JX-JY);//51
				Add(gs+JY,         gs_1+2*JX,        M-2*JX-JY);//52
				Add(gs,            gs_1+2*JX+JY,     M-2*JX-JY);//53
				Add(gs+2*JY,       gs_1+JX,          M-2*JY-JX);//54
				Add(gs+2*JY+JX,    gs_1,             M-2*JY-JX);//55
				Add(gs+JX,         gs_1+2*JY,        M-2*JY-JX);//56
				Add(gs,            gs_1+2*JY+JX,     M-2*JY-JX);//57
				Add(gs+2*JX+JY+1,  gs_1,             M-2*JX-JY-1);//58
				Add(gs+2*JX+1,     gs_1+JY,          M-2*JX-JY-1);//59
				Add(gs+JY+1,       gs_1+2*JX,        M-2*JX-JY-1);//60
				Add(gs+1,          gs_1+2*JX+JY,     M-2*JX-JY-1);//61
				Add(gs+2*JX+JY,    gs_1+1,    	     M-2*JX-JY-1);//62
				Add(gs+2*JX,       gs_1+JY+1,        M-2*JX-JY-1);//63
				Add(gs+JY,         gs_1+2*JX+1,      M-2*JX-JY-1);//64
				Add(gs,            gs_1+2*JX+JY+1,   M-2*JX-JY-1);//65
				Add(gs+JX+2*JY+1,  gs_1,             M-JX-2*JY-1);//66
				Add(gs+JX+1,       gs_1+2*JY,        M-JX-2*JY-1);//67
				Add(gs+2*JY+1,     gs_1+JX,          M-JX-2*JY-1);//68
				Add(gs+1,          gs_1+JX+2*JY,     M-JX-2*JY-1);//69
				Add(gs+JX+2*JY,    gs_1+1,           M-JX-2*JY-1);//70
				Add(gs+JX,         gs_1+2*JY+1,      M-JX-2*JY-1);//71
				Add(gs+2*JY,       gs_1+JX+1,        M-JX-2*JY-1);//72
				Add(gs,            gs_1+JX+2*JY+1,   M-JX-2*JY-1);//73
				Add(gs+JX+JY+2,    gs_1,             M-JX-JY-2);//74
				Add(gs+JX+2,       gs_1+JY,          M-JX-JY-2);//75
				Add(gs+JY+2,       gs_1+JX,          M-JX-JY-2);//76
				Add(gs+2,          gs_1+JX+JY,       M-JX-JY-2);//77
				Add(gs+JX+JY,      gs_1+2,           M-JX-JY-2);//78
				Add(gs+JX,         gs_1+JY+2,        M-JX-JY-2);//79
				Add(gs+JY,         gs_1+JX+2,        M-JX-JY-2);//80
				Add(gs,            gs_1+JX+JY+2,     M-JX-JY-2);//81
				Norm(gs,2.0,M);

				Add(gs+2*JX+2*JY,  gs_1,             M-2*JX-2*JY);//82
				Add(gs,            gs_1+2*JX+2*JY,   M-2*JX-2*JY);//83
				Add(gs+2*JX,       gs_1+2*JY,        M-2*JX-2*JY);//84
				Add(gs+2*JY,       gs_1+2*JX,        M-2*JX-2*JY);//85
				Add(gs+2*JX+2*JY+1,gs_1,             M-2*JX-2*JY-1);//86
				Add(gs+1,          gs_1+2*JX+2*JY,   M-2*JX-2*JY-1);//87
				Add(gs+2*JX+1,     gs_1+2*JY,        M-2*JX-2*JY-1);//88
				Add(gs+2*JY+1,     gs_1+2*JX,        M-2*JX-2*JY-1);//89
				Add(gs+2*JX+2*JY,  gs_1+1,           M-2*JX-2*JY-1);//90
				Add(gs,            gs_1+2*JX+2*JY+1, M-2*JX-2*JY-1);//91
				Add(gs+2*JX,       gs_1+2*JY+1,      M-2*JX-2*JY-1);//92
				Add(gs+2*JY,       gs_1+2*JX+1,      M-2*JX-2*JY-1);//93
				Add(gs+2*JX+2,     gs_1,             M-2*JX-2);//94
				Add(gs+2*JX+JY+2,  gs_1,             M-2*JX-JY-2);//95
				Add(gs+2*JX+2,     gs_1+JY,          M-2*JX-JY-2);//96
				Add(gs+2,          gs_1+2*JX,        M-2*JX-2);//97
				Add(gs+JY+2,       gs_1+2*JX,        M-2*JX-JY-2);//98
				Add(gs+2,          gs_1+2*JX+JY,     M-2*JX-JY-2);//99
				Add(gs+2*JY+2,     gs_1,             M-2*JY-2);//100
				Add(gs+2*JY+JX+2,  gs_1,             M-2*JY-JX-2);//101
				Add(gs+2*JY+2,     gs_1+JX,          M-2*JY-JX-2);//102
				Add(gs+2,          gs_1+2*JY,        M-2*JY-2);//103
				Add(gs+JX+2,       gs_1+2*JY,        M-2*JY-JX-2);//104
				Add(gs+2,          gs_1+2*JY+JX,     M-2*JY-JX-2);//105
				Add(gs+2*JX,       gs_1+2,           M-2*JX-2);//106
				Add(gs+2*JX+JY,    gs_1+2,           M-2*JX-JY-2);//107
				Add(gs+2*JX,       gs_1+JY+2,        M-2*JX-JY-2);//108
				Add(gs,            gs_1+2*JX+2,      M-2*JX-2);//109
				Add(gs+JY,         gs_1+2*JX+2,      M-2*JX-JY-2);//110
				Add(gs,            gs_1+2*JX+JY+2,   M-2*JX-JY-2);//111
				Add(gs+2*JY,       gs_1+2,           M-2*JY-2);//112
				Add(gs+2*JY+JX,    gs_1+2,           M-2*JY-JX-2);//113
				Add(gs+2*JY,       gs_1+JX+2,        M-2*JY-JX-2);//114
				Add(gs,            gs_1+2*JY+2,      M-2*JY-2);//115
				Add(gs+JX,         gs_1+2*JY+2,      M-2*JY-JX-2);//116
				Add(gs,            gs_1+2*JY+JX+2,   M-2*JY-JX-2);//117
				Norm(gs,2.0,M);

				Add(gs+2*JX+2*JY+2,gs_1,             M-2*JX-2*JY-2);//118
				Add(gs+2,          gs_1+2*JX+2*JY,   M-2*JX-2*JY-2);//119
				Add(gs+2*JX+2,     gs_1+2*JY,        M-2*JX-2*JY-2);//120
				Add(gs+2*JY+2,     gs_1+2*JX,        M-2*JX-2*JY-2);//121
				Add(gs+2*JX+2*JY,  gs_1+2,           M-2*JX-2*JY-2);//122
				Add(gs,            gs_1+2*JX+2*JY+2, M-2*JX-2*JY-2);//123
				Add(gs+2*JX,       gs_1+2*JY+2,      M-2*JX-2*JY-2);//124
				Add(gs+2*JY,       gs_1+2*JX+2,      M-2*JX-2*JY-2);//125
				Norm(gs,1.0/512.0,M);

				Times(gs,gs,G1,M);
			}
		}
	}
}


bool LGrad3::ReadRange(int* r, int* H_p, int &n_pos, bool &block, string range, int var_pos, string seg_name, string range_type) {
if (debug) cout <<"ReadRange in LGrad3 " << endl;
	bool success=true;
	vector<string>set;
	vector<string>coor;
	vector<string>xyz;
	In[0]->split(range,';',set);
	coor.clear();
	block=true; In[0]->split(set[0],',',coor);

	if (coor.size()!=3) {cout << "In mon " + 	seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes do not come in set of three: 'x,y,z'" << endl; success=false;}
	else {
		r[0]=In[0]->Get_int(coor[0],0);
		r[1]=In[0]->Get_int(coor[1],0);
		r[2]=In[0]->Get_int(coor[2],0);
	}
	coor.clear(); In[0]->split(set[1],',',coor);

	if (coor.size()!=3) {cout << "In mon " + seg_name+ ", for 'pos 2', in '" + range_type + "', the coordinates do not come in set of three: 'x,y,z'" << endl; success=false;}
	else {
		r[3]=In[0]->Get_int(coor[0],0);
		r[4]=In[0]->Get_int(coor[1],0);
		r[5]=In[0]->Get_int(coor[2],0);
	}
	if (r[0] > r[3]) {cout << "In mon " + seg_name+ ", for 'pos 1', the x-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
	if (r[1] > r[4]) {cout << "In mon " + seg_name+ ", for 'pos 1', the y-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
	if (r[2] > r[5]) {cout << "In mon " + seg_name+ ", for 'pos 1', the z-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
	return success;
}

bool LGrad3::ReadRangeFile(string filename,int* H_p, int &n_pos, string seg_name, string range_type) {
if (debug) cout <<"ReadRangeFile in LGrad3 " << endl;
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
	int px,py,pz,p_i,x,y,z;
	int i=0;
	if (!In[0]->ReadFile(sub[0].append(".").append(filename),content)) {
		success=false;
		return success;
	}

	In[0]->split(content,'#',lines);
	length = lines.size();
	if (length == MX*MY*MZ) { //expect to read 'mask file';
		i=0;
		if (n_pos==0) {
			while (i<length){
				if (In[0]->Get_int(lines[i],0)==1) n_pos++;
				i++;
			};
			if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any unities." << endl;}
		} else {
			i=0; p_i=0;
			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++) for (z=1; z<MZ+1; z++) {
				if (In[0]->Get_int(lines[i],0)==1) {H_p[p_i]=x*JX+y*JY+fjc-1+z; p_i++;}
				i++;
			}
		}
	} else { //expect to read x,y,z
		px=0,py=0,pz=0; i=0;
		if (n_pos==0) n_pos=length;
		else {
			while (i<length) {
				xyz.clear();
				In[0]->split(lines[i],',',xyz);
				length_xyz=xyz.size();
				if (length_xyz!=3) {
					cout << "In mon " + seg_name + " " +range_type+"_filename  the expected 'triple coordinate' structure 'x,y,z' was not found. " << endl;  success = false;
				} else {
					px=In[0]->Get_int(xyz[0],0);
					if (px < 1 || px > MX) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the x-coordinate in "+range_type+"_filename out of bounds: 1.." << MX << endl; success =false;}
					py=In[0]->Get_int(xyz[1],0);
					if (py < 1 || py > MY) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the y-coordinate in "+range_type+"_filename out of bounds: 1.." << MY << endl; success =false;}
					pz=In[0]->Get_int(xyz[2],0);
					if (pz < 1 || pz > MZ) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the y-coordinate in "+range_type+"_filename out of bounds: 1.." << MZ << endl; success =false;}
				}
				H_p[i]=px*JX+py*JY+fjc-1+pz;
				i++;
			}
		}
	}
	return success;
}

bool LGrad3::FillMask(int* Mask, vector<int>px, vector<int>py, vector<int>pz, string filename) {
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
		if (MX*MY*MZ!=length) {
			success=false; cout <<"inputfile for filling delta_range has not expected length in x,y,z-directions" << endl;
		} else {
			for (int x=1; x<MX+1; x++)
			for (int y=1; y<MY+1; y++)
			for (int z=1; z<MZ+1; z++) Mask[x*JX + y*JY + z]=In[0]->Get_int(lines[x*JX + y*JY + z],-1);
		}
	} else  {
		for (int i=0; i<length_px; i++) {
			p=px[i]; if (p<1 || p>MX) {success=false; cout<<" x-value in delta_range out of bounds; " << endl; }
			p=py[i]; if (p<1 || p>MY) {success=false; cout<<" y-value in delta_range out of bounds; " << endl; }
			p=pz[i]; if (p<1 || p>MZ) {success=false; cout<<" z-value in delta_range out of bounds; " << endl; }
			if (success) Mask[px[i]*JX + py[i]*JY + fjc-1+ pz[i]]=1;
		}
	}
	for (int i=0; i<M; i++) if (!(Mask[i]==0 || Mask[i]==1)) {success =false; cout <<"Delta_range does not contain '0' or '1' values. Check delta_inputfile values"<<endl; }
	return success;
}

bool LGrad3::CreateMASK(int* H_MASK, int* r, int* H_P, int n_pos, bool block) {
if (debug) cout <<"CreateMask for LGrad3 " + name << endl;
	bool success=true;
	H_Zero(H_MASK,M);
	if (block) {
		for (int x=r[0]; x<r[3]+1; x++)
		for (int y=r[1]; y<r[4]+1; y++)
		for (int z=r[2]; z<r[5]+1; z++)
			H_MASK[P(x,y,z)]=1;
	} else {
		for (int i = 0; i<n_pos; i++) H_MASK[H_P[i]]=1;
	}
	return success;
}


Real LGrad3::ComputeTheta(Real* phi) {
	Real result=0; remove_bounds(phi);
	Dot(result,phi,L,M);
	return result;
}

void LGrad3::UpdateEE(Real* EE, Real* psi, Real* E) {
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_M_bounds(psi);

	Zero(EE,M);
	AddGradSquare(EE+1,psi,psi+1,psi+2,M-2);
	AddGradSquare(EE+JX,psi,psi+JX,psi+2*JX,M-2*JX);
	AddGradSquare(EE+JY,psi,psi+JY,psi+2*JY,M-2*JY);
	Norm(EE,pf,M);

/* In lattice refinement. this is possibly not working...old version restored.
			YisAminB(EE+1,psi,psi+2,M-2);
			Times(EE,EE,EE,M);
			Norm(EE,pf/4.0,M);
			YisAminB(E+JX,psi,psi+2*JX,M-2*JX);
			Times(E,E,E,M);
			YplusisCtimesX(EE,E,pf/4.0,M);
			YisAminB(E+JY,psi,psi+2*JY,M-2*JY);
			Times(E,E,E,M);
			YplusisCtimesX(EE,E,pf/4.0,M);
*/
}


void LGrad3::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, int* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	int x,y;
#ifndef CUDA
	int z;
#endif

#ifndef CUDA
	Real epsZplus, epsZmin;
#endif
	Real epsXplus, epsXmin, epsYplus,epsYmin;
	//set_M_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);

   if (!fixedPsi0) {

#ifdef CUDA
	C *=6;
	Cp(X,psi,M);
	UpPsi(g+JX+JY+1,psi+JX+JY+1,X+JX+JY+1,eps+JX+JY+1,JX,JY,C,Mask+JX+JY+1,M-2*(JX+JY+1));
	if (fjc==2) cout << "in GPU FJC-choices > 3 not implemented yet " << endl;
#else

	for (x=1; x<MX+1; x++) {
		for (y=1; y<MY+1; y++) {
			epsZplus=eps[x*JX+y*JY]+eps[x*JX+y*JY+1];
			for (z=1; z<MZ+1; z++) {
				epsZmin=epsZplus;
				epsZplus=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z+1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
					epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
				}
			}
		}
	}
	for (x=MX; x>0; x--) {
		for (y=MY; y>0; y--) {
			epsZmin=eps[x*JX+y*JY+MZ+1]+eps[x*JX+y*JY+MZ];
			for (z=MZ; z>0; z--) {
				epsZplus=epsZmin;
				epsZmin=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z-1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
					epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
					g[x*JX+y*JY+z]-=psi[x*JX+y*JY+z];
				}
			}
		}
	}



/*
	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++) {
		X[x*JX+y*JY+z]=(psi[(x-1)*JX+y*JY+z]+psi[(x+1)*JX+y*JY+z]
				 +psi[x*JX+(y-1)*JY+z]+psi[x*JX+(y+1)*JY+z]
			        +psi[x*JX+y*JY+z-1]  +psi[x*JX+y*JY+z+1])/6.0
			        +q[x*JX+y*JY+z]*C/eps[x*JX+y*JY+z]/3.0;
	}

	if (grad_epsilon) for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)  for (z=fjc; z<MZ+fjc; z++){
		X[x*JX+y*JY+z]+=0.25*(eps[(x+1)*JX+y*JY+z]-eps[(x-1)*JX+y*JY+z])*(psi[(x+1)*JX+y*JY+z]-psi[(x-1)*JX+y*JY+z]+
                              eps[x*JX+(y+1)*JY+z]-eps[x*JX+(y-1)*JY+z])*(psi[x*JX+(y+1)*JY+z]-psi[x*JX+(y-1)*JY+z]+
	                       eps[x*JX+y*JY+z-1]  -eps[x*JX+y*JY+z+1])  *(psi[x*JX+y*JY+z-1]  -psi[x*JX+y*JY+z+1])
				  /eps[x*JX+y*JY+z];
	}
	Cp(psi,X,M);
	YisAminB(g,g,psi,M);
*/

#endif
   } else { //fixedPsi0 is true

#ifdef CUDA
	C *=6;
	Cp(X,psi,M);
	UpPsi(g+JX+JY+1,psi+JX+JY+1,X+JX+JY+1,eps+JX+JY+1,JX,JY,C,Mask+JX+JY+1,M-2*(JX+JY+1));
#else
	for (x=1; x<MX+1; x++) {
		for (y=1; y<MY+1; y++) {
			epsZplus=eps[x*JX+y*JY]+eps[x*JX+y*JY+1];
			for (z=1; z<MZ+1; z++) {
				epsZmin=epsZplus;
				epsZplus=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z+1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0)
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					    epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
				           epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					    C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
			}
		}
	}
	for (x=MX; x>0; x--) {
		for (y=MY; y>0; y--) {
			epsZmin=eps[x*JX+y*JY+MZ+1]+eps[x*JX+y*JY+MZ];
			for (z=MZ; z>0; z--) {
				epsZplus=epsZmin;
				epsZmin=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z-1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
					epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
					g[x*JX+y*JY+z]-=psi[x*JX+y*JY+z];
				}
			}
		}
	}


/*  in lattice refinement this part is appartently not working. previous lines is for previous version...

	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++) {
		if (Mask[x*JX+y*JY+z] == 0)
		X[x*JX+y*JY+z]=(psi[(x-1)*JX+y*JY+z]+psi[(x+1)*JX+y*JY+z]
				 +psi[x*JX+(y-1)*JY+z]+psi[x*JX+(y+1)*JY+z]
			        +psi[x*JX+y*JY+z-1]  +psi[x*JX+y*JY+z+1])/6.0
			          +q[x*JX+y*JY+z]*C/eps[x*JX+y*JY+z]/3.0;
	}

	if (grad_epsilon) for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)  for (z=fjc; z<MZ+fjc; z++){
		if (Mask[x*JX+y*JY+z] == 0)
		X[x*JX+y*JY+z]+=0.25*(eps[(x+1)*JX+y*JY+z]-eps[(x-1)*JX+y*JY+z])*(psi[(x+1)*JX+y*JY+z]-psi[(x-1)*JX+y*JY+z]+
                              eps[x*JX+(y+1)*JY+z]-eps[x*JX+(y-1)*JY+z])*(psi[x*JX+(y+1)*JY+z]-psi[x*JX+(y-1)*JY+z]+
			         eps[x*JX+y*JY+z-1]  -eps[x*JX+y*JY+z+1])  *(psi[x*JX+y*JY+z-1]  -psi[x*JX+y*JY+z+1])
				   /eps[x*JX+y*JY+z];
	}
	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)  for (z=fjc; z<MZ+fjc; z++)
	if (Mask[x*JX+y*JY+z] == 0) {
		psi[x*JX+y*JY+z]=X[x*JX+y*JY+z];
		g[x*JX+y*JY+z]-=psi[x*JX+y*JY+z];
	}
*/

#endif
   }
}


void LGrad3::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, int* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
	int x,y;
	#ifndef CUDA
	int z;
	#endif
	Real epsXplus,epsXmin,epsYplus,epsYmin,epsZplus,epsZmin;

	Real C = -e*e/(eps0*k_BT*bond_length);
#ifdef CUDA
	C *=6;
	Cp(X,psi,M);
	UpQ(g+JX+JY+1, q+JX+JY+1, psi+JX+JY+1, eps+JX+JY+1, JX, JY, C, Mask+JX+JY+1, M-2*(JX+JY+1));
#else

	for (x=1; x<MX; x++) {
		for (y=1; y<MY; y++) {
			epsZplus=eps[x*JX+y*JY]+eps[x*JX+y*JY+1];
			for (z=1; z<MZ; z++) {
				epsZmin=epsZplus;
				epsZplus=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z+1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==1) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					                   epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
							     epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]-
							     (epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus)*psi[x*JX+y*JY+z])/C;
					g[x*JX+y*JY+z]=-q[x*JX+y*JY+z];
				}
			}
		}
	}


/* in lattice refinement charge in 3d not working trying to restore....

	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++) {
		if (Mask[x*JX+y*JY+z] == 1)
		q[x*JX+y*JY+z]=-0.5*(psi[(x-1)*JX+y*JY+z]+psi[(x+1)*JX+y*JY+z]
				 +psi[x*JX+(y-1)*JY+z]+psi[x*JX+(y+1)*JY+z]
			        +psi[x*JX+y*JY+z-1]  +psi[x*JX+y*JY+z+1]
			       -6.0*q[x*JX+y*JY+z])*fjc*fjc*eps[x*JX+y*JY+z]/C;
	}

	if (grad_epsilon) for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)  for (z=fjc; z<MZ+fjc; z++){
		if (Mask[x*JX+y*JY+z] == 1)
		q[x*JX+y*JY+z]-=0.25*(eps[(x+1)*JX+y*JY+z]-eps[(x-1)*JX+y*JY+z])*(psi[(x+1)*JX+y*JY+z]-psi[(x-1)*JX+y*JY+z]+
                              eps[x*JX+(y+1)*JY+z]-eps[x*JX+(y-1)*JY+z])*(psi[x*JX+(y+1)*JY+z]-psi[x*JX+(y-1)*JY+z]+
				  eps[x*JX+y*JY+z-1]  -eps[x*JX+y*JY+z+1])  *(psi[x*JX+y*JY+z-1]  -psi[x*JX+y*JY+z+1])*fjc*fjc/C;
	}
	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)  for (z=fjc; z<MZ+fjc; z++)
	if (Mask[x*JX+y*JY+z] == 1) {
		g[x*JX+y*JY+z]=-q[x*JX+y*JY+z];
	}
*/

#endif

}

void LGrad3::set_bounds_x(Real* X,Real* Y,int shifty,int shiftz){
if (debug) cout <<"set_bounds_x (shift in y,z ) in LGrad3 " << endl;
	int y,z;
	int k=0;
	if (BX1>BXM) {
		set_bounds_x(X,shifty,shiftz);
		set_bounds_x(Y,shifty,shiftz);
	} else {

		if (fjc==1) {
			 for (y=1; y<MY+1; y++) for (z=1; z<MZ+1; z++)  {
				X[0+        y*JY+z*JZ] = Y[BX1*JX+(y+shifty)*JY+(z+shiftz)*JZ];
				X[(MX+1)*JX+y*JY+z*JZ] = Y[BXM*JX+(y-shifty)*JY+(z-shiftz)*JZ];
				Y[0        +y*JY+z*JZ] = X[BX1*JX+(y+shifty)*JY+(z+shiftz)*JZ];
				Y[(MX+1)*JX+y*JY+z*JZ] = X[BXM*JX+(y-shifty)*JY+(z-shiftz)*JZ];
			}
		} else {
			for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
				for (k=0; k<fjc; k++) {
					X[k*JX+y*JY+z*JZ] = Y[B_X1[k]*JX+(y+shifty)*JY+(z+shiftz)*JZ];
					Y[k*JX+y*JY+z*JZ] = X[B_X1[k]*JX+(y+shifty)*JY+(z+shiftz)*JZ];

				}
				for (k=0; k<fjc; k++) {
					X[(MX+fjc+k)*JX+y*JY+z*JZ] = Y[B_XM[k]*JX+(y-shifty)*JY+(z-shiftz)*JZ];
					Y[(MX+fjc+k)*JX+y*JY+z*JZ] = X[B_XM[k]*JX+(y-shifty)*JY+(z-shiftz)*JZ];
				}
			}
		}
	}
}
void LGrad3::set_bounds_y(Real* X,Real* Y,int shiftx, int shiftz){
if (debug) cout <<"set_bounds_y (shift in x,z) in LGrad3 " << endl;
	int x,z;
	int k=0;
	if (BY1>BYM) {
		set_bounds_y(X,shiftx,shiftz);
		set_bounds_y(Y,shiftx,shiftz);
	} else {

		if (fjc==1) {
			for (z=1; z<MZ+1; z++) for (x=1; x<MX+1; x++){
				X[x*JX+0+        z*JZ] = Y[(x+shiftx)*JX+BY1*JY+(z+shiftz)*JZ];
				X[x*JX+(MY+1)*JY+z*JZ] = Y[(x-shiftx)*JX+BYM*JY+(z-shiftz)*JZ];
				Y[x*JX+0+        z*JZ] = X[(x+shiftx)*JX+BY1*JY+(z+shiftz)*JZ];
				Y[x*JX+(MY+1)*JY+z*JZ] = X[(x-shiftx)*JX+BYM*JY+(z-shiftz)*JZ];
			}
		} else {
			for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){
				for (k=0; k<fjc; k++) {
					X[x*JX+k*JY+z*JZ] = Y[(x+shiftx)*JX+B_Y1[k]*JY+(z+shiftz)*JZ];
					Y[x*JX+k*JY+z*JZ] = X[(x+shiftx)*JX+B_Y1[k]*JY+(z+shiftz)*JZ];
				}
				for (k=0; k<fjc; k++) {
					X[x*JX+(MY+fjc+k)*JY+z*JZ] = Y[(x-shiftx)*JX+B_YM[k]*JY+(z-shiftz)*JZ];
					Y[x*JX+(MY+fjc+k)*JY+z*JZ] = X[(x-shiftx)*JX+B_YM[k]*JY+(z-shiftz)*JZ];
				}
			}
		}
	}

}

void LGrad3::set_bounds_z(Real* X,Real* Y,int shiftx,int shifty){
if (debug) cout <<"set_bounds_z shift (x,y) in LGrad3 " << endl;
	int x,y;
	int k=0;
	if (BZ1>BZM) { //periodic
		set_bounds_z(X,shiftx,shifty);
		set_bounds_z(Y,shiftx,shifty);
	} else {
		if (fjc==1) {
			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++) {
				X[x*JX+y*JY+   0] = Y[(x+shiftx)*JX+(y+shifty)*JY+BZ1];
				X[x*JX+y*JY+MZ+1] = Y[(x-shiftx)*JX+(y-shifty)*JY+BZM];
				Y[x*JX+y*JY+   0] = X[(x+shiftx)*JX+(y+shifty)*JY+BZ1];
				Y[x*JX+y*JY+MZ+1] = X[(x-shiftx)*JX+(y-shifty)*JY+BZM];
			}
		} else {
			for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++) {
				for (k=0; k<fjc; k++) {
					X[x*JX+y*JY+k] = Y[(x+shiftx)*JX+(y+shifty)*JY+B_Z1[k]];
					Y[x*JX+y*JY+k] = X[(x+shiftx)*JX+(y+shifty)*JY+B_Z1[k]];
				}
				for (k=0; k<fjc; k++) {
					X[x*JX+y*JY+MZ+fjc+k]  = Y[(x-shiftx)*JX+(y-shifty)*JY+B_ZM[k]];
					Y[x*JX+y*JY+MZ+fjc+k]  = X[(x-shiftx)*JX+(y-shifty)*JY+B_ZM[k]];
				}
			}
		}
	}
}

void LGrad3::set_bounds_x(Real* X, int shifty, int shiftz){
if (debug) cout <<"set_bounds_x (shift yz) in LGrad3 " << endl;
	int y,z;
	int k=0;
	//if (BX1>BXM) {shifty=0; shiftz=0;} //periodic
	if (fjc==1) {
		 for (y=1; y<MY+1; y++) for (z=1; z<MZ+1; z++)  {
			X[0        +y*JY+z*JZ] = X[BX1*JX+(y+shifty)*JY+(z+shiftz)*JZ];
			X[(MX+1)*JX+y*JY+z*JZ] = X[BXM*JX+(y-shifty)*JY+(z-shiftz)*JZ];
		}
	} else {//shifts for hexagonal lattice not put in the fjc>1 case...
		for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
			for (k=0; k<fjc; k++) X[k*JX+y*JY+z*JZ] = X[B_X1[k]*JX+(y+shifty)*JY+(z+shiftz)*JZ];
			for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = X[B_XM[k]*JX+(y-shifty)*JY+(z-shiftz)*JZ];
		}
	}
}

void LGrad3::set_bounds_y(Real* X,int shiftx, int shiftz){
if (debug) cout <<"set_bounds_y (shift x,z) in LGrad3 " << endl;
	int x,z;
	int k=0;
	//if (BY1>BYM) {shiftx=0; shiftz=0;} //periodic
	if (fjc==1) {
		for (z=1; z<MZ+1; z++) for (x=1; x<MX+1; x++){
			X[x*JX+    0    +z*JZ] = X[(x+shiftx)*JX+BY1*JY+(z+shiftz)*JZ];
			X[x*JX+(MY+1)*JY+z*JZ] = X[(x-shiftx)*JX+BYM*JY+(z-shiftz)*JZ];
		}
	} else {
		for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){
			for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = X[(x+shiftx)*JX+B_Y1[k]*JY+(z+shiftz)*JZ];
			for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = X[(x-shiftx)*JX+B_YM[k]*JY+(z-shiftz)*JZ];
		}
	}
}

void LGrad3::set_bounds_z(Real* X,int shiftx, int shifty){
if (debug) cout <<"set_bounds_z (shift xy) in LGrad3 " << endl;
	int x,y;
	int k=0;
	//if (BZ1>BZM) {shiftx=0; shifty=0;} //periodic.
	if (fjc==1) {
		for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++) {
			X[x*JX+y*JY+   0] = X[(x+shiftx)*JX+(y+shifty)*JY+BZ1];
			X[x*JX+y*JY+MZ+1] = X[(x-shiftx)*JX+(y-shifty)*JY+BZM];
		}
	} else {
		for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++){
			for (k=0; k<fjc; k++) X[x*JX+y*JY+k] = X[(x+shiftx)*JX+(y+shifty)*JY+B_Z1[k]];
			for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[(x-shiftx)*JX+(y-shifty)*JY+B_ZM[k]];
		}
	}
}

void LGrad3::remove_bounds(Real *X){
if (debug) cout <<"remove_bounds in LGrad3 " << endl;
	int x,y,z;
	int k;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			RemoveBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) RemoveBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); else {
			for (x=0; x<MX+2*fjc; x++) for (y=0; y<MY+2*fjc; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+k] = 0;
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = 0;
			}
			for (y=0; y<MY+2*fjc; y++) for (z=0; z<MZ+2*fjc; z++)  {
				for (k=0; k<fjc; k++) X[k*JX+y*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = 0;
			}
			for (z=0; z<MZ+2*fjc; z++) for (x=0; x<MX+2*fjc; x++){
				for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = 0;
			}
		}
	}
}

void LGrad3::set_bounds(Real* X){
if (debug) cout <<"set_bounds in LGrad3 " << endl;
	int x,y,z;
	int k=0;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			SetBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) {

			//SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);

			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++){
				X[x*JX+y*JY+0]     = X[x*JX+y*JY+BZ1];
				X[x*JX+y*JY+MZ+1]  = X[x*JX+y*JY+BZM];
			}
			for (y=1; y<MY+1; y++) for (z=1; z<MZ+1; z++)  {
				X[0        +y*JY+z*JZ] = X[BX1*JX+y*JY+z*JZ];
				X[(MX+1)*JX+y*JY+z*JZ] = X[BXM*JX+y*JY+z*JZ];
			}
			for (z=1; z<MZ+1; z++) for (x=1; x<MX+1; x++){
				X[x*JX+0        +z*JZ] = X[x*JX+BY1*JY+z*JZ];
				X[x*JX+(MY+1)*JY+z*JZ] = X[x*JX+BYM*JY+z*JZ];
			}

			x=0; {
				for (y=1; y<MY+1; y++){
					X[y*JY+0]     = X[y*JY+BZ1];
					X[y*JY+MZ+1]  = X[y*JY+BZM];
				}
				for (z=1; z<MZ+1; z++){
					X[0        +z*JZ] = X[BY1*JY+z*JZ];
					X[(MY+1)*JY+z*JZ] = X[BYM*JY+z*JZ];
				}
			}
			x=MX+1; {
				for (y=1; y<MY+1; y++){
					X[x*JX+y*JY+0]     = X[x*JX+y*JY+BZ1];
					X[x*JX+y*JY+MZ+1]  = X[x*JX+y*JY+BZM];
				}
				for (z=1; z<MZ+1; z++){
					X[x*JX+0        +z*JZ] = X[x*JX+BY1*JY+z*JZ];
					X[x*JX+(MY+1)*JY+z*JZ] = X[x*JX+BYM*JY+z*JZ];
				}
			}
			y=0; {
				for (x=1; x<MX+1; x++){
					X[x*JX        +0] = X[x*JX+BZ1*JZ];
					X[x*JX+(MZ+1)*JZ] = X[x*JX+BZM*JZ];
				}
			}
			y=MY+1; {
				for (x=1; x<MX+1; x++){
					X[x*JX+y*JY        +0] = X[x*JX+y*JY+BZ1*JZ];
					X[x*JX+y*JY+(MZ+1)*JZ] = X[x*JX+y*JY+BZM*JZ];
				}
			}

			X[0        +0        +0]        =X[BX1*JX+BY1*JY+BZ1*JZ];
			X[(MX+1)*JX+0        +0]        =X[BXM*JX+BY1*JY+BZ1*JZ];
			X[0        +(MY+1)*JY+0]        =X[BX1*JX+BYM*JY+BZ1*JZ];
			X[0        +0        +(MZ+1)*JZ]=X[BX1*JX+BY1*JY+BZM*JZ];
			X[(MX+1)*JX+(MY+1)*JY+0]        =X[BXM*JX+BYM*JY+BZ1*JZ];
			X[(MX+1)*JX+ 0       +(MZ+1)*JZ]=X[BXM*JX+BY1*JY+BZM*JZ];
			X[0        +(MY+1)*JY+(MZ+1)*JZ]=X[BX1*JX+BYM*JY+BZM*JZ];
			X[(MX+1)*JX+(MY+1)*JY+(MZ+1)*JZ]=X[BXM*JX+BYM*JY+BZM*JZ];


		} else {
			for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+k] = X[x*JX+y*JY+B_Z1[k]];
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
			}
			for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
				for (k=0; k<fjc; k++) X[k*JX+y*JY+z*JZ] = X[B_X1[k]*JX+y*JY+z*JZ];
				for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = X[B_XM[k]*JX+y*JY+z*JZ];
			}
			for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){
				for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = X[x*JX+B_Y1[k]*JY+z*JZ];
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = X[x*JX+B_YM[k]*JY+z*JZ];
			}

			for (x=0; x<fjc; x++ ) {
				for (y=fjc; y<MY+fjc; y++){
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+k]         = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
				}
				for (z=fjc; z<MZ+fjc; z++){
					for (int k=0; k<fjc; k++) X[x*JX+k*JY         +z]  = X[x*JX+B_Y1[k]*JY+z*JZ];
					for (int k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z] = X[x*JX+B_YM[k]*JY+z*JZ];
				}
			}
			for (x=MX+fjc; x<MX+2*fjc; x++) {
				for (y=fjc; y<MY+fjc; y++){
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+k]         = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
				}
				for (z=fjc; z<MZ+fjc; z++){
					for (int k=0; k<fjc; k++) X[x*JX+k*JY        +z]  = X[x*JX+B_Y1[k]*JY+z];
					for (int k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z] = X[x*JX+B_YM[k]*JY+z];
				}
			}
			for (y=0; y<fjc; y++) {
				for (x=fjc; x<MX+fjc; x++){
					for (int k=0; k<fjc; k++) X[x*JX +y*JY       +k]    = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX +y*JY +(MZ+fjc+k)] = X[x*JX+y*JY+B_ZM[k]];
				}
			}
			for (y=MY+fjc; y<MY+2*fjc; y++) {
				for (x=fjc; x<MX+fjc; x++){
					for (int k=0; k<fjc; k++) X[x*JX+y*JY        +k]     = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+(MZ+fjc+k)*JZ] = X[x*JX+y*JY+B_ZM[k]];
				}
			}

			for (int k=0; k<fjc; k++) for (int l=0; l<fjc; l++) for (int m=0; m<fjc; m++) {
				X[k*JX         +l*JY         +m*JZ]         =X[B_X1[k]*JX+B_Y1[l]*JY+B_Z1[m]*JZ];
				X[(MX+fjc+k)*JX+l*JY         +m*JZ]         =X[B_XM[k]*JX+B_Y1[l]*JY+B_Z1[m]*JZ];
				X[k*JX         +(MY+fjc+l)*JY+m*JZ]         =X[B_X1[k]*JX+B_YM[l]*JY+B_Z1[m]*JZ];
				X[k*JX         +l*JY         +(MZ+fjc+m)*JZ]=X[B_X1[k]*JX+B_Y1[l]*JY+B_ZM[m]*JZ];
				X[(MX+fjc+k)*JX+(MY+fjc+l)*JY+m*JZ]         =X[B_XM[k]*JX+B_YM[l]*JY+B_Z1[m]*JZ];
				X[(MX+fjc+k)*JX+l*JY         +(MZ+fjc+m)*JZ]=X[B_XM[k]*JX+B_Y1[l]*JY+B_ZM[m]*JZ];
				X[k*JX         +(MY+fjc+l)*JY+(MZ+fjc+m)*JZ]=X[B_X1[k]*JX+B_YM[l]*JY+B_ZM[m]*JZ];
				X[(MX+fjc+k)*JX+(MY+fjc+l)*JY+(MZ+fjc+m)*JZ]=X[B_XM[k]*JX+B_YM[l]*JY+B_ZM[m]*JZ];
			}
		}
	}
}

void LGrad3::set_M_bounds(Real* X){
if (!debug) cout <<"set_bounds (M) in LGrad3 " << endl;
	int x,y,z;
	int k=0;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			SetBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) {
			//SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++){
				X[x*JX+y*JY+0]     = X[x*JX+y*JY+1];
				X[x*JX+y*JY+MZ+1]  = X[x*JX+y*JY+MZ];
			}
			for (y=1; y<MY+1; y++) for (z=fjc; z<MZ+1; z++)  {
				X[0        +y*JY+z*JZ] = X[1*JX+y*JY+z*JZ];
				X[(MX+1)*JX+y*JY+z*JZ] = X[MX*JX+y*JY+z*JZ];
			}
			for (z=1; z<MZ+1; z++) for (x=1; x<MX+1; x++){
				X[x*JX+0        +z*JZ] = X[x*JX+1*JY+z*JZ];
				X[x*JX+(MY+1)*JY+z*JZ] = X[x*JX+MY*JY+z*JZ];
			}

		} else {
			for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+k] = X[x*JX+y*JY+2*fjc-1-k];
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+MZ+fjc-k-1];
			}
			for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
				for (k=0; k<fjc; k++) X[k*JX+y*JY+z*JZ] = X[(2*fjc-1-k)*JX+y*JY+z*JZ];
				for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = X[(MX+fjc-k-1)*JX+y*JY+z*JZ];
			}
			for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){ //corners?
				for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = X[x*JX+(2*fjc-1-k)*JY+z*JZ];
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = X[x*JX+(MY+fjc-k-1)*JY+z*JZ];
			}
		}
	}
}

void LGrad3::remove_bounds(int *X){
if (!debug) cout <<"remove_bounds (int) in LGrad3 " << endl;
	int x,y,z;
	int k;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			RemoveBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) RemoveBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); else {
			for (x=0; x<MX+2*fjc; x++) for (y=0; y<MY+2*fjc; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+k] = 0;
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = 0;
			}
			for (y=0; y<MY+2*fjc; y++) for (z=0; z<MZ+2*fjc; z++)  {
				for (k=0; k<fjc; k++) X[k*JX+y*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = 0;
			}
			for (z=0; z<MZ+2*fjc; z++) for (x=0; x<MX+2*fjc; x++){//corners?
				for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = 0;
			}
		}
	}
}

void LGrad3::set_bounds(int* X){
if (!debug) cout <<"set_bounds (int) in LGrad3 " << endl;
	int x,y,z;
	int k=0;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			SetBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) {
			//SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);

			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++){
				X[x*JX+y*JY+0]     = X[x*JX+y*JY+BZ1];
				X[x*JX+y*JY+MZ+1]  = X[x*JX+y*JY+BZM];
			}
			for (y=1; y<MY+1; y++) for (z=1; z<MZ+1; z++)  {
				X[0        +y*JY+z*JZ] = X[BX1*JX+y*JY+z*JZ];
				X[(MX+1)*JX+y*JY+z*JZ] = X[BXM*JX+y*JY+z*JZ];
			}
			for (z=1; z<MZ+1; z++) for (x=1; x<MX+1; x++){
				X[x*JX+0        +z*JZ] = X[x*JX+BY1*JY+z*JZ];
				X[x*JX+(MY+1)*JY+z*JZ] = X[x*JX+BYM*JY+z*JZ];
			}

			x=0; {
				for (y=1; y<MY+1; y++){
					X[y*JY+0]     = X[y*JY+BZ1];
					X[y*JY+MZ+1]  = X[y*JY+BZM];
				}
				for (z=1; z<MZ+1; z++){
					X[0        +z*JZ] = X[BY1*JY+z*JZ];
					X[(MY+1)*JY+z*JZ] = X[BYM*JY+z*JZ];
				}
			}
			x=MX+1; {
				for (y=1; y<MY+1; y++){
					X[x*JX+y*JY+0]     = X[x*JX+y*JY+BZ1];
					X[x*JX+y*JY+MZ+1]  = X[x*JX+y*JY+BZM];
				}
				for (z=1; z<MZ+1; z++){
					X[x*JX+0        +z*JZ] = X[x*JX+BY1*JY+z*JZ];
					X[x*JX+(MY+1)*JY+z*JZ] = X[x*JX+BYM*JY+z*JZ];
				}
			}
			y=0; {
				for (x=1; x<MX+1; x++){
					X[x*JX        +0] = X[x*JX+BZ1*JZ];
					X[x*JX+(MZ+1)*JZ] = X[x*JX+BZM*JZ];
				}
			}
			y=MY+1; {
				for (x=1; x<MX+1; x++){
					X[x*JX+y*JY        +0] = X[x*JX+y*JY+BZ1*JZ];
					X[x*JX+y*JY+(MZ+1)*JZ] = X[x*JX+y*JY+BZM*JZ];
				}
			}

			X[0        +0        +0]        =X[BX1*JX+BY1*JY+BZ1*JZ];
			X[(MX+1)*JX+0        +0]        =X[BXM*JX+BY1*JY+BZ1*JZ];
			X[0        +(MY+1)*JY+0]        =X[BX1*JX+BYM*JY+BZ1*JZ];
			X[0        +0        +(MZ+1)*JZ]=X[BX1*JX+BY1*JY+BZM*JZ];
			X[(MX+1)*JX+(MY+1)*JY+0]        =X[BXM*JX+BYM*JY+BZ1*JZ];
			X[(MX+1)*JX+ 0       +(MZ+1)*JZ]=X[BXM*JX+BY1*JY+BZM*JZ];
			X[0        +(MY+1)*JY+(MZ+1)*JZ]=X[BX1*JX+BYM*JY+BZM*JZ];
			X[(MX+1)*JX+(MY+1)*JY+(MZ+1)*JZ]=X[BXM*JX+BYM*JY+BZM*JZ];


		}else {
			for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+k] = X[x*JX+y*JY+B_Z1[k]];
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
			}
			for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
				for (k=0; k<fjc; k++) X[k*JX+y*JY+z*JZ] = X[B_X1[k]*JX+y*JY+z*JZ];
				for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = X[B_XM[k]*JX+y*JY+z*JZ];
			}
			for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){
				for (k=0; k<fjc; k++) X[x*JX+k*JY+z*JZ] = X[x*JX+B_Y1[k]*JY+z*JZ];
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = X[x*JX+B_YM[k]*JY+z*JZ];
			}

			for (x=0; x<fjc; x++ ) {
				for (y=fjc; y<MY+fjc; y++){
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+k]         = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
				}
				for (z=fjc; z<MZ+fjc; z++){
					for (int k=0; k<fjc; k++) X[x*JX+k*JY         +z]  = X[x*JX+B_Y1[k]*JY+z*JZ];
					for (int k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z] = X[x*JX+B_YM[k]*JY+z*JZ];
				}
			}
			for (x=MX+fjc; x<MX+2*fjc; x++) {
				for (y=fjc; y<MY+fjc; y++){
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+k]         = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
				}
				for (z=fjc; z<MZ+fjc; z++){
					for (int k=0; k<fjc; k++) X[x*JX+k*JY        +z]  = X[x*JX+B_Y1[k]*JY+z];
					for (int k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z] = X[x*JX+B_YM[k]*JY+z];
				}
			}
			for (y=0; y<fjc; y++) {
				for (x=fjc; x<MX+fjc; x++){
					for (int k=0; k<fjc; k++) X[x*JX +y*JY       +k]    = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX +y*JY +(MZ+fjc+k)] = X[x*JX+y*JY+B_ZM[k]];
				}
			}
			for (y=MY+fjc; y<MY+2*fjc; y++) {
				for (x=fjc; x<MX+fjc; x++){
					for (int k=0; k<fjc; k++) X[x*JX+y*JY        +k]     = X[x*JX+y*JY+B_Z1[k]];
					for (int k=0; k<fjc; k++) X[x*JX+y*JY+(MZ+fjc+k)*JZ] = X[x*JX+y*JY+B_ZM[k]];
				}
			}

			for (int k=0; k<fjc; k++) for (int l=0; l<fjc; l++) for (int m=0; m<fjc; m++) {
				X[k*JX         +l*JY         +m*JZ]         =X[B_X1[k]*JX+B_Y1[l]*JY+B_Z1[m]*JZ];
				X[(MX+fjc+k)*JX+l*JY         +m*JZ]         =X[B_XM[k]*JX+B_Y1[l]*JY+B_Z1[m]*JZ];
				X[k*JX         +(MY+fjc+l)*JY+m*JZ]         =X[B_X1[k]*JX+B_YM[l]*JY+B_Z1[m]*JZ];
				X[k*JX         +l*JY         +(MZ+fjc+m)*JZ]=X[B_X1[k]*JX+B_Y1[l]*JY+B_ZM[m]*JZ];
				X[(MX+fjc+k)*JX+(MY+fjc+l)*JY+m*JZ]         =X[B_XM[k]*JX+B_YM[l]*JY+B_Z1[m]*JZ];
				X[(MX+fjc+k)*JX+l*JY         +(MZ+fjc+m)*JZ]=X[B_XM[k]*JX+B_Y1[l]*JY+B_ZM[m]*JZ];
				X[k*JX         +(MY+fjc+l)*JY+(MZ+fjc+m)*JZ]=X[B_X1[k]*JX+B_YM[l]*JY+B_ZM[m]*JZ];
				X[(MX+fjc+k)*JX+(MY+fjc+l)*JY+(MZ+fjc+m)*JZ]=X[B_XM[k]*JX+B_YM[l]*JY+B_ZM[m]*JZ];
			}
		}
	}
}

Real LGrad3::ComputeGN(Real* G,int Markov, int M){
	Real GN=0;
	int size;

	if (Markov==2) {
		if (lattice_type==simple_cubic) size =6;
		if (lattice_type==hexagonal) size =12;

		for (int k=0; k<size; k++) GN+=WeightedSum(G+k*M);
		GN/=(1.0*size);
	} else GN=WeightedSum(G);
	return GN;

}
void LGrad3::AddPhiS(Real* phi,Real* Gf,Real* Gb,int Markov, int M){
	Real one=1.0;
	if (Markov ==2) {
		if (lattice_type==hexagonal) {
			for (int k=0; k<12; k++) YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,one/12.0,M);
		} else {
			for (int k=0; k<6; k++) YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,one/6.0,M);
		}
	} else AddTimes(phi,Gf,Gb,M);

}

void LGrad3::AddPhiS(Real* phi,Real* Gf,Real* Gb,Real degeneracy,int Markov, int M){
	if (Markov ==2) {
		if (lattice_type==hexagonal) {
			for (int k=0; k<12; k++) YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,degeneracy/12.0,M);
		} else {
			for (int k=0; k<6; k++) YplusisCtimesAtimesB(phi,Gf+k*M,Gb+k*M,degeneracy/6.0,M);
		}
	} else YplusisCtimesAtimesB(phi,Gf,Gb,degeneracy,M);

}


void LGrad3::AddPhiS(Real* phi,Real* Gf,Real* Gb, Real* G1, Real norm, int Markov, int M){
	cout << "Composition phi Alias not implemented in LGrad3 " << endl;
}

void LGrad3::Initiate(Real* G,Real* Gz,int Markov, int M){
	int size;
	if (Markov==2){
		if (lattice_type==simple_cubic) size =6;
		if (lattice_type==hexagonal) {
			size =12;
			//cout <<"BC in Markov==2 and hexagonal lattice-type is not implemented correctly yet. , do not proceed .." <<endl;
		}
		for (int k=0; k<size; k++) Cp(G+k*M,Gz,M);
	} else Cp(G,Gz,M);
}

void LGrad3::Terminate(Real* Gz,Real* G,int Markov, int M){
	int size;
	if (Markov==2){
		Zero(Gz,M);
		if (lattice_type==simple_cubic) size =6;
		if (lattice_type==hexagonal) {
			size =12;
			//cout <<"BC in Markov==2 and hexagonal lattice-type is not implemented correctly yet. , do not proceed .." <<endl;
		}
		for (int k=0; k<size; k++) Add(Gz,G+k*M,M);
		Norm(Gz,1.0/size,M);
	} else Cp(Gz,G,M);
}

bool LGrad3:: PutMask(int* MASK,vector<int>px,vector<int>py,vector<int>pz,int R){
	bool success=true;
	int length =px.size();
	int X,Y,Z;

	for (int i =0; i<length; i++) {
		int xx,yy,zz;
		xx=px[i]; yy=py[i]; zz=pz[i];
		for (int x=xx-R; x<xx+R+1; x++)
		for (int y=yy-R; y<yy+R+1; y++)
		for (int z=zz-R; z<zz+R+1; z++) {
			if ((xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z) <=R*R) {
				X=x; Y=y; Z=z;
				if (x<1) X+=MX;
				if (y<1) Y+=MY;
				if (z<1) Z+=MZ;
				if (x>MX) X-=MX;
				if (y>MY) Y-=MY;
				if (z>MZ) Z-=MZ;
				MASK[P(X,Y,Z)]++;
			}
		}
		for (int x=1; x<MX; x++)
		for (int y=1; y<MY; y++)
		for (int z=1; z<MZ; z++)
		if (MASK[P(x,y,z)]>1) success=false;
	}

	return success;
}


Real LGrad3::DphiDt(Real* g, Real* B_phitot, Real* phiA, Real* phiB, Real* alphaA, Real* alphaB, Real B_A, Real B_B) {
	cout <<"Grad3 DphiDt not implemented yet " << endl;
	return 0;
}
