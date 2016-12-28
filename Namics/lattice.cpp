#include "lattice.h"
Lattice::Lattice(vector<Input*> In_,string name_) {
	In=In_; name=name_; 
	KEYS.push_back("n_layers_x");   KEYS.push_back("n_layers_y"); KEYS.push_back("n_layers_z");
	KEYS.push_back("lowerbound_x"); KEYS.push_back("upperbound_x");
	KEYS.push_back("lowerbound_y"); KEYS.push_back("upperbound_y");
	KEYS.push_back("lowerbound_z"); KEYS.push_back("upperbound_z");
 	KEYS.push_back("bond_length");  KEYS.push_back("lattice_type");  
}
Lattice::~Lattice() {
}

bool Lattice::CheckInput() {
	bool success;
	string Value;
	vector<string> options;
	success = In[0]->CheckParameters("lat",name,KEYS,PARAMETERS,VALUES);
	if (success){
		if (!In[0]->Get_int(GetValue("n_layers_x"),MX,1,1e6,"In 'lat' the parameter 'n_layers_x' is required")) {success=false;}
		if (!In[0]->Get_int(GetValue("n_layers_y"),MY,1,1e6,"In 'lat' the parameter 'n_layers_y' is required")) {success=false;}
		if (!In[0]->Get_int(GetValue("n_layers_z"),MZ,1,1e6,"In 'lat' the parameter 'n_layers_z' is required")) {success=false;}
		Volume=MX*MY*MZ;
		JX=(MX+2)*(MY+2); JY=(MY+2); M = (MX+2)*(MY+2)*(MZ+2);   
		Mx=In[0]->Get_int(GetValue("sub_box_size"),MX/2);
		if (Mx>MX) {cout << "'sub_box_size' can not exceed the size of the main box." << endl; success=false;} else My=Mz=Mx;
		bond_length =  In[0]->Get_int(GetValue("bond_length"),5e-10);
		if (bond_length < 0 || bond_length > 1e-8) {cout <<" bond_length out of range 0..1e-8 " << endl; success=false;}

		options.push_back("simple_cubic"); options.push_back("hexagonal");
		Value=GetValue("lattice_type"); if (Value.length()>0) {
			In[0]->Get_string(Value,lattice_type,options,"Input for 'lattice_type' not recognized."); 
		} else {
			lattice_type  = "simple_cubic";
		} 

		options.clear(); 
		options.push_back("mirror_1"); options.push_back("mirror_2"); options.push_back("periodic"); 
		//options.push_back("shifted_mirror");
		string VALUE1,VALUE2,VALUE3,VALUE4,VALUE5,VALUE6;
		Value.clear();
		Value=GetValue("lowerbound_x"); if (Value.length()>0) {
			In[0]->Get_string(Value,VALUE1,options,"for 'lowerbound_x' boundary condition not recognized.");
			if (VALUE1=="mirror_1") BX1=1; 
			if (VALUE1=="mirror_2") BX1=2;
			if (VALUE1=="periodic") BX1=MX;
		} else {
			VALUE1="periodic";
			BX1=MX;
		} 
		Value.clear();	
		Value=GetValue("lowerbound_y"); if (Value.length()>0) {
			In[0]->Get_string(Value,VALUE2,options,"for 'lowerbound_y' boundary condition not recognized.");
			if (VALUE2=="mirror_1") BY1=1; 
			if (VALUE2=="mirror_2") BY1=2;
			if (VALUE2=="periodic") BY1=MY; 
		} else {
			VALUE2="periodic";
			BY1=MY;
		} 
		Value.clear();	
		Value=GetValue("lowerbound_z"); if (Value.length()>0) {
			In[0]->Get_string(Value,VALUE3,options,"for 'lowerbound_z' boundary condition not recognized.");
			if (VALUE3=="mirror_1") BZ1=1; 
			if (VALUE3=="mirror_2") BZ1=2;
			if (VALUE3=="periodic") BZ1=MZ;  
		} else {
			VALUE3="periodic";
			BZ1=MZ;
		} 
		Value.clear();	
		Value=GetValue("upperbound_x"); if (Value.length()>0) {
			In[0]->Get_string(Value,VALUE4,options,"for 'upperbound_x' boundary condition not recognized."); 
			if (VALUE4=="mirror_1") BXM=MX-1; 
			if (VALUE4=="mirror_2") BXM=MX-2;
			if (VALUE4=="periodic") BXM=1;
		} else {
			VALUE4="periodic";
			BXM=1;
		} 
		Value.clear();	
		Value=GetValue("upperbound_y"); if (Value.length()>0) {
			In[0]->Get_string(Value,VALUE5,options,"for 'upperbound_y' boundary condition not recognized."); 
			if (VALUE5=="mirror_1") BYM=MY-1; 
			if (VALUE5=="mirror_2") BYM=MY-2;
			if (VALUE5=="periodic") BYM=1;
		} else {
			VALUE5="periodic";
			BYM=1;
		} 
		Value.clear();	
		Value=GetValue("upperbound_z"); if (Value.length()>0) { 
			In[0]->Get_string(Value,VALUE6,options,"for 'upperbound_z' boundary condition not recognized."); 
			if (VALUE6=="mirror_1") BZM=MZ-1; 
			if (VALUE6=="mirror_2") BZM=MZ-2;
			if (VALUE6=="periodic") BZM=1;
		} else {
			VALUE6="periodic";
			BZM=1;
		} 	
		if (VALUE1 != VALUE4) {cout <<"In x-direction the boundary conditions do not match:" + VALUE1 << " and " <<  VALUE4 << endl; success=false;}
		if (VALUE2 != VALUE5) {cout <<"In y-direction the boundary conditions do not match:" + VALUE2 << " and " <<  VALUE5 << endl; success=false;}
		if (VALUE3 != VALUE6) {cout <<"In z-direction the boundary conditions do not match:" + VALUE3 << " and " <<  VALUE6 << endl; success=false;}
			
	}
	return success;
}

void Lattice::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

string Lattice::GetValue(string parameter){
	int i=0;
	int length = PARAMETERS.size();
	while (i<length) {
		if (parameter==PARAMETERS[i]) {
			return VALUES[i]; 
		}
		i++;
	}
	return "" ; 
}

void Lattice::AllocateMemory(void) {
	
}
bool Lattice::PrepareForCalculations(void) {
	bool success=true;
	return success; 
}

void Lattice::push(string s, double X) {
	doubles.push_back(s);
	doubles_value.push_back(X); 
}
void Lattice::push(string s, int X) {
	ints.push_back(s);
	ints_value.push_back(X); 
}
void Lattice::push(string s, bool X) {
	bools.push_back(s);
	bools_value.push_back(X); 
}
void Lattice::push(string s, string X) {
	strings.push_back(s);
	strings_value.push_back(X); 	
}
void Lattice::PushOutput() {
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	doubles.clear();
	doubles_value.clear();
	ints.clear();
	ints_value.clear();  
	push("n_layers_x",MX);
	push("n_layers_y",MY);
	push("n_layers_z",MZ);
	string mirror1="mirror_1"; 
	string periodic="periodic";
	string mirror2="mirror_2";
	string s; 
	if (BX1==1) push("lowerbound_x",mirror1);
	if (BXM==MX-1) push("upperbound_x",mirror1);
	if (BX1==2) push("lowerbound_x",mirror2);
	if (BXM==MX-2) push("upperbound_x",mirror2);
	if (BX1==MX) push("lowerbound_x",periodic);
	if (BXM==1) push("upperbound_x",periodic);
	if (BY1==1) push("lowerbound_y",mirror1);
	if (BYM==MY-1) push("upperbound_y",mirror1);
	if (BY1==2) push("lowerbound_y",mirror2);
	if (BYM==MY-2) push("upperbound_y",mirror2);
	if (BY1==MY) push("lowerbound_y",periodic);
	if (BYM==1) push("upperbound_y",periodic);	
	if (BZ1==1) push("lowerbound_z",mirror1);
	if (BZM==MZ-1) push("upperbound_z",mirror1);
	if (BZ1==2) push("lowerbound_z",mirror2);
	if (BZM==MZ-2) push("upperbound_z",mirror2);
	if (BZ1==MZ) push("lowerbound_z",periodic);
	if (BZM==1) push("upperbound_z",periodic);
	push("volume",MX*MY*MZ); 
	push("lattice_type",lattice_type);
	push("bond_length",bond_length); 
}

void Lattice::Side(double *X_side, double *X, int M) {
	Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Norm(X_side,1.0/6.0,M);
}

void Lattice::propagate(double *G, double *G1, int s_from, int s_to) { 
	double *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	Zero(gs,M);
	SetBoundaries(gs_1,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(gs+JX,gs_1,M-JX); Add(gs,gs_1+JX,M-JX);
	Add(gs+JY,gs_1,M-JY); Add(gs,gs_1+JY,M-JY);
	Add(gs+1,gs_1,M-1);  Add(gs,gs_1+1, M-1);
	Norm(gs,1.0/6.0,M); Times(gs,gs,G1,M);
}

void Lattice::remove_bounds(double *X){ 
	RemoveBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
}
 
void Lattice::set_bounds(double *X){  
	SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
}

/* 

All below is commented out. This is here to recover from first program. 

void Sideh(double *X_side, double *X, int M) {
	Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);

	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Add(X_side+JX,X+JY,MM-JX-JY); Add(X_side+JY,X+JX,MM-JY-JX);
	Add(X_side+1,X+JY,MM-JY-1); Add(X_side+JY,X+1,MM-JY-1);
	Add(X_side+1,X+JX,MM-JX-1); Add(X_side+JX,X+1,MM-JX-1);
	Norm(X_side,1.0/12.0,M);
}



void Side(double *X_side, double *X, int M) {
	Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Norm(X_side,1.0/6.0,M);
}

void advanced_average(double *X_side, double *X, int M){
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

void Propagateh(double* G, double* G1, int s_from, int s_to) { //on small boxes
	int MMM=M*n_box;
	double *gs = G+MMM*s_to, *gs_1 = G+MMM*s_from, *g = G1;
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

void PROPAGATEh(double *G, double *G1, int s_from, int s_to) { //on big box
	double *gs = G+MM*(s_to), *gs_1 = G+MM*(s_from), *g = G1;
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



void Propagate(double* G, double* G1, int s_from, int s_to) { //on small boxes
	int MMM=M*n_box;
	double *gs = G+MMM*s_to, *gs_1 = G+MMM*s_from, *g = G1;
	Zero(gs,MMM);
	for (int p=0; p<n_box; p++) SetBoundaries(gs_1+M*p,jx,jy,bx1,bxm,by1,bym,bz1,bzm,Mx,My,Mz);
	Add(gs+jx,gs_1,MMM-jx); Add(gs,gs_1+jx,MMM-jx);
	Add(gs+jy,gs_1,MMM-jy); Add(gs,gs_1+jy,MMM-jy);
	Add(gs+1,gs_1,MMM-1);  Add(gs,gs_1+1, MMM-1);
	Norm(gs,1.0/6.0,MMM); Times(gs,gs,g,MMM);
}

void PROPAGATE(double *G, double *G1, int s_from, int s_to) { //on big box
	double *gs = G+MM*(s_to), *gs_1 = G+MM*(s_from), *g = G1;
	Zero(gs,MM);
	SetBoundaries(gs_1,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(gs+JX,gs_1,MM-JX); Add(gs,gs_1+JX,MM-JX);
	Add(gs+JY,gs_1,MM-JY); Add(gs,gs_1+JY,MM-JY);
	Add(gs+1,gs_1,MM-1);  Add(gs,gs_1+1, MM-1);
	Norm(gs,1.0/6.0,MM); Times(gs,gs,g,MM);
}

#ifdef CUDA
void DistributeG1(double *G1, double *g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	//int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
	//distributeg1<<<n_blocks,block_size>>>(G1,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
        dim3 blocks(ceil(Mx/8.0),ceil(My/8.0),ceil(Mz/8.0));
        dim3 blockdim(8,8,8);
        distributeg1<<<blocks,blockdim>>>(G1,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);

}
#else
void DistributeG1(double *G1, double *g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
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
void CollectPhi(double* phi, double* GN, double* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
 	 dim3 blocks(ceil(Mx/8.0),ceil(My/8.0),ceil(Mz/8.0));
        dim3 blockdim(8,8,8);
        collectphi<<<blocks,blockdim>>>(phi,GN,rho,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
}
#else
void CollectPhi(double* phi, double* GN, double* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int pos_l=-M;
	int pos_x,pos_y,pos_z;
	int Bxp,Byp,Bzp;
	double Inv_H_GNp;
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

//#ifdef CUDA
//		if (charges) {
//			phi_na =(double*)AllOnDev(MM);
//			phi_cl =(double*)AllOnDev(MM);
//			psi_0 =(double*)AllOnDev(MM);
//			psi_side =(double*)AllOnDev(MM);
//			psi=x+MM;
//			q =(double*)AllOnDev(MM);
//		}
#else
//		if (charges) {
//			psi=x+MM;
//			q = new double[MM];
//			psi_0 = new double[MM];
//			psi_side = new double[MM];
//		}
//#endif
*/
