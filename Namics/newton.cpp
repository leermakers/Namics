#include "newton.h"

Newton::Newton(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,vector<System*> Sys_,string name_) {
	In=In_; name=name_; Sys=Sys_; Seg=Seg_; Lat=Lat_; Mol=Mol_; 
	KEYS.push_back("method"); KEYS.push_back("e_info"); KEYS.push_back("s_info");
	KEYS.push_back("delta_max"); KEYS.push_back("m"); KEYS.push_back("i_info");
	KEYS.push_back("iterationlimit" ); KEYS.push_back("tolerance"); KEYS.push_back("store_guess"); KEYS.push_back("read_guess");
	KEYS.push_back("stop_criterion"); 
}
Newton::~Newton() {
	delete [] Aij;
	delete [] Ci;
	delete [] Apij; 
#ifdef CUDA
	cudaFree(x);
	cudaFree(x0);
	cudaFree(g);
	cudaFree(xR);
	cudaFree(x_x0);
#else
	delete[] x;
	delete[] x0;
	delete[] g;
	delete[] xR;
	delete[] x_x0;
#endif
}

void Newton::AllocateMemory() {
	iv = Sys[0]->SysMonList.size() * M;	
	if (method=="DIIS-ext") iv +=M;
 
	Aij = new double[m*m]; for (int i=0; i<m*m; i++) Aij[i]=0; 
	Ci = new double[m]; for (int i=0; i<m; i++) Ci[i]=0;
	Apij = new double[m*m]; for (int i=0; i<m*m; i++) Apij[i]=0;
#ifdef CUDA
	x = (double*)AllOnDev(iv); 
	x0 = (double*)AllOnDev(iv);
	g= (double*)AllOnDev(iv);
	xR= (double*)AllOnDev(m*iv);
	x_x0= (double*)AllOnDev(m*iv);

#else
	x = new double[iv]; 
	x0 = new double[iv];
	g = new double[iv];
	xR = new double[m*iv];
	x_x0 =new double[m*iv];
#endif
	Sys[0]->AllocateMemory(); 
	Zero(x,iv);
}

bool Newton::CheckInput() {
	bool success=true;
	string value;
	MX=Lat[0]->MX; MY=Lat[0]->MY; MZ=Lat[0]->MZ; M=(MX+2)*(MY+2)*(MZ+2);
	JX=(MX+2)*(MY+2); JY=(MY+2); 
	success=In[0]->CheckParameters("newton",name,KEYS,PARAMETERS,VALUES);
	if (success) {
		e_info=In[0]->Get_bool(GetValue("e_info"),true); 
		s_info=In[0]->Get_bool(GetValue("s_info"),false); 
		i_info=In[0]->Get_int(GetValue("i_info"),1);
		iterationlimit=In[0]->Get_int(GetValue("iterationlimit"),1000);
		if (iterationlimit < 0 || iterationlimit>1e6) {iterationlimit = 1000; 
		 cout << "Value of 'iterationlimit' out of range 0..1e6, and value set to default value 1000" <<endl; } 

		delta_max=In[0]->Get_double(GetValue("delta_max"),0.1); 
		if (delta_max < 0 || delta_max>100) {delta_max = 0.1;  cout << "Value of delta_max out of range 0..100, and value set to default value 0.1" <<endl; } 
		tolerance=In[0]->Get_double(GetValue("tolerance"),1e-5);
		if (tolerance < 1e-12 ||tolerance>10) {tolerance = 1e-5;  cout << "Value of tolerance out of range 1e-12..10 Value set to default value 1e-5" <<endl; }  
		if (GetValue("method").size()==0) {method="DIIS";} else {
			vector<string>method_options; 
			method_options.push_back("DIIS");
			method_options.push_back("DIIS-ext"); 
			method_options.push_back("Picard"); 
			if (!In[0]->Get_string(GetValue("method"),method,method_options,"In 'newton' the entry for 'method' not recognized: choose from:")) success=false;
		}
		m=In[0]->Get_int(GetValue("m"),10); 
		if (m < 0 ||m>100) {m=10;  cout << "Value of 'm' out of range 0..100, value set to default value 10" <<endl; }
		store_guess=In[0]->Get_bool(GetValue("store_guess"),true);		
		read_guess=In[0]->Get_bool(GetValue("read_guess"),false);
		if (GetValue("stop_criterion").size() > 0) {
			vector<string>options;
			options.push_back("norm_of_g");
			options.push_back("max_of_element_of_|g|");
			if (!In[0]->Get_string(GetValue("stop_criterion"),stop_criterion,options,"In newton the stop_criterion setting was not recognised")) {success=false; };
		}
	}
	return success; 
}


void Newton::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

string Newton::GetValue(string parameter){
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

void Newton::push(string s, double X) {
	doubles.push_back(s);
	doubles_value.push_back(X); 
}
void Newton::push(string s, int X) {
	ints.push_back(s);
	ints_value.push_back(X); 
}
void Newton::push(string s, bool X) {
	bools.push_back(s);
	bools_value.push_back(X); 
}
void Newton::push(string s, string X) {
	strings.push_back(s);
	strings_value.push_back(X); 	
}
void Newton::PushOutput() {
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	doubles.clear();
	doubles_value.clear();
	ints.clear();
	ints_value.clear();  
	push("method",method);
	push("m",m);
	push("delta_max",delta_max);
	push("residual",residual);
	push("tolerance",tolerance);
	push("iterations",it);
	push("iterationlimit",iterationlimit);
	push("stop_criterion",stop_criterion);
}

int Newton::GetValue(string prop,int &int_result,double &double_result,string &string_result){
	int i=0;
	int length = ints.size();
	while (i<length) {
		if (prop==ints[i]) { 
			int_result=ints_value[i];
			return 1;
		}
		i++;
	}
	i=0;
	length = doubles.size();
	while (i<length) {
		if (prop==doubles[i]) { 
			double_result=doubles_value[i];
			return 2;
		}
		i++;
	}
	i=0;
	length = bools.size();
	while (i<length) {
		if (prop==bools[i]) { 
			if (bools_value[i]) string_result="true"; else string_result="false"; 
			return 3;
		}
		i++;
	}
	i=0;
	length = strings.size();
	while (i<length) {
		if (prop==strings[i]) { 
			string_result=strings_value[i]; 
			return 3;
		}
		i++;
	}
	return 0; 
}

bool Newton::PutU() {
	bool success=true;
	int sysmon_length = Sys[0]->SysMonList.size();
	for (int i=0; i<sysmon_length; i++) {
		double *u=Seg[Sys[0]->SysMonList[i]]->u;
		alpha=Sys[0]->alpha; 
		Cp(Seg[Sys[0]->SysMonList[i]]->u,x+i*M,M); 
		if (method == "Picard") Add(u,alpha,M);
		if (method == "DIIS-ext") Add(u,x+sysmon_length*M,M); 
	}
	return success;
}

void Newton::Ax(double* A, double* X, int N){//From Ax_B; below B is not used: it is assumed to contain a row of unities.
	double* U = new double[N*N];
	double* S = new double[N];
	double* VT = new double[N*N];
	integer MM = (integer)N, NN = (integer)N;
	integer LDA=MM, LDU=MM, LDVT=NN, INFO, LWORK;
	int lwork;
	double WKOPT;
	double* WORK;
	char JOBU='S'; //'S' is nodig om alleen de eerste N colommen in U te schrijven.
	char JOBVT='A';

	LWORK = -1; //grootte hulpgeheugen aanvragen
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT, &WKOPT, &LWORK, &INFO );
	lwork = (int)WKOPT;
	WORK = (double*)malloc( lwork*sizeof(double) );
	LWORK = (integer)lwork; //nu uitrekenen.
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT,WORK, &LWORK, &INFO );
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
void Newton::DIIS(double* x, double* x_x0, double* xR, double* Aij, double* Apij,double* Ci, int k, int m, int iv) {
	double normC=0; int posi;
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

bool Newton::PrepareForCalculations() {
	bool success=true;
	return success; 
}
bool Newton::Solve(void) {
	bool success=true;
	if (read_guess){ //this should definitely go to lattice, but I keep it here because for testing...
		double *X= new double[2*M];
		for (int i=0; i<MX+2; i++) for (int j=0; j<MY+2; j++) for (int k=0; k<MZ+2; k++) {
			if (k<MZ/2) {
				X[JX*i+JY*j+k]= -0.38335;
				X[JX*i+JY*j+k+M]= 0.38335;
			} else {
				X[JX*i+JY*j+k]=0;
				X[JX*i+JY*j+k+M]=0;
			}
		}
#ifdef CUDA
		TransferDataToDevice(X,x,2*M);
#else
		Cp(x,X,2*M);
#endif
	}
	double result1,result2;
	result1=Sum(x,2*M);
	result2=Sum(x,2*M);
	if (result1!=result2) cout <<"GPU is not working properly because consecutive calls to Sum to not give the same answer: " << result1 <<" and " << result2 << endl; 
	if (method=="Picard") success=Iterate_Picard(); else success=Iterate_DIIS(); 
	Sys[0]->CheckResults(); 

	return success; 
}

void Newton:: Message(int it, int iterationlimit,double residual, double tolerance) {
	if (it < iterationlimit/10) cout <<"That was easy." << endl;
	if (it > iterationlimit/10 && it < iterationlimit ) cout <<"That will do." << endl;
	if (it <2) cout <<"You hit the nail on the head." << endl; 
	if (residual > tolerance) { cout << "Iterations failed." << endl;
		if (residual < tolerance/10) cout <<"I almost made it..." << endl;
	}
}

bool Newton:: Iterate_Picard() {
	double chi; 
	alpha=Sys[0]->alpha;	
	bool success=true;
	
	int sysmon_length = Sys[0]->SysMonList.size();
	int mon_length = In[0]->MonList.size();	
	it=0;
	cout <<"Picard has been notified" << endl; 
	residual=1;
	while (residual > tolerance && it < iterationlimit) {
		Cp(x0,x,iv); 
		ComputePhis(); 
		Zero(x,iv);
		for (int i=0; i<sysmon_length; i++) for (int k=0; k<mon_length; k++) { 
                        chi= -1.0*Sys[0]->CHI[Sys[0]->SysMonList[i]*mon_length+k];  //The minus sign here is to change the sign of x! just a trick due to properties of PutAlpha where a minus sing is implemented....
			if (chi!=0) PutAlpha(x+i*M,Sys[0]->phitot,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
		}
		for (int i=0; i<sysmon_length; i++) {
			Lat[0]->remove_bounds(x+i*M);
			Times(x+i*M,x+i*M,Sys[0]->KSAM,M);
		}

		Picard(x,x0,delta_max,iv); 
		YisAminB(g,x,x0,iv); 
		for (int i=0; i<sysmon_length; i++) Times(g+i*M,g+i*M,Sys[0]->KSAM,M);  

		residual=Dot(g,g,iv);
		UpdateAlpha(alpha, Sys[0]->phitot, delta_max, M);
		YisAplusC(g,Sys[0]->phitot,-1.0,M); 
		Lat[0]->remove_bounds(g); 

		residual=residual+Dot(g,g,M);
		residual=sqrt(residual); 
		if(it%i_info == 0){
			printf("it = %i g = %1e \n",it,residual);
		}
		it++; 
	}
	Message(it,iterationlimit,residual,tolerance); 
	return success; 
}

bool Newton:: Iterate_DIIS() {
	Zero(x0,iv);
	it=0; k_diis=1; 
	int k=0;
	if (method=="DIIS-ext") ComputeG_ext(); else ComputeG();
	YplusisCtimesX(x,g,delta_max,iv);
	YisAminB(x_x0,x,x0,iv);
	Cp(xR,x,iv);
	residual = sqrt(Dot(g,g,iv));
	printf("DIIS has been notified\n");
	printf("Your guess = %1e \n",residual);
	while (residual > tolerance && it < iterationlimit) {
		it++;
		Cp(x0,x,iv);
		if (method=="DIIS-ext") ComputeG_ext(); else ComputeG();
		k=it % m; k_diis++; //plek voor laatste opslag
		YplusisCtimesX(x,g,-delta_max,iv);
		Cp(xR+k*iv,x,iv); YisAminB(x_x0+k*iv,x,x0,iv);
		DIIS(x,x_x0,xR,Aij,Apij,Ci,k,m,iv);
		residual = sqrt(Dot(g,g,iv));
		if(it%i_info == 0){
			printf("it = %i g = %1e \n",it,residual);
		}
	}

	Message(it,iterationlimit,residual,tolerance); 
	return it<iterationlimit+1;
} 

void Newton::ComputePhis() {
	PutU(); 	
	Sys[0]->ComputePhis();
}

void Newton::ComputeG(){ 
	alpha=Sys[0]->alpha; 
	double chi; 
	ComputePhis();
	int sysmon_length = Sys[0]->SysMonList.size();
	int mon_length = In[0]->MonList.size();	

	Cp(g,x,iv); Zero(alpha,M);
	for (int i=0; i<sysmon_length; i++) {
		for (int k=0; k<mon_length; k++) { 
                        chi= Sys[0]->CHI[Sys[0]->SysMonList[i]*mon_length+k];  
			if (chi!=0) PutAlpha(g+i*M,Sys[0]->phitot,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
		}
	}
	for (int i=0; i<sysmon_length; i++) {
		Add(alpha,g+i*M,M);
	}

	Norm(alpha,1.0/sysmon_length,M);
	for (int i=0; i<sysmon_length; i++) {
		AddG(g+i*M,Sys[0]->phitot,alpha,M);
		Lat[0]->remove_bounds(g+i*M);
		//Times(g+i*M,g+i*M,Sys[0]->KSAM,M);
	} 
}

void Newton::ComputeG_ext(){ 
	alpha=Sys[0]->alpha;
	double chi; 
	ComputePhis();
	int sysmon_length = Sys[0]->SysMonList.size();
	int mon_length = In[0]->MonList.size();	

	Cp(g,x,iv); 
	for (int i=0; i<sysmon_length; i++) { 
		for (int k=0; k<mon_length; k++){ 
                        chi= Sys[0]->CHI[Sys[0]->SysMonList[i]*mon_length+k];  
			if (chi!=0) PutAlpha(g+i*M,Sys[0]->phitot,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
		}
	}
	for (int i=0; i<sysmon_length; i++){
		Lat[0]->remove_bounds(g+i*M);
		Times(g+i*M,g+i*M,Sys[0]->KSAM,M);
	}
	YisAplusC(g+sysmon_length*M,Sys[0]->phitot, -1.0, M);
	Lat[0]->remove_bounds(g+sysmon_length*M);
	Norm(g+sysmon_length*M,-1,M);//you can improve ......
	//Norm(g,-1,sysmon_length*M);
	
}
