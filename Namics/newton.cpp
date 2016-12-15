class Newton {
public:
	Newton(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,vector<System*>,string);

~Newton();

	string name;
	vector<Input*> In; 
	vector<System*> Sys;
	vector<Segment*> Seg;
	vector<Lattice*> Lat;  
	vector<Molecule*> Mol;
	int iterationlimit,m,i_info;
	int k_diis,it; 
	double tolerance,delta_max;
	double residual;
	bool e_info,s_info; 
	string method;
	bool store_guess;
	bool read_guess;  
	string stop_criterion;
	int iv; 
	int M,MX,MY,MZ,JX,JY;
	double* x;
	double* u;
	double* x0;
	double* g;
	double* xR;
	double* x_x0;
	double* alpha;
	double* Aij;
	double* Ci;
	double* Apij; 
//if properties are added, also read them form input and/or set the default. See CheckInput() below. 
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput();
	void PutParameter(string); 
	string GetValue(string); 
	bool Solve();
	void AllocateMemory(); 
	bool PrepareForCalculations(void);
	void Ax(double* , double* , int );
	void DIIS(double* , double* , double* , double*, double* ,double* , int , int , int );
	void ComputeG(); 
	void ComputeG_alpha();

};
Newton::Newton(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,vector<System*> Sys_,string name_) {
	In=In_; name=name_; Sys=Sys_; Seg=Seg_; Lat=Lat_; Mol=Mol_; 
	KEYS.push_back("method"); KEYS.push_back("e_info"); KEYS.push_back("s_info");
	KEYS.push_back("delta_max"); KEYS.push_back("m"); KEYS.push_back("i_info");
	KEYS.push_back("iterationlimit" ); KEYS.push_back("tolerance"); KEYS.push_back("store_guess"); KEYS.push_back("read_guess");
	KEYS.push_back("stop_criterion"); 
}
Newton::~Newton() {
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
			vector<string>method_options; method_options.push_back("DIIS");method_options.push_back("DIIS-alpha");
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

void Newton::AllocateMemory() {
	iv = Sys[0]->SysMonList.size() * M;	
//define on CPU

#ifdef CUDA
//define on GPU
	x = (double*)AllOnDev(iv); 
	x0 = (double*)AllOnDev(iv);
	g= (double*)AllOnDev(iv);
	xR= (double*)AllOnDev(m*iv);
	x_x0= (double*)AllOnDev(m*iv);
	alpha= (double*)AllOnDev(M);

#else
	x = new double[iv]; 
	x0 = new double[iv];
	g = new double[iv];
	xR = new double[m*iv];
	x_x0 =new double[m*iv];
	alpha = new double[M];
	Aij = new double[m*m]; for (int i=0; i<m*m; i++) Aij[i]=0; 
	Ci = new double[m]; for (int i=0; i<m; i++) Ci[i]=0;
	Apij = new double[m*m]; for (int i=0; i<m*m; i++) Apij[i]=0;
#endif
	Sys[0]->AllocateMemory(); 
	Zero(x,iv);
	u=Seg[0]->u;
//do here initial guess and put it in x. 
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
 	success=Sys[0]->PrepareForCalculations();
	return success; 
}
bool Newton::Solve(void) {
	if (read_guess){ 
		for (int i=0; i<MX+2; i++) for (int j=0; j<MY+2; j++) for (int k=0; k<MZ+2; k++) {
			if (k<MZ/2) {
				x[JX*i+JY*j+k]= -0.38335;
				x[JX*i+JY*j+k+M]= 0.38335;
			} else {
				x[JX*i+JY*j+k]=0;
				x[JX*i+JY*j+k+M]=0;
			}
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
//	Cp(u,H_u,MM*n_seg);
//	if (charges) Cp(psi,H_psi,MM);
#endif
	Zero(x0,iv);
	it=0; k_diis=1; k=0;
	Sys[0]->PutU(x);
	if (method=="DIIS-alpha") ComputeG_alpha(); else ComputeG();
	YplusisCtimesX(x,g,delta_max,iv);
	YisAminB(x_x0,x,x0,iv);
	Cp(xR,x,iv);
	residual = sqrt(Dot(g,g,iv));
	printf("DIIS has been notified\n");
	printf("Your guess = %1e \n",residual);
	while (residual > tolerance && it < iterationlimit) {
		it++;
		Cp(x0,x,iv);
		Sys[0]->PutU(x);
		if (method=="DIIS-alpha") ComputeG_alpha(); else ComputeG();
		k=it % m; k_diis++; //plek voor laatste opslag
		YplusisCtimesX(x,g,-delta_max,iv);
		Cp(xR+k*iv,x,iv); YisAminB(x_x0+k*iv,x,x0,iv);
		DIIS(x,x_x0,xR,Aij,Apij,Ci,k,m,iv);
		Cp(u,x,M); 
		residual = sqrt(Dot(g,g,iv));
		if(it%i_info == 0){
			printf("it = %i g = %1e \n",it,residual);
		}
	}
	cout <<"That will do. " << endl;

	
	//return Helmholtz();


	//success=Sys[0]->ComputePhis(); 
	return it<iterationlimit+1;
} 

void Newton::ComputeG(){ 
#ifdef CUDA
	Zero(phi,4*MM);
#endif 
	Sys[0]->ComputePhis();
	int sysmon_length = Sys[0]->SysMonList.size();
	int mon_length = In[0]->MonList.size();	


	for (int k=0; k<mon_length; k++)	Lat[0]->Side(Seg[k]->phi_side,Seg[k]->phi,M);
	Cp(g,x,iv); Zero(alpha,M);

	for (int i=0; i<sysmon_length; i++) 
		for (int k=0; k<mon_length; k++) 
			PutAlpha(g+k*M,Sys[0]->phitot,Seg[k]->phi_side,Sys[0]->CHI[Sys[0]->SysMonList[i]*mon_length+k],Seg[k]->phibulk,M);

	for (int i=0; i<sysmon_length; i++) Add(alpha,g+i*M,M);

	Norm(alpha,1.0/sysmon_length,M);
	for (int i=0; i<sysmon_length; i++) {
		AddG(g+i*M,Sys[0]->phitot,alpha,M);
		Lat[0]->remove_bounds(g+i*M);
		Times(g+i*M,g+i*M,Sys[0]->KSAM,M);
	} 
}

void Newton::ComputeG_alpha(){ 

#ifdef CUDA
	Zero(phi,4*MM);
#endif
	int sysmon_length = Sys[0]->SysMonList.size();
	int mon_length = In[0]->MonList.size();	
		
	for (int k=0; k<mon_length; k++) Seg[k]->alpha=alpha; 
	
	Sys[0]->ComputePhis();

	YplusisCtimesX(alpha, Sys[0]->phitot, delta_max, M);

	for (int k=0; k<mon_length; k++)	Lat[0]->Side(Seg[k]->phi_side,Seg[k]->phi,M);
	Cp(g,x,iv); 
	for (int i=0; i<sysmon_length; i++) 
		for (int k=0; k<mon_length; k++) 
			PutAlpha(g+k*M,Sys[0]->phitot,Seg[k]->phi_side,Sys[0]->CHI[Sys[0]->SysMonList[i]*mon_length+k],Seg[k]->phibulk,M);
	for (int i=0; i<sysmon_length; i++){
		Lat[0]->remove_bounds(g+i*M);
		Times(g+i*M,g+i*M,Sys[0]->KSAM,M);
	}
}