#include "system.h"
System::System(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,string name_) {
	Seg=Seg_; Mol=Mol_; Lat=Lat_; In=In_; name=name_; 
if (debug) cout << "Constructor for system " << endl;
	KEYS.push_back("calculation_type");
	KEYS.push_back("generate_guess");
	KEYS.push_back("GPU"); 
	int length = In[0]->MonList.size();
	for (int i=0; i<length; i++) KEYS.push_back("guess-"+In[0]->MonList[i]); 
}
System::~System() {
if (debug) cout << "Destructor for system " << endl;	
	free(H_GrandPotentialDensity);
	free(H_FreeEnergyDensity);
	free(H_alpha);
#ifdef CUDA
	cudaFree(phitot);
	cudaFree(alpha);
	cudaFree(GrandPotentialDensity);
	cudaFree(FreeEnergyDensity);
	cudaFree(TEMP);
	cudaFree(KSAM);
#else	
	free(phitot); 
	free(TEMP);
	free(KSAM);
	free(CHI);
#endif
}

void System::AllocateMemory() {
if (debug) cout << "AllocateMemory in system " << endl; 
	int M=Lat[0]->M;
	//H_GN_A = new double[n_box];
	//H_GN_B = new double[n_box];

	H_GrandPotentialDensity = (double*) malloc(M*sizeof(double));
	H_FreeEnergyDensity = (double*) malloc(M*sizeof(double));
	H_alpha = (double*) malloc(M*sizeof(double));
#ifdef CUDA
	phitot = (double*)AllOnDev(M); 
	alpha=(double*)AllOnDev(M);	
	GrandPotentialDensity =(double*)AllOnDev(M);
	FreeEnergyDensity=(double*)AllOnDev(M);
	TEMP =(double*)AllOnDev(M);	
	KSAM =(int*)AllOnDev(M);
#else
	phitot = (double*) malloc(M*sizeof(double));
	alpha= H_alpha;
	KSAM = (int*) malloc(M*sizeof(int));
	FreeEnergyDensity=H_FreeEnergyDensity;
	GrandPotentialDensity = H_GrandPotentialDensity;
	TEMP = (double*) malloc(M*sizeof(double));
#endif
	Zero(KSAM,M);
	n_mol = In[0]->MolList.size(); 
	int i=0;
	Lat[0]->AllocateMemory(); 
	int n_mon=In[0]->MonList.size();
	while (i<n_mon) {Seg[i]->AllocateMemory(); i++;}
	i=0;
	while (i<n_mol) {Mol[i]->AllocateMemory(); i++;}
}

bool System::PrepareForCalculations() {
if (debug) cout <<"PrepareForCalculations in System " << endl; 
	int M=Lat[0]->M;
	bool success=true;
	FrozenList.clear();
	int length = In[0]->MonList.size();
	for (int i=0; i<length; i++) {if (Seg[i]->freedom == "frozen") FrozenList.push_back(i); }
	if (FrozenList.size()+SysMonList.size()+SysTagList.size() !=In[0]->MonList.size()) {cout <<" There are un-used monomers in system. Remove them before starting" << endl; success=false;}

	Zero(KSAM,M); 
	length=FrozenList.size();
	int i=0;
	while (i<length) {
		int*MASK=Seg[FrozenList[i]]->MASK; 
		Add(KSAM,MASK,M);
		i++;
	}
	length=SysTagList.size();
	i=0;
	while (i<length) {
		int* MASK=Seg[SysTagList[i]]->MASK; 
		Add(KSAM,MASK,M);
		i++;
	}

	Invert(KSAM,KSAM,M); 
	n_mol = In[0]->MolList.size();
	success=Lat[0]->PrepareForCalculations(); 
	int n_mon=In[0]->MonList.size();
	for (int i=0; i<n_mon; i++) {success=Seg[i]->PrepareForCalculations(KSAM);}
	for (int i=0; i<n_mol; i++) {success=Mol[i]->PrepareForCalculations();}
	return success; 
}

string System:: GetMonName(int mon_number_I){
if (debug) cout << "GetMonName for system " << endl;
	return Seg[mon_number_I]->name; 
}

bool System::CheckInput(int start) {
if (debug) cout << "CheckInput for system " << endl;
	bool success=true;	
	bool solvent_found=false; tag_segment=-1; solvent=-1;  //value -1 means no solvent defined. tag_segment=-1; 
	double phibulktot=0;  
	success= In[0]->CheckParameters("sys",name,start,KEYS,PARAMETERS,VALUES);
	if (success) {
		success=CheckChi_values(In[0]->MonList.size()); 	
		GPU=In[0]->Get_bool(GetValue("GPU"),false);
		if (GPU) {if (!cuda) cout << "You should compile the program using the CUDA=1 flag: GPU calculations are impossible; proceed with CPU computations..." << endl; GPU=false; }
		if (Lat[0]->gradients<3) {
			if (GPU) cout <<"GPU support is (for the time being) only available for three-gradient calculations " << endl; 
		}
		if (cuda) {if (!GPU) cout <<" program expect that you are going to use the GPU, but the input is not in line with this (either gradients < 3, or GPU != 'true' : compile without CUDA=1 flag." << endl; success=false;}
		int i=0;
		int length = In[0]->MolList.size(); 
		while (i<length) {
			int j=0; 
			int LENGTH=Mol[i]->MolMonList.size(); 
			while (j<LENGTH) {
				if (!In[0]->InSet(SysMonList,Mol[i]->MolMonList[j])) {
					if (Seg[Mol[i]->MolMonList[j]]->freedom!="tagged") 	
					SysMonList.push_back(Mol[i]->MolMonList[j]);
				}
				if (Seg[Mol[i]->MolMonList[j]]->freedom=="tagged"){
					if (In[0]->InSet(SysTagList,Mol[i]->MolMonList[j])) {
						//cout <<"You can not use the 'tag monomer' " + GetMonName(Mol[i]->MolMonList[j]) + " in more than one molecule." << endl; success=false;  
					} else	SysTagList.push_back(Mol[i]->MolMonList[j]);					
				}
				j++;
			}
			if (Mol[i]->freedom=="free") phibulktot +=Mol[i]->phibulk;
			if (Mol[i]->freedom=="solvent") {solvent_found=true; solvent=i; }
			if (Mol[i]->IsTagged()) { 
				tag_segment = Mol[i]->tag_segment; 
				Mol[i]->n=1.0*Seg[tag_segment]->n_pos;
				Mol[i]->theta= Mol[i]->n*Mol[i]->chainlength;  
			}
			i++; 
		}
		if (!solvent_found ||( phibulktot >0.99999999 && phibulktot < 1.0000000001)) {
			cout << "In system '" + name + "' the 'solvent' was not found, while the volume fractions of the bulk do not add up to unity. " << endl;  
			success=false; 
		}
		//find neutralizer when the charge in the bulk does not add to zero.

		vector<string> options;
		options.push_back("micro_emulsion"); options.push_back("micro_phasesegregation");
		CalculationType="";
		if (GetValue("calculation_type").size()>0) {
			if (!In[0]->Get_string(GetValue("calculation_type"),CalculationType,options," Info about calculation_type rejected") ) {success=false;};
			if (CalculationType=="micro_emulsion") {
				if (In[0]->MolList.size()<3) {cout << "In 'calculation_type : micro-emulsion', we expect at least three types of molecules " << endl; success=false;  }
				int length=In[0]->MolList.size();
				int i=0;
				int number_of_solvents=0;
				int number_of_surfactants=0;
				while (i<length) {
					if (Mol[i]->IsTagged()) {number_of_surfactants++;} else {number_of_solvents++; }
					i++; 
				}
				if (number_of_solvents <2 || number_of_surfactants <1) {
					cout << "In 'calculation_type : micro-emulsion', we expect at least two solvents and one surfactant which is tagged. " << endl;
					success=false; 
				}
			}
			if (CalculationType=="micro_phasesegregation") {
				int length=In[0]->MolList.size();
				bool found=false; 
				int i=0;
				while (i<length && !found) {
					if (Mol[i]->MolMonList.size()>1) found=true;
				}
				if (!found) {cout << "In 'calculation_type : micro_phasesegregation', we expect at least copolymer in the system " << endl; success=false;  }
			}			
		}
		options.clear();
		options.push_back("lamellae"); options.push_back("Im3m"); 
		if (GetValue("generate_guess").size()>0) {
			if (!In[0]->Get_string(GetValue("generate_guess"),GuessType,options," Info about 'generate_guess' rejected") ) {success=false;};
			if (GuessType=="lamellae" || GuessType=="Im3m") {
				int length = In[0]->MonList.size();
				int n_guess=0;
				for (int i=0; i<length; i++) {
					string s="guess-"+In[0]->MonList[i];
					if (GetValue(s).size()>0) {
						Seg[i]->guess_u=In[0]->Get_double(GetValue(s),0); 
						if (Seg[i]->guess_u <-2 || Seg[i]->guess_u >2) {
							cout << "Suggested 'guess value' for 'u' for mon '" + In[0]->MonList[i] + "' out of range: current range is -2, .., 2. Value ignored. " << endl;
							Seg[i]->guess_u = 0;
						} else {
							n_guess++;
							if (i==0 && n_guess==1)  MonA=i;
							if (i==1 && n_guess==2)  MonB=i;
						}
					}
				}
				if (n_guess < 2) {
					cout <<"generate_guess needs two valid non-zero (real) values for 'guess_xxx' quantities. Here xxx is a monomomer name; 'generate_guess : " + GuessType + "' is most likely ineffective.  " << endl;
					cout <<"The idea is that the first 'mon' quantity will be the 'oil' and the second 'mon' quantity the 'water' component. Each filling one half-space. Adjust input accordingly. " << endl;    
				}				
			} 

		}
		
	}
	return success; 
}

void System::PutParameter(string new_param) {
if (debug) cout << "PutParameter for system " << endl;
	KEYS.push_back(new_param); 
}

string System::GetValue(string parameter){
if (debug) cout << "GetValue " + parameter + " for system " << endl;
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

void System::push(string s, double X) {
if (debug) cout << "push (double) for system " << endl;
	doubles.push_back(s);
	doubles_value.push_back(X); 
}
void System::push(string s, int X) {
if (debug) cout << "push (int) for system " << endl;
	ints.push_back(s);
	ints_value.push_back(X); 
}
void System::push(string s, bool X) {
if (debug) cout << "push (bool) for system " << endl;
	bools.push_back(s);
	bools_value.push_back(X); 
}
void System::push(string s, string X) {
if (debug) cout << "push (string) for system " << endl;
	strings.push_back(s);
	strings_value.push_back(X); 	
}
void System::PushOutput() {
if (debug) cout << "PushOutput for system " << endl;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	doubles.clear();
	doubles_value.clear();
	ints.clear();
	ints_value.clear();  
	push("e",e);
	push("k_B",k_B);
	push("eps0",eps0);
	push("temperature",T);
	push("free_energy",FreeEnergy);
	push("grand_potential",GrandPotential);
	push("calculation_type",CalculationType);
	push("guess_type",GuessType);
	push("cuda",cuda);
	push("solvent",Mol[solvent]->name);
	int n_mon=In[0]->MonList.size();
	for(int i=0; i<n_mon; i++) for(int j=0; j<n_mon; j++) if (i!=j) Seg[i]->push("chi-" + Seg[j]->name,CHI[i*n_mon+j]); 
	string s="profile;0"; push("alpha",s);
	s="profile;1"; push("GrandPotentialDensity",s);
	s="profile;2"; push("FreeEnergyDensity",s);
	//"profile;3"; push("KSAM",s);
#ifdef CUDA
	TransferDataToHost(H_alpha, alpha, M);
	TransferDataToHost(H_GrandPotentialDensity, GrandPotentialDensity, M);
	TransferDataToHost(H_FreeEnergyDensity, FreeEnergyDensity, M);
#endif
}
	
double* System::GetPointer(string s){
if (debug) cout << "GetPointer for system " << endl;
	vector<string>sub;
	In[0]->split(s,';',sub);
	if (sub[1]=="0") return H_alpha; 
	if (sub[1]=="1") return H_GrandPotentialDensity;
	if (sub[1]=="2") return H_FreeEnergyDensity;
	// (sub[1]=="3") return H_KSAM;
	return NULL; 
}

int System::GetValue(string prop,int &int_result,double &double_result,string &string_result){
if (debug) cout << "GetValue (long) for system " << endl;
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

bool System::CheckChi_values(int n_seg){
if (debug) cout << "CheckChi_values for system " << endl;
	bool success=true;
	CHI = (double*) malloc(n_seg*n_seg*sizeof(double));
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) {
		CHI[i*n_seg+k]=In[0]->Get_double(Seg[i]->GetValue("chi-"+Seg[k]->name),123);
	}
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) if (CHI[i*n_seg+k] == 123) CHI[i*n_seg+k] = CHI[k*n_seg+i];
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) if (CHI[i*n_seg+k] == 123) CHI[i*n_seg+k] = 0;
	for (int i=0; i<n_seg; i++) for (int k=i+1; k<n_seg; k++) if (CHI[i*n_seg+k] != CHI[k*n_seg+i]) {
		cout <<"CHI-value symmetry violated: chi("<<Seg[i]->name<<","<<Seg[k]->name<<") is not equal to chi("<<Seg[k]->name<<","<<Seg[i]->name<<")"<<endl; success=false;
	}
	for (int i=0; i<n_seg; i++) {if (CHI[i*n_seg+i]!=0) {cout <<"CHI-values for 'same'-segments (e.g. CHI(x,x) ) should be zero. " << endl; success = false;} }
 	return success; 
}

bool System::ComputePhis(){
if(debug) cout <<"ComputePhis in system" << endl;
	int M= Lat[0]->M;
	bool success=true;
	//success=PrepareForCalculations();

	double totphibulk=0;		
	double norm=0; 
	Zero(phitot,M); 
	int length=FrozenList.size();
	for (int i=0; i<length; i++) {
		double *phi_frozen=Seg[FrozenList[i]]->phi;  
		Add(phitot,phi_frozen,M); 
	}

	for (int i=0; i<n_mol; i++) {
		success=Mol[i]->ComputePhi();
	}
	for (int i=0; i<n_mol; i++) {
		norm=0;
		int length=Mol[i]->MolMonList.size();
		if (Mol[i]->freedom=="free") {
			norm=Mol[i]->phibulk/Mol[i]->chainlength; 
			totphibulk +=Mol[i]->phibulk; 
			Mol[i]->n=norm*Mol[i]->GN;
		}

		if (Mol[i]->freedom=="restricted") {
			if (Mol[i]->GN>0) {
			norm = Mol[i]->n/Mol[i]->GN; Mol[i]->phibulk=Mol[i]->chainlength*norm;  totphibulk +=Mol[i]->phibulk; 
			} else { norm = 0; cout <<"GN for molecule " << i << " is not larger than zero..." << endl; }
		}
		if (Mol[i]->IsTagged() || Mol[i]->IsPinned()) {
			if (Mol[i]->GN>0) norm=Mol[i]->n/Mol[i]->GN; else {norm=0; cout <<"GN for molecule " << i << " is not larger than zero..." << endl; }
	
		}
		int k=0;
		Mol[i]->norm=norm; 
		while (k<length) {
			double *phi=Mol[i]->phi+k*M;
			double *G1=Seg[Mol[i]->MolMonList[k]]->G1;
			Div(phi,G1,M); if (norm>0) Norm(phi,norm,M);
if (debug) {
double sum; Sum(sum,phi,M); cout <<"Sumphi in mol " << i << " for mon " << k << ": " << sum << endl; 
}
			k++;
		}

	}
	Mol[solvent]->phibulk=1.0-totphibulk; 
	norm=Mol[solvent]->phibulk/Mol[solvent]->chainlength;
	Mol[solvent]->n=norm*Mol[solvent]->GN;
	Mol[solvent]->theta=Mol[solvent]->n*Mol[solvent]->chainlength;  
	int k=0;
	Mol[solvent]->norm=norm;
	length = Mol[solvent]->MolMonList.size();
	while (k<length) {
		double *phi=Mol[solvent]->phi+k*M;
		if (norm>0) Norm(phi,norm,M);
if (debug) {
double sum; Sum(sum,phi,M); cout <<"Sumphi in mol " << solvent << "for mon " << k << ":" << sum << endl; 
}
		k++;
	}
	for (int i=0; i<n_mol; i++) {

		int length=Mol[i]->MolMonList.size();
		int k=0;
		while (k<length) {
			double* phi_mon=Seg[Mol[i]->MolMonList[k]]->phi;
			double* mol_phitot=Mol[i]->phitot;
			double* phi_molmon = Mol[i]->phi + k*M; 
			Add(phi_mon,phi_molmon,M);
			if (!(Seg[Mol[i]->MolMonList[k]]->freedom == "tagged")) Add(phitot,phi_molmon,M); 
			Add(mol_phitot,phi_molmon,M);
			Seg[Mol[i]->MolMonList[k]]->phibulk +=Mol[i]->fraction(Mol[i]->MolMonList[k])*Mol[i]->phibulk; 
			k++; 
		}
		length=SysTagList.size();
		k=0;
		while (k<length) {
			Cp(Seg[SysTagList[k]]->phi,Seg[SysTagList[k]]->MASK,M); 
			k++;
		}
	}
	int n_seg=In[0]->MonList.size();
	for (int i=0; i<n_seg; i++) {
		if (Seg[i]->freedom !="frozen") Lat[0]->set_bounds(Seg[i]->phi); 
		Lat[0]->Side(Seg[i]->phi_side,Seg[i]->phi,M);
	}
	
	return success;  
}

bool System::CheckResults() {
if (debug) cout << "CheckResults for system " << endl;
	bool success=true;	
	FreeEnergy=GetFreeEnergy();
	GrandPotential=GetGrandPotential();
	CreateMu();
cout <<"free energy     = " << FreeEnergy << endl; 
cout <<"grand potential = " << GrandPotential << endl; 
	int n_mol=In[0]->MolList.size();
	double n_times_mu=0;
	for (int i=0; i<n_mol; i++) {
		double Mu=Mol[i]->Mu;
		double n=Mol[i]->n;
		n_times_mu +=  n*Mu; 
	}
cout <<"free energy     (GP + n*mu) = " << GrandPotential + n_times_mu << endl; 
cout <<"Grand potential (F - n*mu)  = " << FreeEnergy - n_times_mu  << endl; 
	
	for (int i=0; i<n_mol; i++) {
		if (Mol[i]->MolAlList.size()>0) {
			Mol[i]->compute_phi_alias=true;
			Mol[i]->ComputePhi(); 
		}
	}

	return success;  
}

double System::GetFreeEnergy(void) {
if (debug) cout << "GetFreeEnergy for system " << endl;
	int M=Lat[0]->M;
	double* F=FreeEnergyDensity;
	double constant=0;
	int n_mol=In[0]->MolList.size();
	//for (int i=0; i<n_mol; i++) Lat[0]->remove_bounds(Mol[i]->phitot);
	int n_mon=In[0]->MonList.size();
	for (int i=0; i<n_mon; i++) {
		//Lat[0]->remove_bounds(Seg[i]->phi); 
		Lat[0]->remove_bounds(Seg[i]->phi_side);
	}

	Zero(F,M);
	for (int i=0; i<n_mol; i++){
		double n=Mol[i]->n;
		double GN=Mol[i]->GN; 
		int N=Mol[i]->chainlength;
		if (Mol[i]->IsTagged()) N--; //assuming there is just one tagged segment per molecule
		double *phi=Mol[i]->phitot; //contains also the tagged segment 
		constant = log(N*n/GN)/N; 
		Cp(TEMP,phi,M); Norm(TEMP,constant,M); Add(F,TEMP,M); 
	}
	int n_sysmon=SysMonList.size();
	for (int j=0; j<n_sysmon; j++) {
		double* phi=Seg[SysMonList[j]]->phi;
		double* u=Seg[SysMonList[j]]->u;
		Times(TEMP,phi,u,M); Norm(TEMP,-1,M); Add(F,TEMP,M);
	}
	for (int j=0; j<n_mon; j++) for (int k=0; k<n_mon; k++) {
		double chi;
		if (Seg[k]->freedom=="frozen") chi=CHI[j*n_mon+k]; else  chi = CHI[j*n_mon+k]/2;
		double *phi_side = Seg[k]->phi_side;
		double *phi = Seg[j]->phi;
		if (Seg[j]->freedom !="frozen") {Times(TEMP,phi,phi_side,M); Norm(TEMP,chi,M); Add(F,TEMP,M);} 
	}

	for (int i=0; i<n_mol; i++){
		constant=0;
		int n_molmon=Mol[i]->MolMonList.size(); 
		for (int j=0; j<n_molmon; j++) for (int k=0; k<n_molmon; k++) {
			double fA=Mol[i]->fraction(Mol[i]->MolMonList[j]);   
			double fB=Mol[i]->fraction(Mol[i]->MolMonList[k]); 
			if (Mol[i]->IsTagged()) {
				int N=Mol[i]->chainlength; 
				if (N>1) {fA=fA*N/(N-1); fB=fB*N/(N-1); } else {fA =0; fB=0;}
			}  
			double chi = CHI[Mol[i]->MolMonList[j]*n_mon+Mol[i]->MolMonList[k]]/2;
			constant -=fA*fB*chi;
		}
		double* phi=Mol[i]->phitot;
		Cp(TEMP,phi,M); Norm(TEMP,constant,M); Add(F,TEMP,M);
	}
	Lat[0]->remove_bounds(F); Times(F,F,KSAM,M); 
	return Lat[0]->WeightedSum(F);
}

double System::GetGrandPotential(void) {
if (debug) cout << "GetGrandPotential for system " << endl;
	int M=Lat[0]->M;
	double* GP =GrandPotentialDensity;
	int n_mol=In[0]->MolList.size();
	int n_mon=In[0]->MonList.size();
	Zero(GP,M);
	for (int i=0; i<n_mol; i++){
		double *phi=Mol[i]->phitot;
		double phibulk = Mol[i]->phibulk;
		int N=Mol[i]->chainlength;
		if (Mol[i]->IsTagged()) N--; //One segment of the tagged molecule is tagged and then removed from GP through KSAM
		Cp(TEMP,phi,M); YisAplusC(TEMP,TEMP,-phibulk,M); Norm(TEMP,1.0/N,M); //GP has wrong sign. will be corrected at end of this routine; 
		Add(GP,TEMP,M); 
	}
	Add(GP,alpha,M);
	int n_sysmon=SysMonList.size();
	for (int j=0; j<n_sysmon; j++)for (int k=0; k<n_sysmon; k++){
		if (!(Seg[SysMonList[j]]->freedom=="tagged" || Seg[SysMonList[j]]->freedom=="tagged"  )) { //not sure about this line...
		double phibulkA=Seg[SysMonList[j]]->phibulk;
		double phibulkB=Seg[SysMonList[k]]->phibulk;
		double chi = CHI[SysMonList[j]*n_mon+SysMonList[k]]/2; 
		double *phi=Seg[SysMonList[j]]->phi; 
		double *phi_side=Seg[SysMonList[k]]->phi_side; 
		Times(TEMP,phi,phi_side,M); YisAplusC(TEMP,TEMP,-phibulkA*phibulkB,M); Norm(TEMP,chi,M); Add(GP,TEMP,M);
	}
	} 
	Norm(GP,-1.0,M); //correct the sign.
	Lat[0]->remove_bounds(GP); Times(GP,GP,KSAM,M);
	return  Lat[0]->WeightedSum(GP);

}

bool System::CreateMu() {
if (debug) cout << "CreateMu for system " << endl;
	bool success=true;
	double constant; 
	int n_mol=In[0]->MolList.size();
	int n_mon=In[0]->MonList.size();
	for (int i=0; i<n_mol; i++) {
		double Mu=0; 
		double NA=Mol[i]->chainlength; 
		if (Mol[i]->IsTagged()) NA=NA-1;
		double n=Mol[i]->n;
		double GN=Mol[i]->GN;
		Mu=log(NA*n/GN)+1;
		constant=0;
		for (int k=0; k<n_mol; k++) {
			double NB = Mol[k]->chainlength;
			if (Mol[k]->IsTagged()) NB=NB-1; 
			double phibulkB=Mol[k]->phibulk;
			constant +=phibulkB/NB;
		}
		Mu = Mu - NA*constant;

		for (int j=0; j<n_mon; j++) for (int k=0; k<n_mon; k++) {
			double chi= CHI[j*n_mon+k]/2;
			double phibulkA=Seg[j]->phibulk;
			double phibulkB=Seg[k]->phibulk;	
			double Fa=Mol[i]->fraction(j);
			double Fb=Mol[i]->fraction(k); 
			if (Mol[i]->IsTagged()) {Fa=Fa*(NA+1)/(NA); Fb=Fb*(NA+1)/NA;}
			Mu = Mu-NA*chi*(phibulkA-Fa)*(phibulkB-Fb);
		}

		Mol[i]->Mu=Mu; 
//cout <<"mol" << i << " = " << Mu << endl;
	}
	return success; 	
}

//TODO //make sure that the tagged positions are not in frozen range. 


//
//#ifdef CUDA
//TransferDataToHost(H_GN_A,GN_A,n_box);
//TransferDataToHost(H_GN_B,GN_B,n_box);
//#endif
//}

//#ifdef CUDA
//	TransferDataToDevice(H_mask, mask, M*n_box);
//	TransferDataToDevice(H_MASK, MASK, MM);
//	TransferDataToDevice(H_KSAM, KSAM, MM);
//	TransferDataToDevice(H_u, u, MM*n_seg);
//	TransferIntDataToDevice(H_Bx, Bx, n_box);
//	TransferIntDataToDevice(H_By, By, n_box);
//	TransferIntDataToDevice(H_Bz, Bz, n_box);
//	TransferIntDataToDevice(H_Px, Px, n_box);
//	TransferIntDataToDevice(H_Py, Py, n_box);
//	TransferIntDataToDevice(H_Pz, Pz, n_box);
//	if (charges) TransferDataToDevice(H_psi,psi,MM);
//#else
//	Cp(u,H_u,MM*n_seg);
//	if (charges) Cp(psi,H_psi,MM);
//#endif

 
