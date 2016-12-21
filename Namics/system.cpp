

class System {
public:
	System(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,string);

~System();

	string name;
	vector<Input*> In; 
	double* CHI;
	vector<Segment*> Seg; 
	vector<Molecule*> Mol;
	vector<Lattice*> Lat; 
	int M; 
	int MX,MY,MZ;
	vector<int> SysMonList; 
	vector<int> FrozenList;
	vector<int> SysTagList; 
	double FreeEnergy;
	double GrandPotential;
	double* phitot; 
	double* KSAM;
	double* H_KSAM;
	double* GrandPotentialDensity;
	double* FreeEnergyDensity;
	double* alpha;
	double* TEMP;
	bool GPU; 
	int n_mol; 
	int solvent; 
	int tag_segment; 
	bool input_error; 
	string calculation_type; 
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput();
	void PutParameter(string); 
	string GetValue(string); 	
	string GetMonName(int );
	bool CheckChi_values(int);
	void AllocateMemory();
	bool PrepareForCalculations();
	bool ComputePhis();
	bool CheckResults();
	double GetFreeEnergy();
	double GetGrandPotential();
	bool CreateMu();
	
};
System::System(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,string name_) {
	Seg=Seg_; Mol=Mol_; Lat=Lat_; In=In_; name=name_; 
	KEYS.push_back("calculation_type");
	KEYS.push_back("GPU"); 	
}
System::~System() {
}

string System:: GetMonName(int mon_number_I){
	return Seg[mon_number_I]->name; 
}

bool System::CheckInput() {
	bool success=true;
	MX=Lat[0]->MX; 
	MY=Lat[0]->MY;
	MZ=Lat[0]->MZ;
	M=(MX+2)*(MY+2)*(MZ+2); 	
	bool solvent_found=false; tag_segment=-1; solvent=-1;  //value -1 means no solvent defined. tag_segment=-1; 
	double phibulktot=0;  
	success= In[0]->CheckParameters("sys",name,KEYS,PARAMETERS,VALUES);
	if (success) {
		success=CheckChi_values(In[0]->MonList.size()); 

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
		options.push_back("micro-emulsion"); options.push_back("membrane");
		if (GetValue("calculation_type").size()>0) {
			if (!In[0]->Get_string(GetValue("calculation_type"),calculation_type,options," Info about calculation_type rejected") ) {success=false;};
			if (calculation_type=="micro-emulsion") {
				if (In[0]->MolList.size()!=3) {cout << "In 'calculation_type : micro-emulsion', we expect three molecules " << endl; success=false;  }
				int length=In[0]->MolList.size();
				int i=0;
				int number_of_solvents=0;
				int number_of_surfactants=0;
				while (i<length) {
					if (Mol[i]->IsTagged()) {number_of_surfactants++;} else {number_of_solvents++; }
					i++; 
				}
				if (number_of_solvents !=2 || number_of_surfactants !=1) {
					cout << "In 'calculation_type : micro-emulsion', we expect exactly two solvents and one surfactant which is tagged. " << endl;
					success=false; 
				}
			}
		}
	}
	return success; 
}

void System::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

string System::GetValue(string parameter){
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

bool System::CheckChi_values(int n_seg){
	bool success=true;
	CHI = new double[n_seg*n_seg]; 
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

void System::AllocateMemory() {
//define on CPU
	H_KSAM=new double[M];
#ifdef CUDA
//define on GPU
	phitot = (double*)AllOnDev(M); 
	KSAM=(double*)AllOnDev(M);
	alpha=(double*)AllOnDev(M);
	FreeEnergyDensity=(double*)AllOnDev(M);
	GrandPotentialDensity =(double*)AllOnDev(M);
	TEMP =(double*)AllOnDev(M);	
#else
	phitot = new double[M]; 
	KSAM = H_KSAM;
	alpha=new double[M];
	FreeEnergyDensity=new double[M];
	GrandPotentialDensity =new double[M];
	TEMP =new double[M];	
#endif
	n_mol = In[0]->MolList.size(); 
	int i=0;
	Lat[0]->AllocateMemory(); 
	int n_mon=In[0]->MonList.size();
	while (i<n_mon) {Seg[i]->AllocateMemory(); i++;}
	i=0;
	while (i<n_mol) {Mol[i]->AllocateMemory(); i++;}
}

bool System::PrepareForCalculations() {
	bool success=true;

	FrozenList.clear();
	int length = In[0]->MonList.size();
	for (int i=0; i<length; i++) {if (Seg[i]->freedom == "frozen") FrozenList.push_back(i); }
	if (FrozenList.size()+SysMonList.size()+SysTagList.size() !=In[0]->MonList.size()) {cout <<" There are un-used monomers in system. Remove them before starting" << endl; success=false;}

	Zero(KSAM,M); 

	length=FrozenList.size();
	int i=0;
	while (i<length) {
		double*MASK=Seg[FrozenList[i]]->MASK; 
		Add(KSAM,MASK,M);
		i++;
	}
	length=SysTagList.size();
	i=0;
	while (i<length) {
		double* MASK=Seg[SysTagList[i]]->MASK; 
		Add(KSAM,MASK,M);
		i++;
	}
	invert(KSAM,KSAM,M); 
	Lat[0]->remove_bounds(KSAM); 
	n_mol = In[0]->MolList.size(); 
	success=Lat[0]->PrepareForCalculations(); 
	int n_mon=In[0]->MonList.size();
	for (int i=0; i<n_mon; i++) {success=Seg[i]->PrepareForCalculations(KSAM);}
	for (int i=0; i<n_mol; i++) {success=Mol[i]->PrepareForCalculations();}
	return success; 
}


bool System::ComputePhis(){
	bool success=true;
	success=PrepareForCalculations();
	double totphibulk=0;		
	double norm=0; 
	Zero(phitot,M); 
	int length=FrozenList.size();
	for (int i=0; i<length; i++) {
		double *phi_frozen=Seg[FrozenList[i]]->GetPhi();  
		Add(phitot,phi_frozen,M); 
	}
	for (int i=0; i<n_mol; i++) {
		success=Mol[i]->ComputePhi();
	}
	for (int i=0; i<n_mol; i++) {
		norm=0;
		int length=Mol[i]->MolMonList.size();
		if (Mol[i]->freedom=="free") {norm=Mol[i]->phibulk/Mol[i]->chainlength; totphibulk +=Mol[i]->phibulk; }
		if (Mol[i]->freedom=="restricted") {norm = Mol[i]->n/Mol[i]->GN; Mol[i]->phibulk=Mol[i]->chainlength*norm;  totphibulk +=Mol[i]->phibulk; }
		if (Mol[i]->IsTagged() || Mol[i]->IsPinned()) {norm=Mol[i]->n/Mol[i]->GN;}
		int k=0;
		while (k<length) {
			double *phi=Mol[i]->phi+k*M;
			double *G1=Seg[Mol[i]->MolMonList[k]]->G1;
			Div(phi,G1,M); if (norm>0) Norm(phi,norm,M);
//cout <<"in mol " << i << "sum phi is " << Sum(phi,M) << endl; 
			k++;
		}

	}
	Mol[solvent]->phibulk=1.0-totphibulk; 
	norm=Mol[solvent]->phibulk/Mol[solvent]->chainlength;
	Mol[solvent]->n=norm*Mol[solvent]->GN; 
	int k=0;
	length = Mol[solvent]->MolMonList.size();
	while (k<length) {
		double *phi=Mol[solvent]->phi+k*M;
		if (norm>0) Norm(phi,norm,M);
//cout <<"in solvent " << Sum(phi,M) << endl; 
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
		Lat[0]->Side(Seg[i]->phi_side,Seg[i]->phi,M); 
	}
	
	return success;  
}

bool System::CheckResults() {
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

	return success;  
}

double System::GetFreeEnergy(void) {
	double* F=FreeEnergyDensity;
	double constant=0;
	int n_mol=In[0]->MolList.size();
	for (int i=0; i<n_mol; i++) Lat[0]->remove_bounds(Mol[i]->phitot);
	int n_mon=In[0]->MonList.size();
	for (int i=0; i<n_mon; i++) {Lat[0]->remove_bounds(Seg[i]->phi); Lat[0]->remove_bounds(Seg[i]->phi_side);}

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
		double chi = CHI[j*n_mon+k]/2;
		double *phi_side = Seg[k]->phi_side;
		double *phi = Seg[j]->phi;
		Times(TEMP,phi,phi_side,M); Norm(TEMP,chi,M); Add(F,TEMP,M); 
	}

	for (int i=0; i<n_mol; i++){
		constant=0;
		int n_molmon=Mol[i]->MolMonList.size(); 
		for (int j=0; j<n_molmon; j++) for (int k=0; k<n_molmon; k++) {
			double fA=Mol[i]->fraction(Mol[i]->MolMonList[j]);   
			double fB=Mol[i]->fraction(Mol[i]->MolMonList[k]); 
			if (Mol[i]->IsTagged()) {int N=Mol[i]->chainlength; fA=fA*N/(N-1); fB=fB*N/(N-1); }  
			double chi = CHI[Mol[i]->MolMonList[j]*n_mon+Mol[i]->MolMonList[k]]/2;
			constant -=fA*fB*chi;
		}
		double* phi=Mol[i]->phitot;
		Cp(TEMP,phi,M); Norm(TEMP,constant,M); Add(F,TEMP,M);
	}
	Lat[0]->remove_bounds(F); Times(F,F,KSAM,M); 
	return Sum(F,M);
}

double System::GetGrandPotential(void) {
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
	return Sum(GP,M); 

}

bool System::CreateMu() {
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


 
