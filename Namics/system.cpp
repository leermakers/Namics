 
class System {
public:
	System(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,string);

~System();

	string name;
	vector<Input*> In; 
	float* CHI;
	vector<Segment*> Seg; 
	vector<Molecule*> Mol;
	vector<Lattice*> Lat; 
	int M; 
	int MX,MY,MZ;
	vector<int> SysMonList; 
	vector<int> FrozenList;
	vector<int> SysTagList; 
	float* phitot; 
	float* KSAM;
	float* H_KSAM;
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
	bool PutU(float*); 

	

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
	float phibulktot=0;  
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
						cout <<"You can not use the 'tag monomer' " + GetMonName(Mol[i]->MolMonList[j]) + " in more than one molecule." << endl; success=false;  
					} else	SysTagList.push_back(Mol[i]->MolMonList[j]);					
				}
				j++;
			}
			if (Mol[i]->freedom=="free") phibulktot +=Mol[i]->phibulk;
			if (Mol[i]->freedom=="solvent") {solvent_found=true; solvent=i; }
			if (Mol[i]->IsTagged()) { 
				tag_segment = Mol[i]->tag_segment; 
cout <<"tag_segment" << tag_segment << " " << Seg[tag_segment]->n_pos<< endl; 
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
				if (!In[0]->MolList.size()==3) {cout << "In 'calculation_type : micro-emulsion', we expect three molecules " << endl; success=false;  }
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
	CHI = new float[n_seg*n_seg]; 
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) {
		CHI[i*n_seg+k]=In[0]->Get_float(Seg[i]->GetValue("chi-"+Seg[k]->name),123);
	}
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) if (CHI[i*n_seg+k] == 123) CHI[i*n_seg+k] = CHI[k*n_seg+i];
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) if (CHI[i*n_seg+k] == 123) CHI[i*n_seg+k] = 0;
	for (int i=0; i<n_seg; i++) for (int k=i+1; k<n_seg; k++) if (CHI[i*n_seg+k] != CHI[k*n_seg+i]) {
		cout <<"CHI-value symmetry violated: chi("<<Seg[i]->name<<","<<Seg[k]->name<<") is not equal to chi("<<Seg[k]->name<<","<<Seg[i]->name<<")"<<endl; success=false;
	}
	for (int i=0; i<n_seg; i++) {if (CHI[i*n_seg+i]!=0) {cout <<"CHI-values for 'same'-segments (e.g. CHI(x,x) ) should be zero. " << endl; success = false;} }
 	return success; 
}

bool System::PutU(float* x) {
	bool success=true;
	int length=SysMonList.size();
	int i=0;
	while (i<length) {
		success=Seg[SysMonList[i]]->PutU(x+M*i);
		i++;
	}
	return success;
}

void System::AllocateMemory() {
//define on CPU
	H_KSAM=new float[M];
#ifdef CUDA
//define on GPU
	phitot = (float*)AllOnDev(M); 
	KSAM=(float*)AllOnDev(M);	
#else
	phitot = new float[M]; 
	KSAM = H_KSAM;
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
		float*MASK=Seg[FrozenList[i]]->MASK; 
		Add(KSAM,MASK,M);
		i++;
	}
	length=SysTagList.size();
	i=0;
	while (i<length) {
		float* MASK=Seg[SysTagList[i]]->MASK; 
		Add(KSAM,MASK,M);
		i++;
	}
	invert(KSAM,KSAM,M); 
	n_mol = In[0]->MolList.size(); 
	success=Lat[0]->PrepareForCalculations(); 
	int n_mon=In[0]->MonList.size();
	for (int i=0; i<n_mon; i++) {success=Seg[i]->PrepareForCalculations(KSAM);}
	for (int i=0; i<n_mol; i++) {success=Mol[i]->PrepareForCalculations();}
	return success; 
}


bool System::ComputePhis(){

	bool success=true;

	PrepareForCalculations();
	float totphibulk=0;		
	float norm=0; 
	Zero(phitot,M); 
	int length=FrozenList.size();
	for (int i=0; i<length; i++) {
		float *phi_frozen=Seg[FrozenList[i]]->GetPhi();  
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
			float *phi=Mol[i]->phi+k*M;
			float *G1=Seg[Mol[i]->MolMonList[k]]->G1;
			Div(phi,G1,M); if (norm>0) Norm(phi,norm,M);
			k++;
		}

	}
	Mol[solvent]->phibulk=1.0-totphibulk; 
	norm=Mol[solvent]->phibulk/Mol[solvent]->chainlength;
	int k=0;
	length = Mol[solvent]->MolMonList.size();
	while (k<length) {
		float *phi=Mol[solvent]->phi+k*M;
		float *G1=Seg[Mol[solvent]->MolMonList[k]]->G1;
		Div(phi,G1,M); if (norm>0) Norm(phi,norm,M);
		k++;
	}
	for (int i=0; i<n_mol; i++) {
		int length=Mol[i]->MolMonList.size();
		int k=0;
		while (k<length) {
			k++; 
		}
	}
	for (int i=0; i<n_mol; i++) {
		int length=Mol[i]->MolMonList.size();
		int k=0;
		while (k<length) {
			float* phi_mon=Seg[Mol[i]->MolMonList[k]]->phi;
			float* phi_molmon = Mol[i]->phi + k*M; 
			Add(phi_mon,phi_molmon,M);
			Add(phitot,phi_mon,M); 
			Seg[i]->phibulk +=Mol[i]->fraction(Mol[i]->MolMonList[k])*Mol[i]->phibulk; 
			k++; 
		}
	}


//Testing time;
	
//normalize.
//make sure that the tagged positions are not in frozen range.  
	
	return success;  
}


 
