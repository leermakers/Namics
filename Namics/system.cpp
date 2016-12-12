 
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
	vector<string> SysMonList; 
	float* phitot; 
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
	bool CheckChi_values(int);
	void PrepareForCalculations();

};
System::System(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,string name_) {
	Seg=Seg_; Mol=Mol_; Lat=Lat_; In=In_; name=name_; 
	KEYS.push_back("calculation_type");
	KEYS.push_back("GPU"); 	
}
System::~System() {
}

bool System::CheckInput() {
	bool success=true;
	bool solvent_found=false; solvent = -1; tag_segment=-1;  //value -1 means no solvent defined. tag_segment=-1; 
	float phibulktot=0;  
	success= In[0]->CheckParameters("sys",name,KEYS,PARAMETERS,VALUES);
	if (success) {
		int i=0;
		int length = In[0]->MolList.size(); 
		while (i<length) {
			int j=0; 
			int LENGTH=Mol[i]->MolMonList.size(); 
			while (j<LENGTH) {
				if (!In[0]->InSet(SysMonList,Mol[i]->MolMonList[j])) SysMonList.push_back(Mol[i]->MolMonList[j]);
				j++;
			}
			if (Mol[i]->GetFreedom()=="free") phibulktot +=Mol[i]->phibulk;
			if (Mol[i]->GetFreedom()=="solvent") {solvent_found=true; solvent=i; }
			if (Mol[i]->IsTagged()) { 
				tag_segment = Mol[i]->tag_segment; 
				Mol[i]->n=1.0*Seg[tag_segment]->n_pos;
				Mol[i]->theta= Mol[i]->n*Mol[i]->chainlength;  
			}
			i++; 
		}
		if (!solvent_found && phibulktot >0.99999999 && phibulktot < 1.0000000001) {
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
	return success; 
}

void System::PrepareForCalculations() {

	MX=Lat[0]->MX; 
	MY=Lat[0]->MY;
	MZ=Lat[0]->MZ;
	M=(MX+2)*(MY+2)*(MZ+2); 	
//define on CPU

#ifdef CUDA
//define on GPU
	phitot = (float*)AllOnDev(M); 	

#else
	phitot = new float[M]; 
#endif

}


 
