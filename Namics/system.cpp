 
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
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput();
	void PutParameter(string); 
	string GetValue(string); 
	bool CheckChi_values(int);

};
System::System(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,string name_) {
	Seg=Seg_; Mol=Mol_; Lat=Lat_; In=In_; name=name_; 
	KEYS.push_back("calculation_type");

	
}
System::~System() {
}
void System::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

bool System::CheckInput() {
	return In[0]->CheckParameters("sys",name,KEYS,PARAMETERS,VALUES);
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
	return " " ; 
}

bool System::CheckChi_values(int n_seg){
	bool success=true;
	CHI = new double[n_seg*n_seg]; for (int i=0; i<n_seg; i++) for(int k=0; k<n_seg; k++) CHI[i*n_seg+k]=123;
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) {
		string s_value;
		double chi_value;
		s_value = Seg[i]->GetValue("chi-"+Seg[k]->name); if (In[0]->Get_double(s_value,chi_value)) CHI[i*n_seg+k]=chi_value;
	}
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) if (CHI[i*n_seg+k] == 123) CHI[i*n_seg+k] = CHI[k*n_seg+i];
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) if (CHI[i*n_seg+k] == 123) CHI[i*n_seg+k] = 0;
	for (int i=0; i<n_seg; i++) for (int k=i+1; k<n_seg; k++) if (CHI[i*n_seg+k] != CHI[k*n_seg+i]) {
		cout <<"CHI-value symmetry violated: chi("<<Seg[i]->name<<","<<Seg[k]->name<<") is not equal to chi("<<Seg[k]->name<<","<<Seg[i]->name<<")"<<endl; success=false;
	}
	return success; 
}


 
