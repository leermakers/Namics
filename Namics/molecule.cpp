
class Molecule { 
public:
	Molecule(vector<Input*>,vector<Segment*>,string,int,int);

~Molecule();

	string name; 
	vector<Input*> In;
	//vector<Lattice*> Lat; 
	vector<Segment*> Seg; 
	int n_mol; 
	int mol_nr; 
	double theta;
	double phibulk; 
	string freedom;
	double n; 
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(void);
	void PutParameter(string);
	bool DeComposition(void); 
	int ChainLength(void); 
	double Theta(void);
	string GetValue(string); 

};
Molecule::Molecule(vector<Input*> In_,vector<Segment*> Seg_, string name_,int molnr,int N_mol) {
	In=In_; Seg=Seg_; name=name_; n_mol=N_mol; mol_nr=molnr; 
	KEYS.push_back("freedom"); 
	KEYS.push_back("chainlength");
	KEYS.push_back("theta");
	KEYS.push_back("phibulk");
	KEYS.push_back("n");
}
Molecule::~Molecule() {
}
void Molecule::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

bool Molecule::CheckInput() {
	return In[0]->CheckParameters("mol",name,KEYS,PARAMETERS,VALUES);
}

string Molecule::GetValue(string parameter) {
	int length = PARAMETERS.size(); 
	int i=0;
	while (i<length) {
		if (PARAMETERS[i]==parameter) { return VALUES[i];}		
		i++;
	} 
	return " "; 
}
 
