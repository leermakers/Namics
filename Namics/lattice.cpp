class Lattice {
public:
	Lattice(vector<Input*>,string);

~Lattice();

	string name;
	vector<Input*> In; 		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput();
	void PutParameter(string); 
	string GetValue(string); 
	int MX,MY,MZ;

};
Lattice::Lattice(vector<Input*> In_,string name_) {
	In=In_; name=name_; 
	KEYS.push_back("n_layers_x"); KEYS.push_back("n_layers_y"); KEYS.push_back("n_layers_z");
 	KEYS.push_back("bond_length");
}
Lattice::~Lattice() {
}
void Lattice::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

bool Lattice::CheckInput() {
	return In[0]->CheckParameters("lat",name,KEYS,PARAMETERS,VALUES);
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
	return " " ; 
}
 
