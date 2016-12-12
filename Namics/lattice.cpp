class Lattice {
public:
	Lattice(vector<Input*>,string);

~Lattice();

	string name;
	vector<Input*> In;
	int MX,MY,MZ;		
	string LBX,LBY,LBZ,UBX,UBY,UBZ;
	int Volume;
	string lattice_type;
	float bond_length; 
	//if you add to this set of properties, you should set the defaults or read frominput as well. got to CheckInput(); 

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput();
	void PutParameter(string); 
	string GetValue(string); 

};
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
		options.push_back("shifted_mirror");
		Value=GetValue("lowerbound_x"); if (Value.length()>0) {
			In[0]->Get_string(Value,LBX,options,"for 'lowerbound_x' boundary condition not recognized."); 
		} else {
			LBX = "periodic";
		} 	
		Value=GetValue("lowerbound_y"); if (Value.length()>0) {
			In[0]->Get_string(Value,LBY,options,"for 'lowerbound_y' boundary condition not recognized."); 
		} else {
			LBY = "periodic";
		} 	
		Value=GetValue("lowerbound_z"); if (Value.length()>0) {
			In[0]->Get_string(Value,LBZ,options,"for 'lowerbound_z' boundary condition not recognized."); 
		} else {
			LBZ = "periodic";
		} 	
		Value=GetValue("upperbound_x"); if (Value.length()>0) {
			In[0]->Get_string(Value,UBX,options,"for 'upperbound_x' boundary condition not recognized."); 
		} else {
			UBX = "periodic";
		} 	
		Value=GetValue("upperbound_y"); if (Value.length()>0) {
			In[0]->Get_string(Value,UBY,options,"for 'upperbound_y' boundary condition not recognized."); 
		} else {
			UBY = "periodic";
		} 	
		Value=GetValue("uppperbound_z"); if (Value.length()>0) {
			In[0]->Get_string(Value,UBZ,options,"for 'upperbound_z' boundary condition not recognized."); 
		} else {
			UBZ = "periodic";
		} 	
		if (LBX != UBX) {cout <<"In X-direction the boundary conditions do not match!" << endl; success=false;}
		if (LBY != UBY) {cout <<"In Y-direction the boundary conditions do not match!" << endl; success=false;}
		if (LBZ != UBZ) {cout <<"In Z-direction the boundary conditions do not match!" << endl; success=false;}
			
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
 
