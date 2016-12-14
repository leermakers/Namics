class Lattice {
public: 
	Lattice(vector<Input*>,string);

~Lattice();

	string name;
	vector<Input*> In;
	int MX,MY,MZ;		 
	int BX1,BY1,BZ1,BXM,BYM,BZM;
	int JX,JY,M;
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
	void AllocateMemory(void); 
	bool PrepareForCalculations(void); 
	void propagate(float*,float*, int, int);
	void remove_bounds(float*);
	void set_bounds(float*); 
	void Side(float *, float *, int);

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
		JX=(MX+2)*(MY+2); JY=(MY+2); M = (MX+2)*(MY+2)*(MZ+2);   
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
		//options.push_back("shifted_mirror");
		string VALUE1,VALUE2,VALUE3,VALUE4,VALUE5,VALUE6;
		Value.clear();
		Value=GetValue("lowerbound_x"); if (Value.length()>0) {
			In[0]->Get_string(Value,VALUE1,options,"for 'lowerbound_x' boundary condition not recognized.");
			if (VALUE1=="mirror_1") BX1=1; 
			if (VALUE1=="mirror_2") BX1=2;
			if (VALUE1=="periodic") BX1=MX;
		} else {
			VALUE1="periodic";
			BX1=MX;
		} 
		Value.clear();	
		Value=GetValue("lowerbound_y"); if (Value.length()>0) {
			In[0]->Get_string(Value,VALUE2,options,"for 'lowerbound_y' boundary condition not recognized.");
			if (VALUE2=="mirror_1") BY1=1; 
			if (VALUE2=="mirror_2") BY1=2;
			if (VALUE2=="periodic") BY1=MY; 
		} else {
			VALUE2="periodic";
			BY1=MY;
		} 
		Value.clear();	
		Value=GetValue("lowerbound_z"); if (Value.length()>0) {
			In[0]->Get_string(Value,VALUE3,options,"for 'lowerbound_z' boundary condition not recognized.");
			if (VALUE3=="mirror_1") BZ1=1; 
			if (VALUE3=="mirror_2") BZ1=2;
			if (VALUE3=="periodic") BZ1=MZ;  
		} else {
			VALUE3="periodic";
			BZ1=MZ;
		} 
		Value.clear();	
		Value=GetValue("upperbound_x"); if (Value.length()>0) {
			In[0]->Get_string(Value,VALUE4,options,"for 'upperbound_x' boundary condition not recognized."); 
			if (VALUE4=="mirror_1") BXM=MX-1; 
			if (VALUE4=="mirror_2") BXM=MX-2;
			if (VALUE4=="periodic") BXM=1;
		} else {
			VALUE4="periodic";
			BXM=1;
		} 
		Value.clear();	
		Value=GetValue("upperbound_y"); if (Value.length()>0) {
			In[0]->Get_string(Value,VALUE5,options,"for 'upperbound_y' boundary condition not recognized."); 
			if (VALUE5=="mirror_1") BYM=MY-1; 
			if (VALUE5=="mirror_2") BYM=MY-2;
			if (VALUE5=="periodic") BYM=1;
		} else {
			VALUE5="periodic";
			BYM=1;
		} 
		Value.clear();	
		Value=GetValue("upperbound_z"); if (Value.length()>0) { 
			In[0]->Get_string(Value,VALUE6,options,"for 'upperbound_z' boundary condition not recognized."); 
			if (VALUE6=="mirror_1") BZM=MZ-1; 
			if (VALUE6=="mirror_2") BZM=MZ-2;
			if (VALUE6=="periodic") BZM=1;
		} else {
			VALUE6="periodic";
			BZM=1;
		} 	
		if (VALUE1 != VALUE4) {cout <<"In x-direction the boundary conditions do not match:" + VALUE1 << " and " <<  VALUE4 << endl; success=false;}
		if (VALUE2 != VALUE5) {cout <<"In y-direction the boundary conditions do not match:" + VALUE2 << " and " <<  VALUE5 << endl; success=false;}
		if (VALUE3 != VALUE6) {cout <<"In z-direction the boundary conditions do not match:" + VALUE3 << " and " <<  VALUE6 << endl; success=false;}
			
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

void Lattice::AllocateMemory(void) {
	
}
bool Lattice::PrepareForCalculations(void) {
	bool success=true;
	return success; 
}

void Lattice::Side(float *X_side, float *X, int M) {
	Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Norm(X_side,1.0/6.0,M);
}

void Lattice::propagate(float *G, float *G1, int s_from, int s_to) { 
	float *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	Zero(gs,M);
	SetBoundaries(gs_1,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(gs+JX,gs_1,M-JX); Add(gs,gs_1+JX,M-JX);
	Add(gs+JY,gs_1,M-JY); Add(gs,gs_1+JY,M-JY);
	Add(gs+1,gs_1,M-1);  Add(gs,gs_1+1, M-1);
	Norm(gs,1.0/6.0,M); Times(gs,gs,G1,M);
}

void Lattice::remove_bounds(float *X){ 
	RemoveBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
}
 
void Lattice::set_bounds(float *X){  
	SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
}
 