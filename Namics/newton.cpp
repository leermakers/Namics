class Newton {
public:
	Newton(vector<Input*>,vector<System*>,string);

~Newton();

	string name;
	vector<Input*> In; 
	vector<System*> Sys; 
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput();
	void PutParameter(string); 
	string GetValue(string); 

};
Newton::Newton(vector<Input*> In_,vector<System*> Sys_,string name_) {
	In=In_; name=name_; Sys=Sys_;   
	KEYS.push_back("method"); KEYS.push_back("n_iterations"); KEYS.push_back("e_info"); KEYS.push_back("s_info");
	KEYS.push_back("delta_max"); KEYS.push_back("m"); KEYS.push_back("i_info");
}
Newton::~Newton() {
}
void Newton::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

bool Newton::CheckInput() {
	return In[0]->CheckParameters("newton",name,KEYS,PARAMETERS,VALUES);
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
	return " " ; 
}
 
