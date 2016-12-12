class Newton {
public:
	Newton(vector<Input*>,vector<System*>,string);

~Newton();

	string name;
	vector<Input*> In; 
	vector<System*> Sys; 
	int iterationlimit,m,i_info;
	float tolerance,delta_max;
	bool e_info,s_info; 
	string method;
	bool store_guess;
	bool read_guess;  
	string stop_criterion;
	int iv; 
	int M; 
	float* x;
	float* x0;
	float* g;
	float* xR;
	float* x_x0;
	float* alpha;
//if properties are added, also read them form input and/or set the default. See CheckInput() below. 
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput();
	void PutParameter(string); 
	string GetValue(string); 
	void PrepareForCalculations(void);

};
Newton::Newton(vector<Input*> In_,vector<System*> Sys_,string name_) {
	In=In_; name=name_; Sys=Sys_;   
	KEYS.push_back("method"); KEYS.push_back("e_info"); KEYS.push_back("s_info");
	KEYS.push_back("delta_max"); KEYS.push_back("m"); KEYS.push_back("i_info");
	KEYS.push_back("iterationlimit" ); KEYS.push_back("tolerance"); KEYS.push_back("store_guess"); KEYS.push_back("read_guess");
	KEYS.push_back("stop_criterion"); 
}
Newton::~Newton() {
}

bool Newton::CheckInput() {
	bool success=true;
	string value;
	success=In[0]->CheckParameters("newton",name,KEYS,PARAMETERS,VALUES);
	if (success) {
		e_info=In[0]->Get_bool(GetValue("e_info"),true); 
		s_info=In[0]->Get_bool(GetValue("s_info"),false); 
		delta_max=In[0]->Get_float(GetValue("delta_max"),0.1); 
		if (delta_max < 0 || delta_max>100) {delta_max = 0.1;  cout << "Value of delta_max out of range 0..100, and value set to default value 0.1" <<endl; } 
		tolerance=In[0]->Get_float(GetValue("tolerance"),1e-5);
		if (tolerance < 1e-12 ||tolerance>10) {tolerance = 1e-5;  cout << "Value of tolerance out of range 1e-12..10 Value set to default value 1e-5" <<endl; }  
		method=In[0]->Get_string(GetValue("method"),"DIIS");
		m=In[0]->Get_int(GetValue("m"),10); 
		if (m < 0 ||m>100) {m=10;  cout << "Value of 'm' out of range 0..100, value set to default value 10" <<endl; }
		store_guess=In[0]->Get_bool(GetValue("store_guess"),true);		
		store_guess=In[0]->Get_bool(GetValue("read_guess"),false);
		if (GetValue("stop_criterion").size() > 0) {
			vector<string>options;
			options.push_back("norm_of_g");
			options.push_back("max_of_element_of_|g|");
			if (!In[0]->Get_string(GetValue("stop_criterion"),stop_criterion,options,"In newton the stop_criterion setting was not recognised")) {success=false; };
		}
	}
	return success; 
}

void Newton::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
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
	return "" ; 
}

void Newton::PrepareForCalculations() {
	M=Sys[0]->M;
	iv=M*Sys[0]->SysMonList.size(); 
	
//define on CPU

#ifdef CUDA
//define on GPU
	x = (float*)AllOnDev(iv); u=x; Zero(x,iv);
	x0 = (float*)AllOnDev(iv);
	g= (float*)AllOnDev(iv);
	xR= (float*)AllOnDev(m*iv);
	x_x0= (float*)AllOnDev(m*iv);
	alpha= (float*)AllOnDev(M);
#else
	x = new float[iv]; Zero(x,iv);
	x0 = new float[iv];
	g = new float[iv];
	xR = new float[m*iv];
	x_x0 =new float[m*iv];
	alpha = new float[M];
#endif

}
 
