
class Output {
public:
	Output(vector<Input*>,string,int,int);

~Output();

	string name;
	vector<Input*> In; 
	int n_output; 
	int output_nr; 
	std::string template_;
	std::string type;
	bool write_bounds;
	bool append; 
		

	std::vector<string> OUT_name;
	std::vector<string> OUT_property;
	std::vector<string> OUT_value;

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckOutInput(void);
	void PutParameter(string); 
	string GetValue(string);
	bool Load(string); 


};
Output::Output(vector<Input*> In_,string name_,int outnr,int N_out) {
	In=In_; name=name_; n_output=N_out; output_nr=outnr; 
	KEYS.push_back("template"); 
	KEYS.push_back("type");
	KEYS.push_back("write_bounds");
	KEYS.push_back("append"); 

}
Output::~Output() {
}
void Output::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

bool Output::Load(string template_) {
return In[0]->LoadItems(template_, OUT_name, OUT_property, OUT_value);
}

bool Output::CheckOutInput() {
	return In[0]->CheckParameters("output",name,KEYS,PARAMETERS,VALUES);
}

string Output::GetValue(string parameter) {
	int length = PARAMETERS.size(); 
	int i=0;
	while (i<length) {
		if (PARAMETERS[i]==parameter) {return VALUES[i];}		
		i++;
	} 
	return " "; 
}

 
