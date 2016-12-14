
class Engine {
public:
	Engine(vector<Input*>,vector<System*>,string);

~Engine();

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
Engine::Engine(vector<Input*> In_,vector<System*> Sys_, string name_) {
	In=In_; name=name_;   Sys=Sys_; 
	KEYS.push_back("MC"); KEYS.push_back("SN_MD");
}
Engine::~Engine() {
}
void Engine::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

bool Engine::CheckInput() {
	return In[0]->CheckParameters("engine",name,KEYS,PARAMETERS,VALUES);
}
 
string Engine::GetValue(string parameter){
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
 
