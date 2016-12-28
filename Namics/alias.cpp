#include "alias.h"

Alias::Alias(vector<Input*> In_, string name_) {
	In=In_; name=name_;   
	KEYS.push_back("value"); 
}
Alias::~Alias() {
}
void Alias::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

bool Alias::CheckInput() {
	bool success=true;
	success= In[0]->CheckParameters("alias",name,KEYS,PARAMETERS,VALUES);
	if (success) {
		if (GetValue("value").size()>0) {
			In[0]->Get_string(GetValue("value"),value);
		}
	}
	return success;
}
 
string Alias::GetValue(string parameter){
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
 
void Alias::push(string s, double X) {
	doubles.push_back(s);
	doubles_value.push_back(X); 
}
void Alias::push(string s, int X) {
	ints.push_back(s);
	ints_value.push_back(X); 
}
void Alias::push(string s, bool X) {
	bools.push_back(s);
	bools_value.push_back(X); 
}
void Alias::push(string s, string X) {
	strings.push_back(s);
	strings_value.push_back(X); 	
}
void Alias::PushOutput() {
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	doubles.clear();
	doubles_value.clear();
	ints.clear();
	ints_value.clear();  
}