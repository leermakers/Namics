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
 
