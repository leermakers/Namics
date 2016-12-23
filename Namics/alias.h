#ifndef ALIASxH
#define ALIASxH
#include "input.h"
class Alias {
public:
	Alias(vector<Input*>,string);

~Alias();
	string value;
	
	string name; 
	vector<Input*> In; 	

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput();
	void PutParameter(string); 
	string GetValue(string); 
};
#endif
