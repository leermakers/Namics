#ifndef ENGINExH
#define ENGINExH
#include "namics.h"
#include "input.h"
#include "system.h"
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
#endif