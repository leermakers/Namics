#ifndef ENGINExH
#define ENGINExH
#include "namics.h"
#include "input.h"
#include "system.h"
#include "newton.h"
class Engine {
public:
	Engine(vector<Input*>,vector<System*>,vector<Newton*>,string);

~Engine();
	void AllocateMemory(); 
	string name; 
	vector<Input*> In; 
	vector<System*> Sys; 	
	vector<Newton*> New; 	
	string brand; 

	vector<string> ints;
	vector<string> doubles;
	vector<string> bools;
	vector<string> strings;
	vector<double> doubles_value;
	vector<int> ints_value;
	vector<bool> bools_value;
	vector<string> strings_value;
	void push(string,double);
	void push(string,int);
	void push(string,bool);
	void push(string,string);
	void PushOutput();
	double* GetPointer(string);
	int GetValue(string,int&,double&,string&);

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput();
	void PutParameter(string); 
	string GetValue(string); 
	bool Doit(void); 
};
#endif
