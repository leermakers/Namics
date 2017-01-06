#ifndef ALIASxH
#define ALIASxH
#include "input.h"
class Alias {
public:
	Alias(vector<Input*>,string);

~Alias();
	void AllocateMemory();
	int value;
	string composition;
	
	string name; 
	vector<Input*> In; 

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
	bool CheckInput(int);
	void PutParameter(string); 
	string GetValue(string); 
};
#endif