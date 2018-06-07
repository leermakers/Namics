#ifndef ALIASxH
#define ALIASxH
#include "namics.h"
#include "input.h"
#include "lattice.h"
#include "tools.h"
class Alias {
public:
	Alias(vector<Input*>,vector<Lattice*>,string);

~Alias();
	void DeAllocateMemory();
	void AllocateMemory(int,int);
	vector<Lattice*> Lat; 
	int value;
	string composition;
	bool active;
	vector<int> frag;
	Real* H_phi;
	Real* phi;
	Real* rho;  
	bool clamp;
	
	string name; 
	vector<Input*> In; 

	vector<string> ints;
	vector<string> Reals;
	vector<string> bools;
	vector<string> strings;
	vector<Real> Reals_value;
	vector<int> ints_value;
	vector<bool> bools_value;
	vector<string> strings_value;
	void push(string,Real);
	void push(string,int);
	void push(string,bool);
	void push(string,string);
	void PushOutput();
	Real* GetPointer(string);
	int GetValue(string,int&,Real&,string&);	
	void PrepareForCalculations();

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string); 
	string GetValue(string); 
};
#endif
