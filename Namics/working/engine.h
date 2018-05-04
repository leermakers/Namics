#ifndef ENGINExH
#define ENGINExH
#include "namics.h"
#include "input.h"
#include "system.h"
#include "newton.h"
#include "mesodyn.h"
#include "monte.h"
class Engine {
public:
	Engine(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,vector<System*>,vector<Newton*>,string);

~Engine();
	void AllocateMemory(); 
	string name; 
	vector<Input*> In; 
	vector<Lattice*> Lat;
	vector<Segment*> Seg;
	vector<Molecule*> Mol;
	vector<System*> Sys; 	
	vector<Newton*> New;
	vector<Mesodyn*> Mdyn;  	
	vector<Monte*> Mont;
	string brand; 

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
	

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string); 
	string GetValue(string); 

};
#endif
