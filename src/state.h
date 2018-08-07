#ifndef STATExH
#define STATExH
#include "namics.h"
#include "input.h"
#include "segment.h"
class State {
public:
	State(vector<Input*>,vector<Segment*>,string);

~State();
	void DeAllocateMemory();
	void AllocateMemory(int,int);
	vector<Input*> In;	
	vector<Segment*> Seg;
	string name; 
 
	Real valence;
	string mon_name;
	int mon_nr; 
	int state_nr;
	Real alphabulk;
	bool fixed;
	bool in_reaction;

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
