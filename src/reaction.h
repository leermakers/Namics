#ifndef REACTIONxH
#define REACTIONxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "state.h"
class Reaction {
public:
	Reaction(vector<Input*>,vector<Segment*>,vector<State*>,string);

~Reaction();
	void DeAllocateMemory();
	void AllocateMemory(int,int);
	vector<State*> Sta; 
	vector<Segment*> Seg;
	vector<int> Sto;
	vector<int> State_nr;
	vector<int> State_in_seg_nr; 
	vector<int> Seg_nr; 
	Real K;
	Real pK;
	string equation;
	
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
	Real* GetPointer(string,int&);
	int* GetPointerInt(string,int&);
	int GetValue(string,int&,Real&,string&);	
	void PrepareForCalculations();

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string); 
	string GetValue(string); 
	Real ChemIntBulk(State*);
	Real pKeff();
	Real Residual_value();
};
#endif
