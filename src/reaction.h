#ifndef REACTIONxH
#define REACTIONxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "state.h"

#include <cmath>

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

	int Var_steps;
	Real Var_step;
	Real Var_end_value;
	Real Var_start_value;
	Real Var_start_search_value;
	int num_of_steps;
	int Var_target;
	int Var_scan_value;
	int Var_search_value;
	string Var_type;
	string scale;
	Real Var_target_value;
	
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
	bool GuessAlpha();
	bool PutAlpha(Real);

	bool PutVarInfo(string,string,Real);
	int PutVarScan(Real,Real,int,string);
	bool ResetInitValue();
	bool UpdateVarInfo(int);
	Real GetError();
	Real GetValue();
	void PutValue(Real);
};
#endif
