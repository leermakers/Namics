#ifndef VARIATExH
#define VARIATExH
#include "tools.h"
#include "input.h"
#include "lattice.h"
#include "segment.h"
#include "molecule.h"
#include "system.h"
class Variate {
public:
	Variate(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,vector<System*>,string);

~Variate();
	void DeAllocateMemory();
	void AllocateMemory();
	vector<Input*> In; 
	vector<Lattice*> Lat; 
	vector<Segment*> Seg;
	vector<Molecule*> Mol;
	vector<System*> Sys;

	string name; 
	int search_nr;
	int scan_nr;
	int target_nr;
	int scanning;
	int searching;
	int targeting;
	int steps;
	Real step;
	Real end_value;
	string scale;
	int num_of_cals;

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
	bool PutVarScan(int);
	bool ResetScanValue();
	void PutValue(Real);
	Real GetValue(void);
	Real GetError(void);

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
