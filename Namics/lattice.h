#ifndef LATTICExH
#define LATTICExH

#include "namics.h"
#include "input.h"
#include "tools.h"

class Lattice {
public: 
	Lattice(vector<Input*>,string);

~Lattice();

	string name;
	vector<Input*> In;
	int MX,MY,MZ;	
	int Mx,My,Mz;	 
	int BX1,BY1,BZ1,BXM,BYM,BZM;
	int JX,JY,M;
	int volume;
	string lattice_type;
	double bond_length; 
	//if you add new properties to this set, you should set the defaults or read value from input; got to CheckInput(). If the quantity has to go to output, also add it to PushOutput(). 

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;

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

	bool CheckInput();
	void PutParameter(string); 
	string GetValue(string); 
	void AllocateMemory(void); 
	bool PrepareForCalculations(void); 
	void propagate(double*,double*, int, int);
	void remove_bounds(double*);
	void set_bounds(double*); 
	void Side(double *, double *, int);
};

#endif
