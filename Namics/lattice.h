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
	int Volume;
	string lattice_type;
	double bond_length; 
	//if you add to this set of properties, you should set the defaults or read frominput as well. got to CheckInput(); 

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
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