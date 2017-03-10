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
	int gradients; 
	string lattice_type;
	string geometry;
	rene offset_first_layer; 
	rene bond_length; 
	rene *L;
	rene lambda; 
	rene *lambda0;
	rene *lambda_1; 
	rene *lambda1;
	//if you add new properties to this set, you should set the defaults or read value from input; got to CheckInput(). If the quantity has to go to output, also add it to PushOutput(). 

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;

	vector<string> BC;
	vector<string> ints;
	vector<string> renes;
	vector<string> bools;
	vector<string> strings;
	vector<rene> renes_value;
	vector<int> ints_value;
	vector<bool> bools_value;
	vector<string> strings_value;
	void push(string,rene);
	void push(string,int);
	void push(string,bool);
	void push(string,string);
	void PushOutput();
	int GetValue(string,int&,rene&,string&);
	rene GetValue(rene*,string);
	rene WeightedSum(rene*); 
	void TimesL(rene*);
	void DivL(rene*);
	void vtk(string, rene*,string);
	void PutProfiles(FILE*,vector<rene*>);

	bool CheckInput(int);
	bool PutM(void);
	void PutParameter(string); 
	string GetValue(string); 
	void DeAllocateMemory(void); 
	void AllocateMemory(void); 
	bool PrepareForCalculations(void); 
	void propagate(rene*,rene*, int, int);
	void remove_bounds(rene*);
	void remove_bounds(int*);
	void set_bounds(rene*);
	void set_bounds(int*); 
	void Side(rene *, rene *, int);
	bool ReadRange(int*, int*, int&, bool&, string, string, string);
	bool ReadRangeFile(string,int* H_p,int&, string, string);
	bool CreateMASK(int*, int*, int*, int, bool);
	bool GenerateGuess(rene*, string, string, rene, rene);
	bool GuessVar(rene*, rene, string, rene, rene);
	bool ReadGuess(string,rene*);
	bool StoreGuess(string,rene*);
};

#endif
