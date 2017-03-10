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
	double offset_first_layer; 
	double bond_length; 
	double *L;
	double lambda; 
	double *lambda0;
	double *lambda_1; 
	double *lambda1;
	//if you add new properties to this set, you should set the defaults or read value from input; got to CheckInput(). If the quantity has to go to output, also add it to PushOutput(). 

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;

	vector<string> BC;
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
	int GetValue(string,int&,double&,string&);
	double GetValue(double*,string);
	double WeightedSum(double*); 
	void TimesL(double*);
	void DivL(double*);
	void vtk(string, double*,string);
	void PutProfiles(FILE*,vector<double*>);

	bool CheckInput(int);
	bool PutM(void);
	void PutParameter(string); 
	string GetValue(string); 
	void DeAllocateMemory(void); 
	void AllocateMemory(void); 
	bool PrepareForCalculations(void); 
	void propagate(double*,double*, int, int);
	void remove_bounds(double*);
	void remove_bounds(int*);
	void set_bounds(double*);
	void set_bounds(int*); 
	void Side(double *, double *, int);
	bool ReadRange(int*, int*, int&, bool&, string, string, string);
	bool ReadRangeFile(string,int* H_p,int&, string, string);
	bool CreateMASK(int*, int*, int*, int, bool);
	bool GenerateGuess(double*, string, string, double, double);
	bool GuessVar(double*, double, string, double, double);
	bool ReadGuess(string,double*);
	bool StoreGuess(string,double*);
};

#endif
