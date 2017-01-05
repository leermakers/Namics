#ifndef NEWTONxH
#define NEWTONxH
#include "tools.h"
#include "namics.h"
#include "input.h"
#include "system.h"
#include "segment.h"
#include "lattice.h"
#include "molecule.h"
class Newton {
public:
	Newton(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,vector<System*>,string);

~Newton();

	string name; 
	vector<Input*> In; 
	vector<System*> Sys;
	vector<Segment*> Seg;
	vector<Lattice*> Lat;  
	vector<Molecule*> Mol;
	int iterationlimit,m,i_info;
	int k_diis,it; 
	double tolerance,delta_max;
	double residual;
	bool e_info,s_info; 
	string method;
	bool store_guess;
	bool read_guess;  
	string stop_criterion;
	int iv; 
	int M,MX,MY,MZ,JX,JY;

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

	double* x;
	double* x0;
	double* g;
	double* xR;
	double* x_x0;
	double* alpha;
	double* Aij;
	double* Ci;
	double* Apij; 
//if properties are added, also read them form input and/or set the default. See CheckInput() below. 
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string); 
	string GetValue(string); 
	bool Solve();
	void AllocateMemory(); 
	bool PrepareForCalculations(void);
	void Ax(double* , double* , int );
	void DIIS(double* , double* , double* , double*, double* ,double* , int , int , int );
	void ComputeG(); 
	void ComputePhis();
	void ComputeG_ext();
	bool Iterate_Picard();
	bool Iterate_DIIS();
	void Message(int, int,double, double); 
	bool PutU();

};
#endif
