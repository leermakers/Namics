#ifndef SEGMENTxH
#define SEGMENTxH
#include "namics.h"
#include "input.h"
#include "lattice.h"
#include "tools.h"
class Segment {
public:
	Segment(vector<Input*>,vector<Lattice*>,string,int,int);

~Segment(); 

	string name; 
	vector<Input*> In;  
	vector<Lattice*> Lat; 
	int n_seg; 
	int seg_nr;  
	double epsilon;
	double valence; 
	double phibulk; 
	string freedom;

	string filename; 
	bool block;
	int n_pos;  
	int MX,MY,MZ;
	int JX,JY;
	int BX1,BXM,BY1,BYM,BZ1,BZM;
	int M; 
	int xl,yl,zl,xh,yh,zh;

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
	double* GetPointer(string); 
	int GetValue(string,int&,double&,string&);


	int* H_Px; 
	int* H_Py;
	int* H_Pz;
	double* H_MASK;
	double* H_u;
	double* H_phi;
	int* Px; 
	int* Py;
	int* Pz;
	double* MASK; 
	double* G1;
	double* phi;
	double* phi_side;
	double* u;  

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
  	std::vector<string> VALUES; 
	bool CheckInput(); 
	void PutChiKEY(string); 
	string GetValue(string); 
	string GetFreedom();
	bool IsFree();
	bool IsPinned();
	bool IsFrozen();
	bool IsTagged(); 
	bool CreateMASK();
	double* GetMASK(); 
	double* GetPhi(); 
	void AllocateMemory();
	bool PrepareForCalculations(double*);  

};

#endif
