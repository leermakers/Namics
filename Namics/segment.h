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
	int* H_Px; 
	int* H_Py;
	int* H_Pz;
	double* H_MASK;
	double* H_G1;
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
