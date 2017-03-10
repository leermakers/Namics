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
	vector<int> px1;
	vector<int> px2;
	vector<int> py1;
	vector<int> py2;
	vector<int> pz1;
	vector<int> pz2;
	vector<int> bx;
	vector<int> by;
	vector<int> bz; 
	int n_seg; 
	int seg_nr;  
	rene epsilon;
	rene valence; 
	rene phibulk; 
	string freedom;
	rene guess_u; 

	string filename; 
	bool block;
	int n_pos;  
	int n_box;
	int* r;
	int mx,my,mz,m;

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
	rene* GetPointer(string); 
	int GetValue(string,int&,rene&,string&);


	int* H_P; 
	int* H_MASK;
	rene* H_u;
	rene* H_phi;
	int* P; 
	int* MASK;
	rene* G1;
	rene* phi;
	rene* phi_side;
	rene* u;  

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
  	std::vector<string> VALUES; 
	bool CheckInput(int); 
	void PutChiKEY(string); 
	string GetValue(string); 
	string GetFreedom();
	bool GetClamp(string); 
	bool IsFree();
	bool IsPinned();
	bool IsFrozen();
	bool IsTagged();
	bool IsClamp(); 
	int* GetMASK(); 
	rene* GetPhi(); 
	void AllocateMemory();
	bool PrepareForCalculations(int*);  

};

#endif
