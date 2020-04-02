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
	vector<int> mx;
	vector<int> my;
	vector<int> mz;
	vector<int> m;
	vector<int> jx;
	vector<int> jy;
	vector<int> n_box;
	vector<string> BC;
	int BX1,BY1,BZ1,BXM,BYM,BZM;
	int *B_X1;
	int *B_Y1;
	int *B_Z1;
	int *B_XM;
	int *B_YM;
	int *B_ZM;
	int JX,JY,JZ,M;
	bool all_lattice;
	bool ignore_sites;
	bool stencil_full;
	int sub_box_on;
	Real volume;
	Real Accesible_volume;

	int gradients;
	string lattice_type;
	string geometry;
	Real offset_first_layer;
	Real bond_length;
	Real *L;
	int Z;
	Real lambda;
	Real *lambda0;
	Real *lambda_1;
	Real *lambda1;
	Real *LAMBDA;
	int fjc, FJC;
	Real *X;
	int VarInitValue;
	string Var_type;
	int Var_target;
	int Var_step;
	int Var_end_value;
	//if you add new properties to this set, you should set the defaults or read value from input; got to CheckInput(). If the quantity has to go to output, also add it to PushOutput().

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;

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
	int P(int,int,int);
	int P(int,int);
	int P(int);

	int GetValue(string,int&,Real&,string&);
	Real GetValue(Real*,string);
	Real WeightedSum(Real*);
	Real Moment(Real*,int);
	void TimesL(Real*);
	void DivL(Real*);
	void vtk(string, Real*,string,bool);
	void PutProfiles(FILE*,vector<Real*>,bool);

	bool CheckInput(int);
	bool PutM(void);
	bool PutSub_box(int,int,int,int);
	void PutParameter(string);
	string GetValue(string);
	void DeAllocateMemory(void);
	void AllocateMemory(void);
	bool PrepareForCalculations(void);
	void propagate(Real*,Real*, int, int,int);

	template <typename T>
	void remove_bounds(T*);

	template <typename T>
	void set_bounds(T*);

	void Side(Real *, Real *, int);
	bool ReadRange(int*, int*, int&, bool&, string, string, string);
	bool ReadRangeFile(string,int* H_p,int&, string, string);
	bool FillMask(int*, vector<int>, vector<int>, vector<int>, string);
	bool CreateMASK(int*, int*, int*, int, bool);
	bool GenerateGuess(Real*, string, string, Real, Real);
	bool GuessVar(Real*, Real, string, Real, Real);
	void DistributeG1(Real*, Real*, int*, int*, int*, int);
	void CollectPhi(Real*, Real*, Real*, int*, int*, int*, int);
	void ComputeGN(Real*, Real*, int*, int*, int*, int*, int*, int*, int, int);
	Real ComputeTheta(Real*);
	void UpdateEE(Real*, Real*,Real*);
	void UpdatePsi(Real*, Real*, Real* , Real*, int*,bool,bool);
	void UpdateQ(Real*,Real*,Real*,Real*,int*,bool);
	bool ReadGuess(string, Real* ,string&, vector<string>&, vector<string>&, bool&, int&, int&, int&, int&, int);
	bool StoreGuess(string,Real*,string, vector<string>,vector<string>, bool,int);
	bool PutVarInfo(string,string,Real);
	bool UpdateVarInfo(int);
	bool ResetInitValue();
	int PutVarScan(int, int);
};
#endif
