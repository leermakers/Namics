#ifndef OUTPUTxH
#define OUTPUTxH
#include "namics.h"
#include "input.h"
#include "lattice.h"
#include "segment.h"
#include "molecule.h"
#include "system.h"
#include "newton.h"
#include "alias.h"
#include <limits.h>
#include <unistd.h>

class Output {
public:
	Output(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,vector<System*>,vector<Newton*>,string,int,int);

~Output();

	string name;
	vector<Input*> In;
	vector<Lattice*> Lat;
	vector<Segment*> Seg;
	vector<Molecule*> Mol;
	vector<System*> Sys;
	vector<Newton*> New;
	int n_output;
	int start;
	int output_nr;
	bool write_bounds;
	bool append;
	bool input_error;
	string output_folder;
	string bin_folder;
	bool use_output_folder;

  	vector<string> ints;
 	vector<string> Reals;
 	vector<string> bools;
 	vector<string> strings;
  	vector<Real> Reals_value;
	vector<int > SizeVectorReal;
	vector<int > SizeVectorInt;
	vector<Real*> PointerVectorReal;
	vector<int*> PointerVectorInt;
	vector<int> pointer_size;
  	vector<int> ints_value;
  	vector<bool> bools_value;
  	vector<string> strings_value;

	std::vector<string> OUT_key;
	std::vector<string> OUT_name;
	std::vector<string> OUT_prop;

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string);
	string GetValue(string);
	bool Load();
	void vtk(string, Real *);
	void density();
	void printlist();
	void WriteOutput(int);
	int GetValue(string, string, string, int&, Real&, string&);
	Real* GetPointer(string, string, string, int&);
	int* GetPointerInt(string, string, string, int&);
	int GetValue(string, string, int& , Real& , string&);
	void push(string, Real);
	void push(string, int);
	void push(string, bool);
	void push(string, string);

};
#endif
