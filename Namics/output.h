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
#include "engine.h"
class Output {
public:
	Output(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,vector<System*>,vector<Newton*>,vector<Engine*>,string,int,int);

~Output();

	string name;
	vector<Input*> In; 
	vector<Lattice*> Lat;
	vector<Segment*> Seg; 
	vector<Molecule*> Mol;
	vector<System*> Sys;
	vector<Newton*> New;
	vector<Engine*> Eng;
	int n_output; 
	int output_nr; 
	bool write_bounds;
	bool append; 
	bool input_error; 
		

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
	void vtk(string, rene *);
	void density();
	void printlist();
	void WriteOutput(int);
	int GetValue(string, string, string, int&, rene&, string&);
	rene* GetPointer(string, string, string);

};
#endif
