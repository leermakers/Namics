#ifndef OUTPUTxH
#define OUTPUTxH
#include "namics.h"
#include "input.h"
#include "lattice.h"
#include "molecule.h"
#include "system.h"
class Output {
public:
	Output(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,vector<System*>,string,int,int);

~Output();

	string name;
	vector<Input*> In; 
	vector<Lattice*> Lat;
	vector<Segment*> Seg; 
	vector<Molecule*> Mol;
	vector<System*> Sys;
	int n_output; 
	int output_nr; 
	std::string template_;
	std::string type;
	bool write_bounds;
	bool append; 
	bool input_error; 
		

	std::vector<string> OUT_name;
	std::vector<string> OUT_property;
	std::vector<string> OUT_value;

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckOutInput(void);
	void PutParameter(string); 
	string GetValue(string);
	bool Load(string); 
	void vtk(string, double *);
	void density();
	void printlist();


};
#endif