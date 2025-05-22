#ifndef MICROEMULSIONxH
#define MICROEMULSIONxH
#include "namics.h"
#include "segment.h"
#include "state.h"
#include "reaction.h"
#include "molecule.h"
#include "lattice.h"
#include "system.h"
#include "solve_scf.h"
#include "output.h"
#include "variate.h"


class Microemulsion {
public:
	Microemulsion(vector<Input*>,vector<Output*>,vector<Lattice*>,vector<Segment*>,vector<State*>,vector<Reaction*>,vector<Molecule*>,vector<System*>,vector<Solve_scf*>,vector<Variate*>,string);

~Microemulsion();
  	const string name;
  	const vector<Input*> In;
  	vector<Output*> Out;
  	const vector<Lattice*> Lat;
  	const vector<Molecule*> Mol;
  	const vector<Segment*> Seg;
  	const vector<State*> Sta;
	const vector<Reaction*> Rea;
  	const vector<System*> Sys;
 	const vector<Solve_scf*> New;
 	const vector<Variate*> Var;
	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	int surfactant;
	int co_solvent;
	int start;
	Real * X;
	string METHOD;
	vector<string> MONLIST;
	vector<string> STATELIST;
	bool CHARGED;
	int MX;
	int MY;
	int MZ;
	int fjc_old;
	int search_nr;
	int ets_nr;
	int etm_nr;
	int target_nr;
	int bm_nr;
	int subloop;

	bool CheckInput(int);
	bool Doit(Real*,string,vector<string>,vector<string>,bool,int,int,int,int,int,int,int,int,int,int);
	string GetValue(string);
	bool FixedPoint();
	bool WriteResults();
};
#endif
