#ifndef BATExH
#define BATExH
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


class Balancedtensionless {
public:
	Balancedtensionless(vector<Input*>,vector<Output*>,vector<Lattice*>,vector<Segment*>,vector<State*>,vector<Reaction*>,vector<Molecule*>,vector<System*>,vector<Solve_scf*>,vector<Variate*>,string);

~Balancedtensionless();
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
 	CP ControlParameter;
	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	Real *phi_oil;
	int surfactant;
	int water;
	int monT;
	int monH;
	int monW;
	Real chi_start;
	Real chi_step;
	Real chi_end;
	int n_steps;
	Real GuessS;
	Real GuessC;
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
	bool kal_append;
	Real j_tolerance;
	int j_info;
	Real GS0;
	Real GC0;
	Real ini_oil, ini_water, ini_surf, ini_chi;
	bool previous_guess;
	Real caution_factor;
	bool HequalW;
	Real J0m;

	bool CheckInput(int);
	bool Doit(Real*,string,vector<string>,vector<string>,bool,int,int,int,int,int,int,int,int,int,int,bool&);
	string GetValue(string);
	Real FixedPoint(Real);
	bool WriteResults();
	bool PutChi(Real);
	Real zero_J0(Real,Real);

};
#endif
