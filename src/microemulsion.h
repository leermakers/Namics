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
 	CP ControlParameter;
	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	Real *phi_oil;
	int surfactant;
	int co_solvent;
	int oil;
	int water;
	int monA;
	int monB;
	int monC;
	int monD;
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
	Real g_tolerance;
	Real j_tolerance;
	int g_info;
	int j_info;
	bool compute_kappa;
	string co_solvent_freedom;
	Real GS0;
	Real GC0;
	Real ini_oil, ini_water, ini_surf, ini_chi;
	bool previous_guess;

	bool CheckInput(int);
	bool Doit(Real*,string,vector<string>,vector<string>,bool,int,int,int,int,int,int,int,int,int,int,bool&);
	string GetValue(string);
	bool FixedPoint(Real, Real);
	bool WriteResults();
	bool PutChi(Real);
	Real get_gamma(Real, Real);
	Real zero_gamma(Real, Real);
	Real zero_J0(Real,Real,Real);
	Real ConvertSurfactantXtoT(Real);
	Real ConvertCoSolventXtoT(Real);
	Real ConvertSurfactantTtoX(Real);
	Real ConvertCoSolventTtoX(Real);
	bool SlipInSphericalCoordinates();
};
#endif
