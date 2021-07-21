#ifndef SOLVE_BFGSxH
#define SOLVE_BFGSxH
#include "solve_scf.h"


class Solve_BFGS : public Solve_scf 
{
public:

	Solve_BFGS(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<State*>,vector<Reaction*>,vector<Molecule*>,vector<System*>,vector<Variate*>,string);

	~Solve_BFGS();

	

	Real operator()(Vector& , Vector& );


};
#endif
