#ifndef MOL_COMBxH
#define MOL_COMBxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "alias.h"
#include "lattice.h"
#include "tools.h"
class mol_comb : public Molecule
{
	public: mol_comb(vector<Input*>,vector<Lattice*>,vector<Segment*>,string);
	~mol_comb();

	Real fraction(int);
	bool ComputePhi();
	bool GoBackAndForth2ndO();
	bool GoBackAndForth();
};

#endif
