#ifndef MOL_BRANCHEDxH
#define MOL_BRANCHEDxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "alias.h"
#include "lattice.h"
#include "tools.h"
class mol_branched : public Molecule
{
	public: mol_branched(vector<Input*>,vector<Lattice*>,vector<Segment*>,string);
	~mol_branched();

	Real* ForwardBra(int generation, int &s);
	void BackwardBra(Real*, int, int&);
	Real* ForwardBra2ndO(int generation, int &s);
	void BackwardBra2ndO(Real*, int, int&);

	bool ComputePhi();

};

#endif
