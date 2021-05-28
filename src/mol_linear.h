#ifndef MOL_LINEARxH
#define MOL_LINEARxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "alias.h"
#include "lattice.h"
#include "tools.h"
class mol_linear : public Molecule
{
	public: mol_linear(vector<Input*>,vector<Lattice*>,vector<Segment*>,string);
	~mol_linear();

	bool ComputePhi();

};

#endif
