#ifndef MOL_WATERxH
#define MOL_WATERxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "alias.h"
#include "lattice.h"
#include "tools.h"
class mol_water : public Molecule
{
	public: mol_water(vector<Input*>,vector<Lattice*>,vector<Segment*>,string);
	~mol_water();
	
	void AddToGP(Real*);
	void AddToF(Real*);
	Real GetPhib1();
	bool ComputePhi();

};

#endif
