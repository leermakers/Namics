#ifndef MOL_ASYM_DENDxH
#define MOL_ASYM_DENDxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "alias.h"
#include "lattice.h"
#include "tools.h"
class mol_asym_dend : public Molecule
{
	public: mol_asym_dend(vector<Input*>,vector<Lattice*>,vector<Segment*>,string);
	~mol_asym_dend();

	Real fraction(int);
	bool ComputePhi();
	Real* Forward();
	bool Backward(int,int,int); 
	Real* Forward2ndO();
	bool Backward2ndO(int,int,int); 

};

#endif
