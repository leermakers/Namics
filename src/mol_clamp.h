#ifndef MOL_CLAMPxH
#define MOL_CLAMPxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "alias.h"
#include "lattice.h"
#include "tools.h"
class mol_clamp : public Molecule
{
	public: mol_clamp(vector<Input*>,vector<Lattice*>,vector<Segment*>,string);
	~mol_clamp();

	//Real fraction(int);
	//Real* ForwardBra(int generation, int &s);
	//void BackwardBra(Real*, int, int&);
	//bool ComputePhi(Real*,int); 
	bool ComputePhi();
	//bool ComputeClampLin();
	//bool ComputePhiLin();
	//bool ComputePhiBra();
	//bool ComputePhiDendrimer();
	//bool ComputePhiComb();
	//bool ComputePhiRing();
};

#endif
