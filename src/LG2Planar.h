#ifndef LG2PLANAR_H
#define LG2PLANAR_H
#include "LGrad2.h"
class LG2Planar : public LGrad2 
{
	public:	LG2Planar(vector<Input*> In_,string name_);
	~LG2Planar();

	void ComputeLambdas(void);
	void Side(Real *, Real *, int);
	void propagateF(Real*, Real*, Real*, int, int, int);
	void propagateB(Real*, Real*, Real*, int, int, int);
	void propagate(Real*,Real*, int, int,int);
	void UpdateEE(Real*, Real*,Real*);
	void UpdatePsi(Real*, Real*, Real* , Real*, int*,bool,bool);
	void UpdateQ(Real*,Real*,Real*,Real*,int*,bool);
	Real DphiDt(Real*,Real*,Real*,Real*,Real*,Real*,Real,Real);
};
#endif

