#ifndef LG1PLANAR_H
#define LG1PLANAR_H
#include "LGrad1.h"
class LG1Planar : public LGrad1
{
	public:	LG1Planar(vector<Input*> In_,string name_);
	~LG1Planar();

	void ComputeLambdas(void);
	Real MomentPlanar(Real*,int,Real);
	void Side(Real *, Real *, int);
	void propagate(Real*,Real*, int, int,int);
	void propagateF(Real*,Real*, Real*, int, int,int);
	void propagateB(Real*,Real*, Real*, int, int,int);
	void UpdateEE(Real*, Real*,Real*);
	void UpdatePsi(Real*, Real*, Real* , Real*, Real*,bool,bool);
	void UpdateQ(Real*,Real*,Real*,Real*,Real*,bool);
	bool PutMask(Real* ,vector<int>,vector<int>,vector<int>,int);
	Real DphiDt(Real*,Real*,Real*,Real*,Real*,Real*,Real,Real);
};
#endif

