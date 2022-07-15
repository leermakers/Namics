#ifndef LG1PLANAR_H
#define LG1PLANAR_H
#include "LGrad1.h"
class LG1Planar : public LGrad1
{
	public:	LG1Planar(vector<Input*> In_,string name_);
	~LG1Planar();

	void ComputeLambdas(void);
	void Side(Real *, Real *, int);
	void propagate(Real*,Real*, int, int,int);
	void propagateF(Real*,Real*, Real*, int, int,int);
	void propagateB(Real*,Real*, Real*, int, int,int);
	void UpdateEE(Real*, Real*,Real*);
	void UpdatePsi(Real*, Real*, Real* , Real*, int*,bool,bool);
	void UpdateQ(Real*,Real*,Real*,Real*,int*,bool);
	bool PutMask(int* ,vector<int>,vector<int>,vector<int>,int);
	Real DphiDt(Real*,Real*,Real*,Real*,Real*,Real*,Real,Real);
};
#endif

