#ifndef LGRAD1_H
#define LGRAD1_H
class LGrad1 : public Lattice 
{
	public: LGrad1(vector<Input*> In_,string name_);
	virtual ~LGrad1();

	
	bool PutM();
	void TimesL(Real*);
	void DivL(Real*);
	Real Moment(Real*,Real,int);
	Real WeightedSum(Real*);
	void vtk(string, Real*,string,bool);
	void PutProfiles(FILE*,vector<Real*>,bool);
	bool ReadRange(int*, int*, int&, bool&, string, int, string, string);
	bool ReadRangeFile(string,int* H_p,int&, string, string);
	bool FillMask(int*, vector<int>, vector<int>, vector<int>, string);
	bool CreateMASK(int*, int*, int*, int, bool);
	Real ComputeTheta(Real*);
	void remove_bounds(Real*);
	void set_bounds(Real*);
	void set_bounds(Real*,Real*);
	void remove_bounds(int*);
	void set_bounds(int*);

	virtual void ComputeLambdas(void);
	virtual void UpdateEE(Real*, Real*,Real*);
	virtual void UpdatePsi(Real*, Real*, Real* , Real*, int*,bool,bool);
	virtual void UpdateQ(Real*,Real*,Real*,Real*,int*,bool);
	virtual void Side(Real *, Real *, int);
	virtual void propagate(Real*,Real*, int, int,int);
	virtual void propagateF(Real*,Real*,Real*, int, int,int);
	virtual void propagateB(Real*,Real*,Real*, int, int,int);
	virtual Real ComputeGN(Real*,int);
	virtual void AddPhiS(Real*,Real*,Real*,int);
	virtual void AddPhiS(Real*,Real*,Real*,Real*,Real,int);
	virtual void Initiate(Real*,Real*,int);

};
#endif

