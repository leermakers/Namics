#ifndef LGRAD2_H
#define LGRAD2_H
class LGrad2 : public Lattice 
{
	public:	LGrad2(vector<Input*> In_,string name_);
	~LGrad2();

	void ComputeLambdas(void);
	bool PutM();
	void TimesL(Real*);
	void DivL(Real*);
	Real Moment(Real*,Real,int);
	Real WeightedSum(Real*);
	void vtk(string, Real*,string,bool);
	void PutProfiles(FILE*,vector<Real*>,bool);
	void Side(Real *, Real *, int);
	void propagate(Real*,Real*, int, int,int);
	void LReflect(Real*,Real*,Real*);
	void UReflect(Real*,Real*,Real*);
	virtual void propagateF(Real*,Real*,Real*, int, int,int);
	virtual void propagateB(Real*,Real*,Real*, int, int,int);
	bool ReadRange(int*, int*, int&, bool&, string, int, string, string);
	bool ReadRangeFile(string,int* H_p,int&, string, string);
	bool FillMask(int*, vector<int>, vector<int>, vector<int>, string);
	bool CreateMASK(int*, int*, int*, int, bool);
	Real ComputeTheta(Real*);
	void UpdateEE(Real*, Real*,Real*);
	void UpdatePsi(Real*, Real*, Real* , Real*, int*,bool,bool);
	void UpdateQ(Real*,Real*,Real*,Real*,int*,bool);
	void remove_bounds(Real*);
	void set_bounds(Real*);
	void remove_bounds(int*);
	void set_bounds(int*);
	Real ComputeGN(Real*, int);
	void AddPhiS(Real*,Real*,Real*,int);
	void AddPhiS(Real*,Real*,Real*,Real*, Real, int);
	void Initiate(Real*,Real*,int);
	void set_bounds_x(Real*,Real*);
	void set_bounds_y(Real*,Real*);
	void set_bounds_x(Real*);
	void set_bounds_y(Real*);
};
#endif

