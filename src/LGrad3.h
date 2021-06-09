#ifndef LGRAD3_H
#define LGRAD3_H
class LGrad3 : public Lattice 
{
	public:	LGrad3(vector<Input*> In_,string name_);
	~LGrad3();

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
	void propagateF(Real*,Real*,Real*, int, int,int);
	void propagateB(Real*,Real*,Real*, int, int,int);
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
	Real ComputeGN(Real*,int);
	void AddPhiS(Real*,Real*,Real*,int);
	void AddPhiS(Real*,Real*,Real*,Real*,Real,int);
	void Initiate(Real*,Real*,int);
	void set_bounds_x(Real*,Real*,int,int);
	void set_bounds_y(Real*,Real*,int,int);
	void set_bounds_z(Real*,Real*,int,int);
	void set_bounds_x(Real*,int,int);
	void set_bounds_y(Real*,int,int);
	void set_bounds_z(Real*,int,int);
};
#endif

