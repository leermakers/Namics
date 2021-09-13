#ifndef MOL_PREVIEWxH
#define MOL_PREVIEWxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "alias.h"
#include "lattice.h"
#include "tools.h"
class mol_preview {
public:
	mol_preview(vector<Input*>,vector<Lattice*>,vector<Segment*>,string);
virtual ~mol_preview();

	string name;
	bool all_molecule;
	vector<Input*> In;
	vector<Segment*> Seg;
	vector<Lattice*> Lat;
	vector<Alias*> Al;
	vector<int> MolMonList;
	vector<int> MolAlList;
	int start;
	int n_mol;
	int mol_nr;
	int n_box;
	int var_al_nr;
	Real Mu;
	Real theta;
	Real theta_range,n_range;
	Real thata_in_range;
	Real phibulk;
	string freedom;
	MoleculeType MolType;
	Real n;
	Real GN,GN1,GN2;
	Real norm;
	Real phi1,phiM,width,Dphi,pos_interface,phi_av;
	int Markov;

	int chainlength,N;
	bool save_memory;
	bool compute_phi_alias;
	bool sym_dend;
	string composition;
	vector<int> Gnr; //generation-number
	vector<int> first_s;
	vector<int> last_s;
	vector<int> first_b;
	vector<int> last_b;
	vector<int> first_a;
	vector<int> last_a;
	vector<int> n_arm;
	vector<int> mon_nr;
	vector<int> n_mon;
	vector<int> d_mon; //degeneracy of mon : for dendrimer case.
	vector<int> molmon_nr;
	vector<int> memory;
	vector<int> last_stored;
	vector<Real> mu_state;
	vector<Real> block;
	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;

	int *Bx;
	int *By;
	int *Bz;
	int *Px1;
	int *Py1;
	int *Pz1;
	int *Px2;
	int *Py2;
	int *Pz2;
	int *H_Bx;
	int *H_By;
	int *H_Bz;
	int *H_Px1;
	int *H_Py1;
	int *H_Pz1;
	int *H_Px2;
	int *H_Py2;
	int *H_Pz2;
	Real *phi;
	//Real *G1;
	//Real *u;
	Real *rho;
	Real *H_phi;
	Real *H_mask1;
	Real *H_mask2;
	//Real* H_u;
	Real *H_gn;
	Real *mask1;
	Real *mask2;
	int *R_mask;
	Real *gn;
	Real *g1;
	Real *phitot;
	Real *H_phitot;
	Real *Gg_f;
	Real *Gg_b;
	Real *Gs;
	Real *UNITY;
	int tag_segment;
	int Var_steps;
	Real Var_step;
	Real Var_end_value;
	Real Var_start_value;
	Real Var_start_search_value;
	int num_of_steps;
	int Var_target;
	int Var_scan_value;
	int Var_search_value;
	string Var_type;
	string scale;
	Real Var_target_value;
	int n_generations;
	int GetChainlength(void);
	bool MakeMonList(void);
	int GetAlNr(string s);
	int GetMonNr(string s);
	string GetValue(string);
	bool DeleteAl();
	bool CheckInput(int);
	bool ExpandAlias(vector<string>,string&);
	bool ExpandBrackets(string&);
	bool Interpret(string,int);
	bool GenerateTree(string,int,int&,vector<int>,vector<int>);
	bool Decomposition(string);
};

#endif
