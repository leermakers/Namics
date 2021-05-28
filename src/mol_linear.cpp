#include "molecule.h"
#include "mol_linear.h"


mol_linear::mol_linear(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {}


mol_linear::~mol_linear() {
}


bool mol_linear::ComputePhi() {
if (debug) cout <<"ComputePhi in mol_linear " << endl;
	int b0 = first_b[0];
	int bN = last_b[0];
	int M=Lat[0]->M; 
	bool success=true;
	int s=0;
	Real* Glast=NULL;
 
	if (Markov ==2) 
	for (int k = b0; k<=bN ; ++k) Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,P,0,M);
	else 
	for (int k = b0; k<=bN ; ++k) Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,0,M);
 
	GN=Lat[0]->ComputeGN(Glast,M); 
 
	s--;
	if (save_memory) {
		//Cp(Gg_b,Seg[mon_nr[last_b[0]]]->G1,M); //FL
		Lat[0]->Initiate(Gg_b,Seg[mon_nr[last_b[0]]]->G1,M);
		//Cp(Gg_b+M,Seg[mon_nr[last_b[0]]]->G1,M); //FL
		Lat[0]->Initiate(Gg_b+M*size,Seg[mon_nr[last_b[0]]]->G1,M); 
	} 
 
	if (Markov==2) 	
	for (int k = bN ; k >= b0 ; k--) propagate_backward(Seg[mon_nr[k]]->G1,s,k,P,0,M);
	else 
	for (int k = bN ; k >= b0 ; k--) propagate_backward(Seg[mon_nr[k]]->G1,s,k,0,M);

 
	return success;
}




