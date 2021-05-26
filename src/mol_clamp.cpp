#include "molecule.h"
#include "mol_clamp.h"


mol_clamp::mol_clamp(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {}


mol_clamp::~mol_clamp() {
}


bool mol_clamp::ComputePhi(){
	if (debug) cout <<"ComputePhi for mol_clamp " + name << endl;
	bool success=true;
	int M=Lat[0]->M;
	int m=0;
	if (freedom=="clamped") m=Lat[0]->m[Seg[mon_nr[0]]->clamp_nr];
	int blocks=mon_nr.size();
	Zero(rho,m*n_box*MolMonList.size());
	int s=1;
	if (save_memory) {
		Cp(Gs,mask1,m*n_box);
	} else {
		Cp(Gg_f,mask1,m*n_box);
	}
	for (int i=1; i<blocks-1; i++) {
		Lat[0]->DistributeG1(Seg[mon_nr[i]]->G1,g1,Bx,By,Bz,n_box);

		propagate_forward(g1,s,i,0,m*n_box);
	}
	if (save_memory) {
		int k=last_stored[blocks-2];
		int N=memory[n_mon.size()-1];
		Lat[0]->propagate(Gg_f,mask2,k,N-1,m*n_box);
		Lat[0]->ComputeGN(gn,Gg_f,H_Bx,H_By,H_Bz,H_Px2,H_Py2,H_Pz2,N-1,n_box);
	} else {
		Lat[0]->propagate(Gg_f,mask2,s-1,s,m*n_box);
		Lat[0]->ComputeGN(gn,Gg_f,H_Bx,H_By,H_Bz,H_Px2,H_Py2,H_Pz2,chainlength-1,n_box);
	}
	s=chainlength-1;
	Cp(Gg_b+(s%2)*m*n_box,mask2,m*n_box);
	if (save_memory) Cp(Gg_b+((s-1)%2)*m*n_box,Gg_b+(s%2)*m*n_box,m*n_box);
	s--;
	for (int i=blocks-2; i>0; i--) {
		Lat[0]->DistributeG1(Seg[mon_nr[i]]->G1,g1,Bx,By,Bz,n_box);

		propagate_backward(g1,s,i,0,m*n_box);
	}
	for (size_t i = 1; i < MolMonList.size(); i++ )
	{
		Lat[0]->CollectPhi(phi+M*i,gn,rho+m*n_box*i,Bx,By,Bz,n_box);
	}

	return success;
}







