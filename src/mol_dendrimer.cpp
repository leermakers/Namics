#include "molecule.h"
#include "mol_dendrimer.h"


mol_dend::mol_dend(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {}


mol_dend::~mol_dend() {
}




Real mol_dend::fraction(int segnr){
if (debug) cout <<"fraction for mol_test " + name << endl;
	int Nseg=0;
	int length = mon_nr.size();
	int i=0;
	while (i<length) {
		if (segnr==mon_nr[i]) Nseg += n_mon[i]*d_mon[i];
		i++;
	}
	return 1.0*Nseg/chainlength;
}





bool mol_dend::ComputePhi() {
	if (debug) cout <<"ComputePhiDendrimer for mol_test " + name << endl;
	int N;
	int M=Lat[0]->M;
	Real* GS = new Real[3*M];

	bool success=true;
	int n_g=first_a.size();
	int slast=last_s[n_g-1];
	int s=slast;

	for (int g=n_g-1; g>=0; g--) {
			Cp(GS+2*M,UNITY,M);
			int b0=first_b[g], bN=last_b[g];
			for (int b=bN; b>=b0; b--) {
				N= n_mon[b];
				for (int k=0; k<N; k++) {
					if (s<slast) {
						Lat[0] ->propagate(Gg_f,Seg[mon_nr[b]]->G1,s+1,s,M);
					} else {
						Cp(Gg_f+slast*M,Seg[mon_nr[b]]->G1,M);
					}
					s--;
				}
			}

			Lat[0] ->propagate(Gg_f,Seg[mon_nr[b0-1]]->G1,s+1,s,M);
			Cp(GS,Gg_f+(s+1)*M,M);
			Lat[0] ->propagate(GS,UNITY,0,1,M);
			for (int k=0; k<n_arm[g]-1; k++) Times(Gg_f+s*M,Gg_f+s*M,GS+M,M); //Times(GS+2*M,GS+2*M,GS+M,M);
			s--;
	}

	GN=Lat[0]->WeightedSum(Gg_f);
	Times(rho+molmon_nr[0]*M, Gg_f, Seg[mon_nr[0]]->G1, M);

	s=0;
	Cp(GS+2*M,UNITY,M);
	Cp(GS,Gg_f+M,M);
	Lat[0] ->propagate(GS,UNITY,0,1,M);
	for (int k=0; k<n_arm[0]-1; k++) Times(GS+2*M,GS+2*M,GS+M,M);
	Times(Gg_b,GS+2*M,Seg[mon_nr[0]]->G1,M);
	for (int g=0; g<n_g; g++) {
		Cp(GS+2*M,UNITY,M);
		int b0=first_b[g], bN=last_b[g];
		for (int b=b0; b<=bN; b++) {
			N= n_mon[b];
			for (int k=0; k<N; k++) {
				Lat[0] ->propagate(Gg_b,Seg[mon_nr[b]]->G1,s%2,(s+1)%2,M);
				s++;
				Times(GS, Gg_f+s*M, Gg_b+(s%2)*M, M); Norm(GS, d_mon[b],M);
				Add(rho+molmon_nr[b]*M,GS,M);
			}
		}
		if (s<slast) {
			Cp(GS+2*M,UNITY,M);
			Cp(GS,Gg_f+(s+2)*M,M);
			Lat[0] ->propagate(GS,UNITY,0,1,M);
			for (int k=0; k<n_arm[g+1]-1; k++) Times(GS+2*M,GS+2*M,GS+M,M);
			Lat[0] ->propagate(Gg_b,Seg[mon_nr[bN+1]]->G1,s%2,(s+1)%2,M);
			s++;
			Times(GS, Gg_f+(s)*M, Gg_b+(s%2)*M, M);
			Norm(GS, d_mon[bN+1],M);
			Add(rho+molmon_nr[bN+1]*M,GS,M);
			Times(Gg_b+(s%2)*M,Gg_b+(s%2)*M,GS+2*M,M);
		}
	}
	delete [] GS;
	return success;
}



