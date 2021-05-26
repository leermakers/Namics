#include "molecule.h"
#include "mol_comb.h"


mol_comb::mol_comb(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {}


mol_comb::~mol_comb() {
}


Real mol_comb::fraction(int segnr){
if (debug) cout <<"fraction for mol_comb " + name << endl;
	int Nseg=0;
	int length = mon_nr.size();
	int i=0;
	while (i<length) {
		if (segnr==mon_nr[i]) Nseg += n_mon[i]*d_mon[i]; 
		i++;
	}
	return 1.0*Nseg/chainlength;
}



bool mol_comb::ComputePhi() {
	if (debug) cout <<"ComputePhi for mol_comb " + name << endl;
	int N;
	int M=Lat[0]->M;
	Real* GS = new Real[3*M];

	bool success=true;
	int slast=0,sfirst=0;
	int s;
	int n_arms=n_arm[0];
	int b0,bN,g;


	s=last_s[1]; //first do the arm;
	slast=s;
	b0=first_b[1]; bN=last_b[1];
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
	sfirst = s+1; //keep the reference to beginning of arm.

	s=last_s[n_arms + 2]; //then start at end of backbone.
	slast=s;
	b0=first_b[n_arms+2]; bN=last_b[n_arms+2];
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


	for (g=n_arms+1; g>=2; g--) { //do the repeat
		b0=first_b[g]; bN=last_b[g];
		for (int b=bN; b>=b0; b--) {
			N= n_mon[b];
			for (int k=0; k<N; k++) {
				Lat[0] ->propagate(Gg_f,Seg[mon_nr[b]]->G1,s+1,s,M);
				s--;
			}
		}


		Lat[0] ->propagate(Gg_f,Seg[mon_nr[b0-1]]->G1,s+1,s,M);
		Cp(GS,Gg_f+sfirst*M,M);
		Lat[0] ->propagate(GS,UNITY,0,1,M);
		Times(Gg_f+s*M,Gg_f+s*M,GS+M,M);
		s--;
	}

	b0=first_b[0]; bN=last_b[0];
	Lat[0] ->propagate(Gg_f,Seg[mon_nr[bN]]->G1,s+1,last_s[0],M);
	s= last_s[0];
	for (int b=bN; b>=b0; b--) {
		N= n_mon[b];
		for (int k=0; k<N; k++) {
			if (s<last_s[0]) {
				Lat[0] ->propagate(Gg_f,Seg[mon_nr[b]]->G1,s+1,s,M);
			}
			s--;
		}
	}

	GN=Lat[0]->WeightedSum(Gg_f);
	Times(rho+molmon_nr[0]*M, Gg_f, Seg[mon_nr[0]]->G1, M);
	s=-1;

	b0=first_b[0]; bN=last_b[0];
	for (int b=b0; b<=bN; b++) {
		N= n_mon[b];
		for (int k=0; k<N; k++) {
			if (s>-1) {
				Lat[0] ->propagate(Gg_b,Seg[mon_nr[b]]->G1,s%2,(s+1)%2,M);
				AddTimes(rho+molmon_nr[b]*M, Gg_f+(s+1)*M,Gg_b+((s+1)%2)*M, M);
			} else Cp(Gg_b,Seg[mon_nr[b]]->G1,M);
			s++;
		}
	}

	for (g=2; g<=n_arms+1; g++) {
		Lat[0] ->propagate(Gg_b,Seg[mon_nr[first_b[g]-1]]->G1,s%2,(s+1)%2,M);
		Cp(GS,Gg_b+(s+1)%2*M,M); //opslag
		Lat[0] ->propagate(Gg_f,UNITY,first_s[g],first_s[g]-1,M);
 		Times(Gg_b+((sfirst-1)%2)*M,Gg_f+(first_s[g]-1)*M,GS,M);
		b0=first_b[1]; bN=last_b[1]; s=sfirst-1;
		for (int b=b0; b<=bN; b++) {
			N= n_mon[b];
			for (int k=0; k<N; k++) {
				Lat[0] ->propagate(Gg_b,Seg[mon_nr[b]]->G1,s%2,(s+1)%2,M);
				AddTimes(rho+molmon_nr[b]*M, Gg_f+(s+1)*M,Gg_b+((s+1)%2)*M, M);
				s++;

			}
		}
		Cp(GS+M,Gg_f+sfirst*M,M);
		Lat[0] ->propagate(GS,UNITY,1,2,M);
		Times(Gg_b+((first_s[g]-1)%2)*M,GS+2*M,GS,M);
		Cp(GS,Gg_f+first_s[g]*M,M); Lat[0]->propagate(GS,Seg[mon_nr[first_b[g]-1]]->G1,0,1,M);
		AddTimes(rho+molmon_nr[first_b[g]-1]*M, GS+M,Gg_b+((first_s[g]-1)%2)*M, M);
		b0=first_b[g];
		bN=last_b[g];
		s=first_s[g]-1;
		for (int b=b0; b<=bN; b++) {
			N= n_mon[b];
			for (int k=0; k<N; k++) {
				Lat[0] ->propagate(Gg_b,Seg[mon_nr[b]]->G1,s%2,(s+1)%2,M);
				AddTimes(rho+molmon_nr[b]*M, Gg_f+(s+1)*M,Gg_b+((s+1)%2)*M, M);
				s++;
			}
		}
	}

	b0=first_b[n_arms+2]; bN=last_b[n_arms+2];
	for (int b=b0; b<=bN; b++) {
		N= n_mon[b];
		for (int k=0; k<N; k++) {
			Lat[0] ->propagate(Gg_b,Seg[mon_nr[b]]->G1,s%2,(s+1)%2,M);
			AddTimes(rho+molmon_nr[b]*M, Gg_f+(s+1)*M,Gg_b+((s+1)%2)*M, M);
			s++;
		}
	}

	delete [] GS;
	return success;
}


