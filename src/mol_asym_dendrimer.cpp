#include "molecule.h"
#include "mol_asym_dendrimer.h"


mol_asym_dend::mol_asym_dend(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {}


mol_asym_dend::~mol_asym_dend() {
}


Real mol_asym_dend::fraction(int segnr){
if (debug) cout <<"fraction for mol_dend " + name << endl;
	int Nseg=0;
	int length = mon_nr.size();
	int i=0;
	while (i<length) {
		if (segnr==mon_nr[i]) Nseg += n_mon[i]*d_mon[i];
		i++;
	}
	return 1.0*Nseg/chainlength;
}

Real* mol_asym_dend::Forward2ndO() {
if (debug) cout <<"Forward2ndO for mol_asym_dend " + name << endl;

	int N;
	int M=Lat[0]->M;
	Real* GS = new Real[3*M];

	int n_g=first_a.size();
	int slast=0;
	int s=0;

	for (int g=n_g-1; g>=0; g--) {
		s=last_s[last_a[g]];
		Cp(GS+2*M,UNITY,M);
		for (int a=last_a[g]; a>=first_a[g]; a--) {
			int b0=first_b[a], bN=last_b[a];
			for (int b=bN; b>=b0; b--) {
				N= n_mon[b];
				for (int k=0; k<N; k++) {
					if (g==n_g-1 && s==last_s[a]) {
						Lat[0]->Initiate(Gg_f+s*M*size,Seg[mon_nr[b]]->G1,Markov,M);
					} else {
						if (b==bN &&k==0) {
							Lat[0]->Terminate(GS,Gg_f+slast*M*size,Markov,M);
							Lat[0]->propagate(GS,Seg[mon_nr[b]]->G1,0,1,M);
							Lat[0]->Initiate(Gg_f+s*M*size,GS+M,Markov,M);
						} else {
							Lat[0] ->propagateF(Gg_f,Seg[mon_nr[b]]->G1,P,s+1,s,M);
						}
					}
					s--;
				}
			}

			Lat[0]->Terminate(GS,Gg_f+(s+1)*M*size,Markov,M);
			Lat[0]->propagate(GS,UNITY,0,1,M);
			for (int k=0; k<n_arm[a]; k++) Times(GS+2*M,GS+2*M,GS+M,M);
		}
		Times(GS+2*M,GS+2*M,Seg[mon_nr[first_b[first_a[g]]-1]]->G1,M);
		Lat[0]->Initiate(Gg_f+s*M*size,GS+2*M,Markov,M);
		slast = s;
	}
	delete [] GS;
	return Gg_f;
}

bool mol_asym_dend::Backward2ndO(int g,int n_repeats, int ss) {
if (debug) cout <<"Backward2ndO for mol_asym_dend " + name << endl;

	int M=Lat[0]->M;
	int N;
	Real* GX = new Real[M*size];
	Real* GS = new Real[2*M];
	int n_g=first_a.size();
	bool success=true;
	int s=0;

	if (g==0) {
		Lat[0]->Initiate(Gg_b,Seg[mon_nr[0]]->G1,Markov,M);
	} else {
		s= first_s[first_a[g]]-1;
		Lat[0]->Terminate(GS,Gg_b+((s-1)%2)*M*size,Markov,M);
		Lat[0]->propagate(GS,Seg[mon_nr[first_b[first_a[g]]-1]]->G1,0,1,M);
		Lat[0]->Initiate(Gg_b+(s%2)*M*size,GS+M,Markov,M);
	}

	Times(GX,Gg_b+(s%2)*M*size, Gg_f+s*M*size,M*size);
	Lat[0]->AddPhiS(rho+molmon_nr[first_b[first_a[g]]-1]*M,Gg_f+s*M*size,Gg_b+(s%2)*M*size,n_repeats,Markov,M);

	for (int k=0; k<size; k++) Div(GX+k*M,Seg[mon_nr[first_b[first_a[g]]-1]]->G1,M);

	for (int a=first_a[g]; a<=last_a[g]; a++) {
		s=first_s[a];
		Lat[0]->Terminate(GS,Gg_f+s*M*size,Markov,M);
		Lat[0]->propagate(GS,UNITY,0,1,M);
		Lat[0]->Terminate(GS,GX,Markov,M);
		Div(GS,GS+M,M);

		int b0=first_b[a], bN=last_b[a];
		for (int b=b0; b<=bN; b++) {

			N= n_mon[b];
			for (int k=0; k<N; k++) {
				if (k==0 && b==b0) {
						//Lat[0]->Terminate(GS,Gg_b+((s-1)%2)*M*size,Markov,M);
					Lat[0]->propagate(GS,Seg[mon_nr[b]]->G1,0,1,M);
					Lat[0]->Initiate(Gg_b+(s%2)*M*size,GS+M,Markov, M);
				} else {
					Lat[0]->propagateB(Gg_b,Seg[mon_nr[b]]->G1,P,(s-1)%2,s%2,M);
				}
				Lat[0]->AddPhiS(rho+molmon_nr[b]*M,Gg_f+s*M*size,Gg_b+(s%2)*M*size,n_arm[a]*n_repeats,Markov,M);
				s++;
			}
		}
		Cp(Gg_b+(s%2)*M*size,Gg_b+((s-1)%2)*M*size,M*size);
		if (g<n_g-1) Backward2ndO(g+1,n_arm[a]*n_repeats,s-1);
	}

	delete [] GS; delete [] GX;
	return success;
}



Real* mol_asym_dend::Forward() {
if (debug) cout <<"Forward for mol_asym_dend " + name << endl;

	int N;
	int M=Lat[0]->M;
	Real* GS = new Real[3*M];

	int n_g=first_a.size();
	int slast=0;
	int s=0;

	for (int g=n_g-1; g>=0; g--) {
		s=last_s[last_a[g]];
		Cp(GS+2*M,UNITY,M);
		for (int a=last_a[g]; a>=first_a[g]; a--) {
			int b0=first_b[a], bN=last_b[a];
			for (int b=bN; b>=b0; b--) {
				N= n_mon[b];
				for (int k=0; k<N; k++) {
					if (g==n_g-1 && s==last_s[a]) {
						Lat[0]->Initiate(Gg_f+s*M*size,Seg[mon_nr[b]]->G1,Markov,M);
					} else {
						if (b==bN &&k==0) {
							Lat[0] ->propagate(Gg_f,Seg[mon_nr[b]]->G1,slast,s,M);
						} else {
							Lat[0] ->propagate(Gg_f,Seg[mon_nr[b]]->G1,s+1,s,M);
						}
					}
					s--;
				}
			}

			Lat[0]->Terminate(GS,Gg_f+(s+1)*M,Markov,M);
			Lat[0]->propagate(GS,UNITY,0,1,M);
			for (int k=0; k<n_arm[a]; k++) Times(GS+2*M,GS+2*M,GS+M,M);
		}
		Times(GS+2*M,GS+2*M,Seg[mon_nr[first_b[first_a[g]]-1]]->G1,M);
		Lat[0]->Initiate(Gg_f+s*M,GS+2*M,Markov,M);
		slast = s;
	}

	delete [] GS;
	return Gg_f;
}

bool mol_asym_dend::Backward(int g,int n_repeats, int ss) {
if (debug) cout <<"Backward for mol_asym_dend " + name << endl;

	int M=Lat[0]->M;
	int N;
	Real* GX = new Real[M];
	Real* GS = new Real[2*M];
	int n_g=first_a.size();
	bool success=true;
	int s=0;

	if (g==0) {
		Lat[0]->Initiate(Gg_b,Seg[mon_nr[0]]->G1,Markov,M);
	} else {
		s= first_s[first_a[g]]-1 ;
		Lat[0]->propagate(Gg_b,Seg[mon_nr[first_b[first_a[g]]-1]]->G1,(s-1)%2,s%2,M);
	}
	Times(GX,Gg_b+(s%2)*M, Gg_f+s*M,M);
	Lat[0]->AddPhiS(rho+molmon_nr[first_b[first_a[g]]-1]*M,Gg_f+s*M,Gg_b+(s%2)*M,n_repeats,Markov,M);

	Div(GX,Seg[mon_nr[first_b[first_a[g]]-1]]->G1,M);

	for (int a=first_a[g]; a<=last_a[g]; a++) {
		s=first_s[a];
		Cp(GS,Gg_f+s*M,M);
		Lat[0]->propagate(GS,UNITY,0,1,M);
		Cp(Gg_b+((s-1)%2)*M,GX,M);
		Div(Gg_b+((s-1)%2)*M,GS+M,M);

		int b0=first_b[a], bN=last_b[a];
		for (int b=b0; b<=bN; b++) {

			N= n_mon[b];
			for (int k=0; k<N; k++) {
				Lat[0]->propagate(Gg_b,Seg[mon_nr[b]]->G1,(s-1)%2,s%2,M);
				Lat[0]->AddPhiS(rho+molmon_nr[b]*M,Gg_f+s*M,Gg_b+(s%2)*M,n_arm[a]*n_repeats,Markov,M);
				s++;
			}
		}
		Cp(Gg_b+(s%2)*M,Gg_b+((s-1)%2)*M,M);
		if (g<n_g-1) Backward(g+1,n_arm[a]*n_repeats,s-1);
	}


	delete [] GS; delete [] GX;
	return success;
}


bool mol_asym_dend::ComputePhi() {
	int M=Lat[0]->M;
	int s;
	int generation=0;
	if (debug) cout <<"ComputePhi for mol_asym_dend " + name << endl;
	bool success=true;
	if (Markov ==2) {
		GN=Lat[0]->ComputeGN(Forward2ndO(),Markov,M);
		s=0; Backward2ndO(generation,1,s);
	} else {
		GN=Lat[0]->ComputeGN(Forward(),Markov,M);
		s=0; Backward(generation,1,s);
	}
	return success;
}




