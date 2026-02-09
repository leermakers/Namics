#include "molecule.h"
#include "mol_dendrimer.h"


mol_dend::mol_dend(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {}


mol_dend::~mol_dend() {
}


Real mol_dend::fraction(int segnr){
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


bool mol_dend::BackAndForth2ndO() {
if (debug) cout <<"BackAndForth2ndO for mol_dend " + name << endl;

	int N;
	int M=lat->M;
#ifdef CUDA
	Real* GS = (Real*)AllOnDev(3*M);
#else
	Real* GS = new Real[3*M];
#endif

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
					if (b==bN && k==0) {
						lat->Terminate(GS,Gg_f+(s+1)*M*size,Markov,M);
						lat->propagate(GS,Seg[mon_nr[b]]->G1,0,1,M);
						lat->Initiate(Gg_f+s*M*size,GS+M,Markov,M);
					} else {
						lat->propagateF(Gg_f,Seg[mon_nr[b]]->G1,P,s+1,s,M);
					}
				} else {
					lat->Initiate(Gg_f+slast*M*size,Seg[mon_nr[b]]->G1,Markov,M);
				}
				s--;
			}
		}

		lat->Terminate(GS,Gg_f+(s+1)*M*size,Markov,M);
		lat->propagate(GS,UNITY,0,1,M);
		lat->propagate(GS,Seg[mon_nr[b0-1]]->G1,0,2,M);
		for (int k=0; k<n_arm[g]-1; k++) Times(GS+2*M,GS+2*M,GS+M,M); //Times(GS+2*M,GS+2*M,GS+M,M);
		lat->Initiate(Gg_f+s*M*size,GS+2*M,Markov,M);
		s--;
	}

	GN=lat->ComputeGN(Gg_f,Markov,M);
	lat->Initiate(Gg_b,Seg[mon_nr[0]]->G1,Markov,M);
	lat->AddPhiS(rho+molmon_nr[0]*M,Gg_f,Gg_b,Markov,M);

	s=0;
	Cp(GS+2*M,UNITY,M);
	//lat->Terminate(GS,Gg_f+M*size,Markov,M); //GS heeft  al Gg_f+M*size
	lat->propagate(GS,UNITY,0,1,M);
	for (int k=0; k<n_arm[0]-1; k++) Times(GS+2*M,GS+2*M,GS+M,M);
	Times(GS+2*M,GS+2*M,Seg[mon_nr[0]]->G1,M);
	lat->Initiate(Gg_b,GS+2*M,Markov,M);
	for (int g=0; g<n_g; g++) {
		int b0=first_b[g], bN=last_b[g];
		for (int b=b0; b<=bN; b++) {
			N= n_mon[b];
			for (int k=0; k<N; k++) {
				if (k==0&& b==b0) {
					lat->Terminate(GS,Gg_b+(s%2)*M*size,Markov,M);
					lat->propagate(GS,Seg[mon_nr[b]]->G1,0,1,M);
					lat->Initiate(Gg_b+((s+1)%2)*M*size,GS+M,Markov,M);

				} else {
					lat->propagateB(Gg_b,Seg[mon_nr[b]]->G1,P,s%2,(s+1)%2,M);
				}
				s++;
				lat->AddPhiS(rho+molmon_nr[b]*M, Gg_f+s*M*size, Gg_b+(s%2)*M*size,d_mon[b],Markov, M);
			}
		}
		if (s<slast) {
			Cp(GS+2*M,UNITY,M);
			lat->Terminate(GS,Gg_f+(s+2)*M*size,Markov,M);
			lat->propagate(GS,UNITY,0,1,M);
			for (int k=0; k<n_arm[g+1]-1; k++) Times(GS+2*M,GS+2*M,GS+M,M);
			lat->Terminate(GS,Gg_b+(s%2)*M*size,Markov,M);
			lat->propagate(GS,Seg[mon_nr[bN+1]]->G1,0,1,M);
			lat->Initiate(Gg_b+((s+1)%2)*M*size,GS+M,Markov,M);
			s++;
			lat->AddPhiS(rho+molmon_nr[bN+1]*M, Gg_f+(s)*M*size, Gg_b+(s%2)*M*size,d_mon[bN+1],Markov, M);
			Times(GS+M,GS+M,GS+2*M,M);
			lat->Initiate(Gg_b+(s%2)*M*size,GS+M,Markov,M);
		}
	}
#ifdef CUDA
	cudaFree(GS);
#else
	delete [] GS;
#endif
	return success;
}


bool mol_dend::BackAndForth() {
if (debug) cout <<"BackAndForth for mol_dend " + name << endl;

	int N;
	int M=lat->M;
#ifdef CUDA
	Real* GS = (Real*)AllOnDev(3*M);
#else
	Real* GS = new Real[3*M];
#endif

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
						lat->propagate(Gg_f,Seg[mon_nr[b]]->G1,s+1,s,M);
					} else {
						lat->Initiate(Gg_f+slast*M*size,Seg[mon_nr[b]]->G1,Markov,M);
					}
					s--;
				}
			}

			lat->Terminate(GS,Gg_f+(s+1)*M*size,Markov,M);
			lat->propagate(GS,UNITY,0,1,M);
			lat->propagate(GS,Seg[mon_nr[b0-1]]->G1,0,2,M);
			for (int k=0; k<n_arm[g]-1; k++) Times(GS+2*M,GS+2*M,GS+M,M); //Times(GS+2*M,GS+2*M,GS+M,M);
			lat->Initiate(Gg_f+s*M*size,GS+2*M,Markov,M);
			s--;
	}

	GN=lat->ComputeGN(Gg_f,Markov,M);
	lat->Initiate(Gg_b,Seg[mon_nr[0]]->G1,Markov,M);
	lat->AddPhiS(rho+molmon_nr[0]*M,Gg_f,Gg_b,Markov,M);

	s=0;
	Cp(GS+2*M,UNITY,M);
	lat->Terminate(GS,Gg_f+M*size,Markov,M);
	lat->propagate(GS,UNITY,0,1,M);
	for (int k=0; k<n_arm[0]-1; k++) Times(GS+2*M,GS+2*M,GS+M,M);
	Times(GS+2*M,GS+2*M,Seg[mon_nr[0]]->G1,M);
	lat->Initiate(Gg_b,GS+2*M,Markov,M);
	for (int g=0; g<n_g; g++) {
		Cp(GS+2*M,UNITY,M);
		int b0=first_b[g], bN=last_b[g];
		for (int b=b0; b<=bN; b++) {
			N= n_mon[b];
			for (int k=0; k<N; k++) {
				lat->propagate(Gg_b,Seg[mon_nr[b]]->G1,s%2,(s+1)%2,M);
				s++,
				lat->AddPhiS(rho+molmon_nr[b]*M, Gg_f+s*M*size, Gg_b+(s%2)*M*size,d_mon[b],Markov,M); 			}
		}
		if (s<slast) {
			Cp(GS+2*M,UNITY,M);
			lat->Terminate(GS,Gg_f+(s+2)*M*size,Markov,M);
			lat->propagate(GS,UNITY,0,1,M);
			for (int k=0; k<n_arm[g+1]-1; k++) Times(GS+2*M,GS+2*M,GS+M,M);
			lat->propagate(Gg_b,Seg[mon_nr[bN+1]]->G1,s%2,(s+1)%2,M);
			s++;
			lat->AddPhiS(rho+molmon_nr[bN+1]*M, Gg_f+(s)*M*size, Gg_b+(s%2)*M*size,d_mon[bN+1],Markov, M);
			lat->Terminate(GS,Gg_b+(s%2)*M*size,Markov,M);
			Times(GS,GS,GS+2*M,M);
			lat->Initiate(Gg_b+(s%2)*M,GS,Markov,M);
		}
	}
#ifdef CUDA
	cudaFree(GS);
#else
	delete [] GS;
#endif
	return success;
}


bool mol_dend::ComputePhi() {
	if (debug) cout <<"ComputePhi for mol_dend " + name << endl;
	bool success=true;
	if (Markov ==2) success=BackAndForth2ndO(); else success=BackAndForth();
	return success;
}




