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

bool mol_comb::GoBackAndForth2ndO() { 
	int N;
	int M=Lat[0]->M;
	Real* GB = new Real[2*M]; Zero(GB,2*M);

	bool success=true;
	int slast=0,sfirst=0;
	int s;
	int n_arms=n_arm[0];
	int b0,bN,g; 

	s=last_s[1]+1; //first do the arm;
	slast=s;
	b0=first_b[1]; bN=last_b[1];
	for (int b=bN; b>=b0; b--) {
		N= n_mon[b];
		for (int k=0; k<N; k++) {
			if (s<slast) {
				Lat[0] ->propagateF(Gg_f,Seg[mon_nr[b]]->G1,P,s+1,s,M);
			} else {
				Lat[0]->Initiate(Gg_f+slast*M*size,Seg[mon_nr[b]]->G1,M);
			}
			s--;
		}
	}

	Lat[0] ->Terminate(GB+M,Gg_f+(s+1)*M*size,M);
	Lat[0] ->propagate(GB,UNITY,1,0,M); //side chain is freely joined


	sfirst = s+1; //keep the reference to beginning of arm.
	s=last_s[n_arms + 2]; //then start at end of backbone.
	slast=s;
	b0=first_b[n_arms+2]; bN=last_b[n_arms+2];
	for (int b=bN; b>=b0; b--) {
		N= n_mon[b];
		for (int k=0; k<N; k++) {
			if (s<slast) {
				Lat[0] ->propagateF(Gg_f,Seg[mon_nr[b]]->G1,P,s+1,s,M);
			} else {
				Lat[0]->Initiate(Gg_f+slast*M*size,Seg[mon_nr[b]]->G1,M);
			}
			s--;
		}
	}

	for (g=n_arms+1; g>=2; g--) { //do the repeat
		b0=first_b[g]; bN=last_b[g];
		for (int b=bN; b>=b0; b--) {
			N= n_mon[b];
			for (int k=0; k<N; k++) {  //for the spacer
				Lat[0] ->propagateF(Gg_f,Seg[mon_nr[b]]->G1,P,s+1,s,M);
				s--;
			}
		}
		if (g==2) {
			Lat[0] ->propagateF(Gg_f,Seg[mon_nr[b0-1]]->G1,P,s+1,(last_s[0]+1),M); //go to the branching point
			for (int k=0; k<size; k++) Times(Gg_f+(last_s[0]+1)*M*size+k*M,Gg_f+(last_s[0]+1)*M*size+k*M,GB,M); //connect side
		} else {
			Lat[0] ->propagateF(Gg_f,Seg[mon_nr[b0-1]]->G1,P,s+1,s,M); //go to the branching point
			for (int k=0; k<size; k++) Times(Gg_f+s*M*size+k*M,Gg_f+s*M*size+k*M,GB,M); //connect side
		}
		s--;
	}
	b0=first_b[0]; bN=last_b[0];s= last_s[0];
	
	for (int b=bN; b>=b0; b--) { //terminal part of main chain
		N= n_mon[b];
		for (int k=0; k<N; k++) {
			Lat[0] ->propagateF(Gg_f,Seg[mon_nr[b]]->G1,P,s+1,s,M);
			s--;
		}
	}
	GN=Lat[0]->ComputeGN(Gg_f,M);
	success = GN>=0;
	Lat[0]->Initiate(Gg_b,Seg[mon_nr[0]]->G1,M);
	Lat[0]->AddPhiS(rho+molmon_nr[0]*M, Gg_f, Gg_b, M);
	s=-1;
	b0=first_b[0]; bN=last_b[0];
	for (int b=b0; b<=bN; b++) { //start with main chain
		N= n_mon[b];
		for (int k=0; k<N; k++) {
			if (s<0) {
			} else {
				Lat[0]->propagateB(Gg_b,Seg[mon_nr[b]]->G1,P,s%2,(s+1)%2,M);
				Lat[0]->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M*size,Gg_b+((s+1)%2)*M*size, M);
			}
			s++;
		}
	}

	for (g=2; g<=n_arms+1; g++) { //only do the spacer and keep semiflexibility
		Lat[0] ->propagateB(Gg_b,Seg[mon_nr[first_b[g]-1]]->G1,P,s%2,(s+1)%2,M);
		Lat[0] ->AddPhiS(rho+molmon_nr[first_b[g]-1]*M,Gg_b+((s+1)%2)*M*size,Gg_f+(s+1)*M*size,M); //this should add phi's for the branch point
		Times(Gg_f+(s+1)*M*size,Gg_b+((s+1)%2)*M*size,Gg_f+(s+1)*M*size,M*size); 
		for (int k=0; k<size; k++) {
			Div(Gg_f+(s+1)*M*size+k*M,GB,M);
			Div(Gg_f+(s+1)*M*size+k*M,Seg[mon_nr[first_b[g]-1]]->G1 ,M); //getting ready for phi side computation!
			Times(Gg_b+((first_s[g]-1)%2)*M*size+k*M,Gg_b+((s+1)%2)*M*size+k*M,GB,M); //getting ready to propagate along backbone	

		}
		
		b0=first_b[g];
		bN=last_b[g];
		s=first_s[g]-1;
		for (int b=b0; b<=bN; b++) { //spacer
			N= n_mon[b];
			for (int k=0; k<N; k++) {
				Lat[0] ->propagateB(Gg_b,Seg[mon_nr[b]]->G1,P,s%2,(s+1)%2,M);
				Lat[0] ->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M*size,Gg_b+((s+1)%2)*M*size, M);
				s++;
			}
		}
	}
	b0=first_b[n_arms+2]; bN=last_b[n_arms+2];
	for (int b=b0; b<=bN; b++) {
		N= n_mon[b];
		for (int k=0; k<N; k++) { //terminal part of main chain
			Lat[0] ->propagateB(Gg_b,Seg[mon_nr[b]]->G1,P,s%2,(s+1)%2,M);
			Lat[0] ->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M*size,Gg_b+((s+1)%2)*M*size, M);
			s++;
		}
	}

	for (g=2; g<=n_arms+1; g++) { //computing phi's for sides
		s=sfirst-1;
 		if (g==2) {
			Lat[0]->Terminate(GB,Gg_f+(last_s[0]+1)*M*size,M);
		} else {
			Lat[0]->Terminate(GB,Gg_f+(first_s[g]-1)*M*size,M);

		}		
		
		b0=first_b[1]; bN=last_b[1]; 
		
		for (int b=b0; b<=bN; b++) {
			N= n_mon[b];
			for (int k=0; k<N; k++) {
				if (k==0&&b==b0) {
					Lat[0] ->propagate(GB,Seg[mon_nr[b]]->G1,0,1,M);
					Lat[0] ->Initiate(Gg_b+((s+1)%2)*M*size,GB+M,M);
					Lat[0] ->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M*size,Gg_b+((s+1)%2)*M*size, M);

				} else {
					Lat[0] ->propagateB(Gg_b,Seg[mon_nr[b]]->G1,P,s%2,(s+1)%2,M);
					Lat[0] ->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M*size,Gg_b+((s+1)%2)*M*size, M);
				}
				s++;
			}
		}
	}

	delete [] GB;
	return success;

}


bool mol_comb::GoBackAndForth() { //somewhat different way to deal with branching point but it's okay as well. 
	bool success = true;
	int N;
	int M=Lat[0]->M;
	Real* GS = new Real[3*M];

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
	success = GN>=0;
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


bool mol_comb::ComputePhi() {
	if (debug) cout <<"ComputePhi for mol_comb " + name << endl;
	bool success=true;
	if (Markov==2) success = GoBackAndForth2ndO(); else success = GoBackAndForth(); 	
	return success;
}


