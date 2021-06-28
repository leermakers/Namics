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
	if (ring) i++; 
	while (i<length) {
		if (segnr==mon_nr[i]) Nseg += n_mon[i]*d_mon[i]; 
		i++;
	}
	return 1.0*Nseg/chainlength;
}

bool mol_comb::GoBackAndForth2ndO(Real *Mask,Real* G0) { 
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
				if (ring) 
					Lat[0]->Initiate(Gg_f+slast*M*size,G0,M);
				else
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
	if (ring) {
		for (int k=0; k<size; k++) Times(Gg_f+k*M,Gg_f+k*M,Mask,M);

	}
	GN+=Lat[0]->ComputeGN(Gg_f,M);
	success = GN>=0;
	if (ring) {
		Lat[0]->Initiate(Gg_b,G0,M);
	} else {
		Lat[0]->Initiate(Gg_b,Seg[mon_nr[0]]->G1,M);

	}
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
	if (ring) {
		for (int b=b0; b<bN; b++) { //last block not needed in 'ring' because this segment is already accounted for 
			N= n_mon[b];
			for (int k=0; k<N; k++) { //terminal part of main chain
				Lat[0] ->propagateB(Gg_b,Seg[mon_nr[b]]->G1,P,s%2,(s+1)%2,M);
				Lat[0] ->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M*size,Gg_b+((s+1)%2)*M*size, M);
				s++;
			}
		}

	} else {
		for (int b=b0; b<=bN; b++) {
			N= n_mon[b];
			for (int k=0; k<N; k++) { //terminal part of main chain
				Lat[0] ->propagateB(Gg_b,Seg[mon_nr[b]]->G1,P,s%2,(s+1)%2,M);
				Lat[0] ->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M*size,Gg_b+((s+1)%2)*M*size, M);
				s++;
			}
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

bool mol_comb::GoBackAndForth(Real *Mask,Real* G0) { 
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
				Lat[0] ->propagate(Gg_f,Seg[mon_nr[b]]->G1,s+1,s,M);
			} else {
				Lat[0]->Initiate(Gg_f+slast*M*size,Seg[mon_nr[b]]->G1,M);
			}
			s--;
		}
	}

	Lat[0] ->Terminate(GB+M,Gg_f+(s+1)*M,M);
	Lat[0] ->propagate(GB,UNITY,1,0,M); //side chain is freely joined


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
				if (ring) 
					Lat[0]->Initiate(Gg_f+slast*M*size,G0,M);
				else
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
				Lat[0] ->propagate(Gg_f,Seg[mon_nr[b]]->G1,s+1,s,M);
				s--;
			}
		}
		if (g==2) {
			Lat[0] ->propagate(Gg_f,Seg[mon_nr[b0-1]]->G1,s+1,(last_s[0]+1),M); //go to the branching point
			Times(Gg_f+(last_s[0]+1)*M,Gg_f+(last_s[0]+1)*M,GB,M); //connect side
		} else {
			Lat[0] ->propagate(Gg_f,Seg[mon_nr[b0-1]]->G1,s+1,s,M); //go to the branching point
			Times(Gg_f+s*M,Gg_f+s*M,GB,M); //connect side
		}
		s--;
	}
	b0=first_b[0]; bN=last_b[0];s= last_s[0];
	
	for (int b=bN; b>=b0; b--) { //terminal part of main chain
		N= n_mon[b];
		for (int k=0; k<N; k++) {
			Lat[0] ->propagate(Gg_f,Seg[mon_nr[b]]->G1,s+1,s,M);
			s--;
		}
	}
	if (ring) {
		Times(Gg_f,Gg_f,Mask,M);

	}
	GN+=Lat[0]->ComputeGN(Gg_f,M);
	success = GN>=0;
	if (ring) {
		Lat[0]->Initiate(Gg_b,G0,M);
	} else {
		Lat[0]->Initiate(Gg_b,Seg[mon_nr[0]]->G1,M);

	}
	Lat[0]->AddPhiS(rho+molmon_nr[0]*M, Gg_f, Gg_b, M);
	s=-1;
	b0=first_b[0]; bN=last_b[0];
	for (int b=b0; b<=bN; b++) { //start with main chain
		N= n_mon[b];
		for (int k=0; k<N; k++) {
			if (s<0) {
				
			} else {
				Lat[0]->propagate(Gg_b,Seg[mon_nr[b]]->G1,s%2,(s+1)%2,M);
				Lat[0]->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M,Gg_b+((s+1)%2)*M, M);
			}
			s++;
		}
	}

	for (g=2; g<=n_arms+1; g++) { //only do the spacer and keep semiflexibility
		Lat[0] ->propagate(Gg_b,Seg[mon_nr[first_b[g]-1]]->G1,s%2,(s+1)%2,M);
		Lat[0] ->AddPhiS(rho+molmon_nr[first_b[g]-1]*M,Gg_b+((s+1)%2)*M*size,Gg_f+(s+1)*M*size,M); //this should add phi's for the branch point
		Times(Gg_f+(s+1)*M*size,Gg_b+((s+1)%2)*M*size,Gg_f+(s+1)*M*size,M*size); 
		Div(Gg_f+(s+1)*M,GB,M);
			Div(Gg_f+(s+1)*M,Seg[mon_nr[first_b[g]-1]]->G1 ,M); //getting ready for phi side computation!
			Times(Gg_b+((first_s[g]-1)%2)*M,Gg_b+((s+1)%2)*M,GB,M); //getting ready to propagate along backbone	
		
		b0=first_b[g];
		bN=last_b[g];
		s=first_s[g]-1;
		for (int b=b0; b<=bN; b++) { //spacer
			N= n_mon[b];
			for (int k=0; k<N; k++) {
				Lat[0] ->propagate(Gg_b,Seg[mon_nr[b]]->G1,s%2,(s+1)%2,M);
				Lat[0] ->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M,Gg_b+((s+1)%2)*M, M);
				s++;
			}
		}
	}
	b0=first_b[n_arms+2]; bN=last_b[n_arms+2];
	if (ring) {
		for (int b=b0; b<bN; b++) { //last block not needed in 'ring' because this segment is already accounted for 
			N= n_mon[b];
			for (int k=0; k<N; k++) { //terminal part of main chain
				Lat[0] ->propagate(Gg_b,Seg[mon_nr[b]]->G1,s%2,(s+1)%2,M);
				Lat[0] ->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M,Gg_b+((s+1)%2)*M, M);
				s++;
			}
		}

	} else {
		for (int b=b0; b<=bN; b++) {
			N= n_mon[b];
			for (int k=0; k<N; k++) { //terminal part of main chain
				Lat[0] ->propagate(Gg_b,Seg[mon_nr[b]]->G1,s%2,(s+1)%2,M);
				Lat[0] ->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M,Gg_b+((s+1)%2)*M, M);
				s++;
			}
		}
	}

	for (g=2; g<=n_arms+1; g++) { //computing phi's for sides
		s=sfirst-1;
 		if (g==2) {
			Lat[0]->Terminate(GB,Gg_f+(last_s[0]+1)*M,M);
		} else {
			Lat[0]->Terminate(GB,Gg_f+(first_s[g]-1)*M,M);
		}		
		
		b0=first_b[1]; bN=last_b[1]; 
		
		for (int b=b0; b<=bN; b++) {
			N= n_mon[b];
			for (int k=0; k<N; k++) {
				if (k==0&&b==b0) {
					Lat[0] ->propagate(GB,Seg[mon_nr[b]]->G1,0,1,M);
					Lat[0] ->Initiate(Gg_b+((s+1)%2)*M,GB+M,M);
					Lat[0] ->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M,Gg_b+((s+1)%2)*M, M);

				} else {
					Lat[0] ->propagate(Gg_b,Seg[mon_nr[b]]->G1,s%2,(s+1)%2,M);
					Lat[0] ->AddPhiS(rho+molmon_nr[b]*M, Gg_f+(s+1)*M,Gg_b+((s+1)%2)*M, M);
				}
				s++;
			}
		}
	}

	delete [] GB;
	return success;

}


bool mol_comb::ComputePhi() {
	if (debug) cout <<"ComputePhi for mol_comb " + name << endl;
	bool success=true;	
	int M=Lat[0]->M;
	Real* G0 = new Real[M]; Zero(G0,M);
	Real* Mask = new Real[M]; Zero(Mask,M);
	int gradients=Lat[0]->gradients;

	int MX=Lat[0]->MX;
	int MY=Lat[0]->MY;
	int MZ=Lat[0]->MZ;
	int JX=Lat[0]->JX;
	int JY=Lat[0]->JY;
	int JZ=Lat[0]->JZ;
	bool doit;
	GN=0;
	int x=0,y=0,z=0,i=0,i_old=0;
	vector<int> Ds;
	vector<int> SegPinned;


	if (ring) {
	
		if (IsPinned()) {
			
			int length=MolMonList.size();
			int i=0;
			while (i<length) {
        			if (Seg[MolMonList[i]]->GetFreedom()=="pinned") SegPinned.push_back(MolMonList[i]);
				i++;
			}
			length = mon_nr.size();
			
			int pinnedlength=SegPinned.size();
			for (int j=0; j<pinnedlength; j++) {
				int Ds1=0,Ds2=0;
				bool found=false;
				for (int i=0; i<length; i++) {
					if (found) {
						Ds2 +=n_mon[i];
					} else {
						if (mon_nr[i]==SegPinned[j]) found=true; else Ds1 +=n_mon[i];
					}
				}
				if (Ds1<Ds2) Ds.push_back(Ds1); else Ds.push_back(Ds2); //shortest distance along contour to 'ends' of the chain
				//cout <<"Ds for "<<SegPinned[j] << " : " << Ds[j] << endl; 
			}
		}

		switch (gradients) {
			case 3:
				for (z=1; z<MZ+1; z++)
			case 2:
				for (y=1; y<MY+1; y++)
			case 1:
				for (x=1; x<MX+1; x++) {
					i_old=i; i=JX*x+JY*y+JZ*z;
					Mask[i_old]=0;Mask[i]=1;
					doit =true;
					int pinnedlength=SegPinned.size();
					for (int j=0; j<pinnedlength; j++) {
						if (doit) doit=Seg[SegPinned[j]]->CanBeReached(x,y,z,Ds[j]*Lat[0]->fjc); //make sure that the pinned positions can be reached.
					}
					if (Seg[mon_nr[0]]->G1[i]>0 && doit) {
						Times(G0,Mask,Seg[mon_nr[0]]->G1,M);

						if (Markov==2) success = GoBackAndForth2ndO(Mask,G0); else success = GoBackAndForth(Mask,G0);
					}
	
				}
		}	


	}

	if (Markov==2) success = GoBackAndForth2ndO(Mask,G0); else success = GoBackAndForth(Mask,G0); 
	
	delete [] G0; delete [] Mask;

	return success;
}


