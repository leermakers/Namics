#include "molecule.h"
#include "mol_linear.h"


mol_linear::mol_linear(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {}


mol_linear::~mol_linear() {
}



bool mol_linear::ComputePhi() {
if (debug) cout <<"ComputePhi in mol_linear " << endl;
	int b0 = first_b[0];
	int bN = last_b[0];
	int M=lat->M;
	bool success=true;
	int unity=0;
	int s=0;
	int N;
#ifdef CUDA
	Real* Mask = (Real*)AllManagedOnDev(M); Zero(Mask,M);
#else
	Real* Mask = new Real[M]; Zero(Mask,M);
#endif
	int gradients=lat->gradients;
	int MX=lat->MX;
	int MY=lat->MY;
	int MZ=lat->MZ;
	int JX=lat->JX;
	int JY=lat->JY;
	int JZ=lat->JZ;
	int x=0,y=0,z=0,i=0,i_old=0;
	vector<int> Ds;
	vector<int> SegPinned;
	Real* Glast=NULL;
	bool doit;

	if (ring) {
#ifdef CUDA
		Real* G0 = (Real*)AllManagedOnDev(M); Zero(G0,M);
#else
		Real* G0 = new Real[M]; Zero(G0,M);
#endif

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
		GN=0;
		switch (gradients) {
			case 3:
				for (z=1; z<MZ+1; z++)
			case 2:
				for (y=1; y<MY+1; y++)
			case 1:
				for (x=1; x<MX+1; x++) {
				i_old=i; i=JX*x+JY*y+JZ*z;
				s=0;
				Mask[i_old]=0;Mask[i]=1;
				doit =true;
				int pinnedlength=SegPinned.size();
				for (int j=0; j<pinnedlength; j++) {
					if (doit) doit=Seg[SegPinned[j]]->CanBeReached(x,y,z,Ds[j]*lat->fjc); //make sure that the pinned positions can be reached.
				}
				if (Seg[mon_nr[0]]->G1[i]>0 && doit) {
					//cout <<" x = " << x << endl;

					Times(G0,Mask,Seg[mon_nr[0]]->G1,M);
					lat->Initiate(Gg_f,G0,Markov,M); //initialisatie
					s++;
					if (Markov==2)
						for (int b = b0+1; b<=bN ; ++b) Glast=propagate_forward(Seg[mon_nr[b]]->G1,s,b,P,0,M);
					else
						for (int b = b0+1; b<=bN ; ++b) Glast=propagate_forward(Seg[mon_nr[b]]->G1,s,b,0,M);

					for (int k=0; k<size; k++) Times(Glast+k*M,Glast+k*M,Mask,M);
					GN+=lat->ComputeGN(Glast,Markov,M);
  					s=chainlength;
					lat->Initiate(Gg_b+(s%2)*M*size,G0,Markov,M);
					lat->AddPhiS(rho+molmon_nr[0]*M,Gg_b+(s%2)*M*size,Gg_f+s*M*size,Markov,M);
					s--;
					N=n_mon[bN-1];
					for (int k=0; k<N; k++) { //first block done here because of initialisation in propagate
						if (Markov==2)
							lat->propagateB(Gg_b,Seg[mon_nr[bN-1]]->G1,P,(s+1)%2,s%2,M);
						else
							lat->propagate(Gg_b,Seg[mon_nr[bN-1]]->G1,(s+1)%2,s%2,M);
						lat->AddPhiS(rho+molmon_nr[bN-1]*M, Gg_f+s*M*size, Gg_b+(s%2)*M*size, Markov,M);
						s--;
					}
					if (Markov==2)
						for (int b = bN-2; b > b0 ; b--) propagate_backward(Seg[mon_nr[b]]->G1,s,b,P,unity,M);
					else
						for (int b = bN-2; b > b0 ; b--) propagate_backward(Seg[mon_nr[b]]->G1,s,b,0,M);
				}
			}
		}
#ifdef CUDA
		cudaFree(G0);
#else
		delete [] G0;
#endif
	} else {

		if (Markov ==2)
			for (int b = b0; b<=bN ; ++b) Glast=propagate_forward(Seg[mon_nr[b]]->G1,s,b,P,0,M);
		else
			for (int b = b0; b<=bN ; ++b) Glast=propagate_forward(Seg[mon_nr[b]]->G1,s,b,0,M);

		GN=lat->ComputeGN(Glast,Markov,M);

		s--;
		if (save_memory) {
			lat->Initiate(Gg_b,Seg[mon_nr[last_b[0]]]->G1,Markov, M);
			lat->Initiate(Gg_b+M*size,Seg[mon_nr[last_b[0]]]->G1,Markov,M);
		}
		if (Markov==2)
			for (int b = bN ; b >= b0 ; b--) propagate_backward(Seg[mon_nr[b]]->G1,s,b,P,unity,M);
		else
			for (int b = bN ; b >= b0 ; b--) propagate_backward(Seg[mon_nr[b]]->G1,s,b,0,M);
	}

#ifdef CUDA
	cudaFree(Mask);
#else
	delete [] Mask;
#endif
	return success;
}




