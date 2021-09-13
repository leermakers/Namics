#include "molecule.h"
#include "mol_branched.h"


mol_branched::mol_branched(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {

}


mol_branched::~mol_branched() { }

void mol_branched::BackwardBra2ndO(Real* G_start, int generation,int &unity, int &s){//not yet robust for GPU computations: GS and GX need to be available on GPU
if (debug) cout <<"BackwardBra2ndO in mol_branched " << endl;
	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=Lat[0]->M;
	Real* GS= (Real*) malloc(3*M*sizeof(Real));
	Real* GB= (Real*) malloc(2*size*M*sizeof(Real));
	int ss=0;
	for (int k = bN ; k >= b0 ; k--){
		if (Gnr[k]!=generation) {
			Br.clear(); Gb.clear();
			while (Gnr[k] != generation){
				Br.push_back(Gnr[k]);
				if (save_memory) {
					Gb.push_back(Gg_f+last_stored[k]*M*size);
				} else {
					Gb.push_back(Gg_f+last_s[Gnr[k]]*M*size);
				}

				ss=first_s[Gnr[k]];
				k-=(last_b[Gnr[k]]-first_b[Gnr[k]]+1) ;
			}
			Br.push_back(generation); ss--;
			if (save_memory) {
				Gb.push_back(Gg_f+last_stored[k]*M*size);
			} else {
				Gb.push_back(Gg_f+ss*M*size); //backbone last before branch point.

			}
			int length = Br.size();
			Real* GX= (Real*) malloc(length*M*sizeof(Real));

			for (int i=0; i<length; i++) {
				Lat[0]->Terminate(GX+i*M,Gb[i],Markov,M);
			}

			Cp(GB,Gg_b+((s+1)%2)*M*size,M*size); //Upto the branch point; no sides connected

			for (int i=0; i<length; i++) {
				Cp(GS+2*M,UNITY,M);
				for (int j=0; j<length; j++) {
					if (i !=j) {
						if (j==length-1) { //linking main chain
							Cp(Gg_b,Gb[j],M*size);
							Lat[0]->propagateF(Gg_b,UNITY,P,0,1,M); //connect main chain including semiflexibility
							Times(GB+M*size,GB,Gg_b+M*size,M*size);
						} else { //linking sides
							Cp(GS,GX+j*M,M);
							Lat[0]->propagate(GS,UNITY,0,1,M);
							Times(GS+2*M,GS+2*M,GS+M,M);
						}
					}
				}
				if (i<length-1) {
					for (int t=0; t<size; t++) Times(Gg_b + t*M,GB +M*size + t*M,GS+2*M,M); //freely jointed onto branch
					Cp(Gg_b+M*size,Gg_b,M*size);
					unity=-1;
					BackwardBra2ndO(Gg_b,Br[i],unity,s);

				} else { //prepare for main chain propagation
					for (int t=0; t<size; t++) Times(Gg_b + t*M,GB + t*M,GS+2*M,M); //freely jointed onto main chain
					Cp(Gg_b+M*size,Gg_b,M*size);
				}
			}
			free(GX);
			k++;
		} else {
			if (generation==0) {
				if (ring) {
					if (k==bN) {
						Lat[0]->Initiate(Gg_b + (s%2)*M*size,G_start,Markov,M);
						Lat[0]->AddPhiS(rho+molmon_nr[k]*M, Gg_f+s*M*size, Gg_b+(s%2)*M*size,Markov, M);
						s--;
					} else {
						if (k==bN-1) {
							N= n_mon[k];
							for (int b=0; b<N; b++) {
								Lat[0]->propagateB(Gg_b,Seg[mon_nr[k]]->G1,P,(s+1)%2,s%2,M);
								Lat[0]->AddPhiS(rho+molmon_nr[k]*M, Gg_f+s*M*size, Gg_b+(s%2)*M*size, Markov,M);
								s--;
							}
						} else
						if (k>b0) {
							propagate_backward(Seg[mon_nr[k]]->G1,s,k,P,unity,M); unity=0;
						}
					}
				} else {
					propagate_backward(Seg[mon_nr[k]]->G1,s,k,P,unity,M);unity=0;
				}

			} else {
				propagate_backward(Seg[mon_nr[k]]->G1,s,k,P,unity,M); unity=0;
			}
		}
	}
	free(GS); free(GB);
}


Real* mol_branched::ForwardBra2ndO(Real* G0, int generation, int &s) {
if (debug) cout <<"ForwardBra2nd0 in mol_branched " << endl;
	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=Lat[0]->M;
	Real* GS= (Real*) malloc(3*M*sizeof(Real));
	Real* GB= (Real*) malloc(2*size*M*sizeof(Real));

	Real* Glast=NULL;
	for (int k = b0; k<=bN ; ++k) {
		if (b0<k && k<bN) {
			if (Gnr[k]==generation ){
				Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,P,generation,M);
			} else {
				Br.clear(); Gb.clear();
				Cp(GB,Glast,M*size);
				while (Gnr[k] !=generation) { //collect information from branches.
					Br.push_back(Gnr[k]);
					Gb.push_back(ForwardBra2ndO(G0,Gnr[k],s));
					k+=(last_b[Gnr[k]]-first_b[Gnr[k]]+1);
				}
				int length=Br.size();

				Lat[0]->propagateF(GB,Seg[mon_nr[k]]->G1,P,0,1,M); //propagate main chain to branch point; keep semiflexibility

				Cp(GS+2*M,UNITY,M);
				for (int i=0; i<length; i++) {
					Lat[0]->Terminate(GS,Gb[i],Markov,M);
					Lat[0]->propagate(GS,UNITY,0,1,M);
					Times(GS+2*M,GS+2*M,GS+M,M);
				}
				for (int t=0; t<size; t++) Times(GB+M*size+t*M,GB+M*size+t*M,GS+2*M,M); //all side freely jointed
				if (save_memory) {
					Cp(Gs,GB+M*size,M*size); Cp(Gs+M*size,GB+M*size,M*size);
					Cp(Gg_f+(memory[k]-1)*M*size,GB+M*size,M*size); //correct because in this block there is just one segment.
				} else {
					Cp(Gg_f+s*M*size,GB+M*size,M*size);
				}
				s++;
			}
		} else {
			if (k==0 && ring) {
				Lat[0]->Initiate(Gg_f,G0,Markov,M); s++;
			} else {
				Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,P,generation,M);
			}

		}
	}
	free(GS); free(GB);
	return Glast;
}



void mol_branched::BackwardBra(Real* G_start, int generation, int &s){//not yet robust for GPU computations: GS and GX need to be available on GPU
if (debug) cout <<"BackwardBr in mol_branched " << endl;

	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=Lat[0]->M;
	int N;
	Real* GS= (Real*) malloc(4*M*sizeof(Real));
	int ss=0;
	for (int k = bN ; k >= b0 ; k--){
		if (Gnr[k]!=generation) {
			Br.clear(); Gb.clear();
			while (Gnr[k] != generation){
				Br.push_back(Gnr[k]);
				if (save_memory) {
					Gb.push_back(Gg_f+last_stored[k]*M);
				} else {
					Gb.push_back(Gg_f+last_s[Gnr[k]]*M);
				}

				ss=first_s[Gnr[k]];
				k-=(last_b[Gnr[k]]-first_b[Gnr[k]]+1) ;
			}
			Br.push_back(generation); ss--;
			if (save_memory) {
				Gb.push_back(Gg_f+last_stored[k]*M);
			} else {
				Gb.push_back(Gg_f+ss*M);
			}
			int length = Br.size();
			Real* GX= (Real*) malloc(length*M*sizeof(Real));
			for (int i=0; i<length; i++) Cp(GX+i*M,Gb[i],M);
			Cp(GS+3*M,Gg_b+((s+1)%2)*M,M);
			for (int i=0; i<length; i++) {
				Cp(GS+2*M,GS+3*M,M);
				for (int j=0; j<length; j++) {
					if (i !=j) {
						Cp(GS,GX+j*M,M);
						Lat[0]->propagate(GS,UNITY,0,1,M);
						Times(GS+2*M,GS+2*M,GS+M,M);
					}
				}
				Cp(Gg_b,GS+2*M,M);
				Cp(Gg_b+M,GS+2*M,M);
				if (i<length-1) {
					BackwardBra(G_start,Br[i],s);
				}
			}
			free(GX);
			k++;
		} else {
			if (generation==0) {
				if (ring) {
					if (k==bN) {
						Lat[0]->Initiate(Gg_b + (s%2)*M,G_start,Markov,M);
						Lat[0]->AddPhiS(rho+molmon_nr[k]*M, Gg_f+s*M, Gg_b+(s%2)*M,Markov, M);
						s--;
					} else {
						if (k==bN-1) {
							N= n_mon[k];
							for (int b=0; b<N; b++) {
								Lat[0]->propagate(Gg_b,Seg[mon_nr[k]]->G1,(s+1)%2,s%2,M);
								Lat[0]->AddPhiS(rho+molmon_nr[k]*M, Gg_f+s*M, Gg_b+(s%2)*M, Markov,M);
								s--;
							}
						} else
						if (k>b0) propagate_backward(Seg[mon_nr[k]]->G1,s,k,generation,M);
					}
				} else propagate_backward(Seg[mon_nr[k]]->G1,s,k,generation,M);

			} else
				propagate_backward(Seg[mon_nr[k]]->G1,s,k,generation,M);
		}

	}
	free(GS);
}

Real* mol_branched::ForwardBra(Real* G0, int generation, int &s) {
if (debug) cout <<"ForwardBra in mol_branched " << endl;
	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=Lat[0]->M;
	Real* GS= (Real*) malloc(3*M*sizeof(Real));

	Real* Glast=NULL;
	for (int k = b0; k<=bN ; ++k) {
		if (b0<k && k<bN) {
			if (Gnr[k]==generation ){
				Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,generation,M);
			} else {
				Br.clear(); Gb.clear();
				Cp(GS,Glast,M);
				while (Gnr[k] !=generation) {
					Br.push_back(Gnr[k]);
					Gb.push_back(ForwardBra(G0,Gnr[k],s));
					k+=(last_b[Gnr[k]]-first_b[Gnr[k]]+1);
				}
				int length=Br.size();
				Lat[0]->propagate(GS,Seg[mon_nr[k]]->G1,0,2,M);

				for (int i=0; i<length; i++) {
					Cp(GS,Gb[i],M);
					Lat[0]->propagate(GS,UNITY,0,1,M);
					Times(GS+2*M,GS+2*M,GS+M,M);
				}
				if (save_memory) {
					Cp(Gs,GS+2*M,M); Cp(Gs+M,GS+2*M,M);
					Cp(Gg_f+(memory[k]-1)*M,GS+2*M,M); //correct because in this block there is just one segment.
				} else {
					Cp(Gg_f+s*M,GS+2*M,M);
				}
				s++;
			}
		} else {
			if (k==0 && ring) {
				Cp(Gg_f,G0,M); s++;
			} else {
				Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,generation,M);
			}

		}
	}
	free(GS);
	return Glast;
}



bool mol_branched::ComputePhi() {
if (debug) cout <<"ComputePhi in mol_branched " << endl;
	int M=Lat[0]->M;
	bool success=true;
	int generation=0;
	int unity=0;
	int s=0;

	int gradients=Lat[0]->gradients;
	int MX=Lat[0]->MX;
	int MY=Lat[0]->MY;
	int MZ=Lat[0]->MZ;
	int JX=Lat[0]->JX;
	int JY=Lat[0]->JY;
	int JZ=Lat[0]->JZ;
	int x=0,y=0,z=0,i=0,i_old=0;
	vector<int> Ds;
	vector<int> SegPinned;
	Real* G0;
	Real* Mask;
	bool doit;

	Real* G;
	if (ring) {
		G0 = new Real[M]; Zero(G0,M);
		Mask = new Real[M]; Zero(Mask,M);
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
						if (doit) doit=Seg[SegPinned[j]]->CanBeReached(x,y,z,Ds[j]*Lat[0]->fjc); //make sure that the pinned positions can be reached.
					}
					if (Seg[mon_nr[0]]->G1[i]>0 && doit) {
						Times(G0,Mask,Seg[mon_nr[0]]->G1,M);

						if (Markov == 2) {
							G=ForwardBra2ndO(G0,generation,s);
						} else {
							G=ForwardBra(G0,generation,s);
						}
						for (int k=0; k<size; k++) Times(G+k*M,G+k*M,Mask,M); //to make sure GN is computed correctly.
						GN+=Lat[0]->ComputeGN(G,Markov,M);
						s--;
						if (save_memory) {
							Lat[0]->Initiate(Gg_b,G0,Markov,M);
							Lat[0]->Initiate(Gg_b+M*size,G0,Markov,M);
						} //toggle; initialize on both spots the same G1, so that we always get proper start.
						if (Markov == 2) {
							BackwardBra2ndO(G0,generation,unity,s);
						} else {
							BackwardBra(G0,generation,s);
						}
					}
				}
		}
		delete [] G0; delete [] Mask;
	} else {
		if (Markov == 2) {
			G=ForwardBra2ndO(Seg[mon_nr[last_b[0]]]->G1,generation,s);
		} else {
			G=ForwardBra(Seg[mon_nr[last_b[0]]]->G1,generation,s);
		}
		GN=Lat[0]->ComputeGN(G,Markov,M);
		s--;
		if (save_memory) {
			Lat[0]->Initiate(Gg_b,Seg[mon_nr[last_b[0]]]->G1,Markov,M);
			Lat[0]->Initiate(Gg_b+M*size,Seg[mon_nr[last_b[0]]]->G1,Markov,M);
		} //toggle; initialize on both spots the same G1, so that we always get proper start.
		if (Markov == 2) {
			BackwardBra2ndO(Seg[mon_nr[last_b[0]]]->G1,generation,unity,s);
		} else {
			BackwardBra(Seg[mon_nr[last_b[0]]]->G1,generation,s);
		}
	}


	return success;
}



