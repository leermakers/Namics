#include "molecule.h"
#include "mol_branched.h"


mol_branched::mol_branched(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {}


mol_branched::~mol_branched() {
}

void mol_branched::BackwardBra(Real* G_start, int generation, int &s){//not yet robust for GPU computations: GS and GX need to be available on GPU
	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=Lat[0]->M;
	//Real* GS = new Real[4*M];
	Real* GS= (Real*) malloc(4*M*sizeof(Real));
	int ss=0;
	for (int k = bN ; k >= b0 ; k--){
		if (k>b0 && k<bN) {
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
						BackwardBra(Gg_b,Br[i],s);
					}
				}
				free(GX);
				k++;
			} else {
				propagate_backward(Seg[mon_nr[k]]->G1,s,k,generation,M);
			}

		} else {
			propagate_backward(Seg[mon_nr[k]]->G1,s,k,generation,M);
		}
	}
	free(GS);
}

Real* mol_branched::ForwardBra(int generation, int &s) {
if (debug) cout <<"ForwardBra in mol_test " << endl;
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
					Gb.push_back(ForwardBra(Gnr[k],s));
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
			Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,generation,M);
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
	int s=0;
	Real* G=ForwardBra(generation,s);
	//Lat[0]->remove_bounds(G);
	GN=Lat[0]->WeightedSum(G);
	s--;
	if (save_memory) {Cp(Gg_b,Seg[mon_nr[last_b[0]]]->G1,M); Cp(Gg_b+M,Seg[mon_nr[last_b[0]]]->G1,M);} //toggle; initialize on both spots the same G1, so that we always get proper start.
	BackwardBra(Seg[mon_nr[last_b[0]]]->G1,generation,s);
if (debug) cout <<" ComputePhiBra in Mol " << endl;
	return success;
}




