#include "molecule.h"
#include "mol_water.h"

//Leermakers, Rabinovich, Balabaev, Phys Rev. E 67, 011910
mol_water::mol_water(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {
phib1=0;
}


mol_water::~mol_water() {
}


void mol_water::AddToGP(Real* GP) {
if (debug) cout <<"AddToGP in mol_water " << endl;
 	int M=Lat[0]->M;
	Real* G=Seg[MolMonList[0]]->G1;

	for (int i=0; i<M; i++) {
		GP[i]+=(phi[i]-phibulk)-phib1*G[i]/(1-Kw*phib1*G[i]) + phib1/(1-Kw*phib1);
	}
}

void mol_water::AddToF(Real* F) {
if (debug) cout <<"AddToGP in mol_water " << endl;
 	int M=Lat[0]->M;
	Real* G=Seg[MolMonList[0]]->G1;

	for (int i=0; i<M; i++) {
		F[i]-=(phi[i])-phib1*G[i]/(1-Kw*phib1*G[i]);	}
		//F[i]-=(phi[i]-phibulk)-phib1*G[i]/(1-Kw*phib1*G[i]) + phib1/(1-Kw*phib1);	}

}




Real mol_water::GetPhib1() {
if (debug) cout <<"GetPhib1 in mol_water " << endl;
	if (phibulk <0) { cout <<"problem in computation of phib1 for MolType water." << endl; return 0;}
	phib1=1/Kw +(1-sqrt(4*Kw*phibulk+1))/(2*Kw*Kw*phibulk);

	return phib1;
}


bool mol_water::ComputePhi() {
if (debug) cout <<"ComputePhi in mol_water " << endl;
	bool success=true;
	int M =Lat[0]->M;
	if (phib1>0) {
		Real* G=Seg[MolMonList[0]]->G1;
		for (int i=0; i<=M; i++) {
			rho[i]=phib1*G[i]/pow((1-Kw*phib1*G[i]),2);
		}
		GN=Lat[0]->ComputeGN(rho,Markov,M)/phib1;
	}
	return success;
}




