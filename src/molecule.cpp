#include "molecule.h"


Molecule::Molecule(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) {
	In=In_; Seg=Seg_; name=name_;  Lat=Lat_;
if (debug) cout <<"Constructor for Mol " + name << endl;
	KEYS.push_back("freedom");
	KEYS.push_back("composition");
	KEYS.push_back("ring");
	KEYS.push_back("theta");
	KEYS.push_back("phibulk");
	KEYS.push_back("n");
	KEYS.push_back("save_memory");
	KEYS.push_back("restricted_range");
	KEYS.push_back("compute_width_interface");
	KEYS.push_back("Kw");
	KEYS.push_back("Markov");
	KEYS.push_back("k_stiff");
	KEYS.push_back("phi_LB_x");
	KEYS.push_back("phi_UB_x");
	KEYS.push_back("B");

	width=0;
	phi_LB_X=0.0;
	phi_UB_X=0.0;
	phi1=0;
	phiM=0;
	Dphi=0;
	pos_interface=0;
	ring=false;
	all_molecule=false;
	Markov =1;
	Var_scan_value=-1;
	Var_search_value=-1;
	Var_target=-1;
	FillRangesList.clear();
	Filling=false;
	save_memory=false;
	J=0; Delta_MU=0; B=1;

}

Molecule::~Molecule() {
	DeAllocateMemory();
}

void Molecule :: DeAllocateMemory(){
if (debug) cout <<"DeallocateMemory for Mol " + name << endl;
	if (!all_molecule) return;
	//if (H_phi!=NULL)
	free(H_phi);
	free(H_phitot);
	if (Markov==2) { free(P);}
	if (freedom=="clamped") {
		free(H_Bx);
		free(H_By);
		free(H_Bz);
		free(H_Px1);
		free(H_Py1);
		free(H_Pz1);
		free(H_Px2);
		free(H_Py2);
		free(H_Pz2);
		free(H_mask1);
		free(H_mask2);
		free(H_gn);
	}
#ifdef CUDA
	if (freedom=="clamped") {
		cudaFree(Bx);
		cudaFree(By);
		cudaFree(Bz);
		cudaFree(Px1);
		cudaFree(Py1);
		cudaFree(Pz1);
		cudaFree(Px2);
		cudaFree(Py2);
		cudaFree(Pz2);
		cudaFree(mask1);
		cudaFree(mask2);
		cudaFree(gn);
		cudaFree(Gg_f);
		cudaFree(Gg_b);
		cudaFree(g1);
		cudaFree(rho);
		cudaFree(phi);
	} else {
		cudaFree(Gg_f);
		cudaFree(Gg_b);
		cudaFree(phi);
	}
	cudaFree(phitot);
	if (save_memory) cudaFree(Gs);
#else
	if (freedom=="clamped"){
		free(rho);
		free(g1);
	}
	free(Gg_f);
	free(Gg_b);
	if (save_memory) free(Gs);
	free(UNITY);

#endif
	all_molecule=false;
}

bool Molecule::DeleteAl() {
	bool success = true;
	int length_al = MolAlList.size();
	if (length_al>0) {
		for (int k=0; k<length_al; k++) delete Al[k];
		Al.clear();
	}
	Gnr.clear();
	first_s.clear();
	last_s.clear();
	first_b.clear();
	last_b.clear();
	mon_nr.clear();
	n_mon.clear();
	molmon_nr.clear();
	memory.clear();
	last_stored.clear();
	MolAlList.clear();
	MolMonList.clear();
	return success;
}

void Molecule:: AllocateMemory() {
if (debug) cout <<"AllocateMemory in Mol " + name << endl;
	DeAllocateMemory();
	int M=Lat[0]->M;
	int m=0;
	if (Markov==2){// && Lat[0]->lattice_type == simple_cubic) {
		int FJC = Lat[0]->FJC;
		P = (Real*) malloc(FJC*sizeof(Real));
		Real Q=0;
		KStiff=k_stiff;
		for (int k=0; k<FJC-1; k++) {
			P[k]=exp(-0.5*KStiff*(k*PIE/(FJC-1))*(k*PIE/(FJC-1)) ); //alternative to put U(theta)=-k(1-cos(theta))
			if (k>0) {
				if (Lat[0]->lattice_type==hexagonal) Q+= 2*P[k]; else Q+= 4*P[k]; //alternative is to use u_bend = Kstiff(1-cos(theta)), persistence length is l_p = b/ln <cos (theta)>
			} else Q=P[k];
		}
		P[FJC-1]=0; //Q+=P[FJC-1];
		//if (Lat[0]->lattice_type==hexagonal) Q=2*Q-P[0]; else Q=4*Q-3*P[0];
		if (Lat[0]->lattice_type==hexagonal&& !Lat[0]->stencil_full) Q*=2.0;
		for (int k=0; k<FJC; k++) { P[k]/=Q;
			cout << "P["<<k<<"] = " << P[k] << endl;
		}
	}

/*
	if (Markov==2 && Lat[0]->lattice_type == hexagonal) {
		int FJC = Lat[0]->FJC;
		P = (Real*) malloc(2*sizeof(Real)); //assuming only default k_stiff value for P's so that P array is small.
		Real Q=0;
		KStiff=k_stiff;
		P[0]=exp(-0.5*KStiff*(PIE/3.0)*(PIE/3.0) ); Q+= 3*P[0]; //alternative is to use P[0]=exp(-KStiff*(1-cos(theta))
		P[1]=exp(-0.5*KStiff*(2.0*PIE/3.0)*(2.0*PIE/3.0) ); Q+= 6*P[1];
		P[FJC-1]=0;
		for (int k=0; k<Lat[0]->FJC-1; k++) { P[k]/=Q;
			cout << "P["<<k<<"] = " << P[k] << endl;
		}
	}
*/

	if (freedom=="clamped") {
		m=Lat[0]->m[Seg[mon_nr[0]]->clamp_nr];
		n_box = Seg[mon_nr[0]]->n_box;
	}

	if (save_memory) {
		int length_ = mon_nr.size();
		for (int i=0; i<length_; i++) last_stored.push_back(0);
		for (int i=0; i<length_; i++) {
			int n=int(pow(n_mon[i]*2,1.0/2.0)+0.5); //in sfbox apparently the pow 1/3 is reached. Here it fails for unknown reasons.
			//if (n_mon[i]<120) n++;
			//if (n_mon[i]<60) n++;
			if (n>n_mon[i]) n=n_mon[i]; //This seems to me to be enough...needs a check though..
			if (i==0) memory.push_back(n); else memory.push_back(n+memory[i-1]);
			//cout <<"n " << n << endl;
		}
	}
	N=0;
	if (save_memory) {
		N=memory[n_mon.size()-1];
	} else {
		int length_ = mon_nr.size();
		for (int i=0; i<length_; i++) {N+=n_mon[i];}
	}

	H_phi = (Real*) malloc(M*MolMonList.size()*sizeof(Real)); H_Zero(H_phi,M*MolMonList.size());
	H_phitot = (Real*) malloc(M*sizeof(Real)); H_Zero(H_phitot,M);
	if (freedom=="clamped") {
		H_Bx=(int*) malloc(n_box*sizeof(int));
		H_By=(int*) malloc(n_box*sizeof(int));
		H_Bz=(int*) malloc(n_box*sizeof(int));
		H_Px1=(int*) malloc(n_box*sizeof(int));
		H_Py1=(int*) malloc(n_box*sizeof(int));
		H_Pz1=(int*) malloc(n_box*sizeof(int));
		H_Px2=(int*) malloc(n_box*sizeof(int));
		H_Py2=(int*) malloc(n_box*sizeof(int));
		H_Pz2=(int*) malloc(n_box*sizeof(int));
		H_mask1=(Real*) malloc(n_box*m*sizeof(Real)); H_Zero(H_mask1,m*n_box);
		H_mask2=(Real*) malloc(n_box*m*sizeof(Real)); H_Zero(H_mask2,m*n_box);
		H_gn = (Real*) malloc(n_box*sizeof(Real));
	}
#ifdef CUDA
	if (freedom=="clamped") {
		Bx=(int*)AllIntOnDev(n_box);
		By=(int*)AllIntOnDev(n_box);
		Bz=(int*)AllIntOnDev(n_box);
		Px1=(int*)AllIntOnDev(n_box);
		Py1=(int*)AllIntOnDev(n_box);
		Pz1=(int*)AllIntOnDev(n_box);
		Px2=(int*)AllIntOnDev(n_box);
		Py2=(int*)AllIntOnDev(n_box);
		Pz2=(int*)AllIntOnDev(n_box);
		mask1=(Real*)AllOnDev(n_box*m);
		mask2=(Real*)AllOnDev(n_box*m);
		gn =(Real*)AllOnDev(n_box);
		Gg_f=(Real*)AllOnDev(m*n_box*N); Zero(Gg_f,m*n_box*N);
		Gg_b=(Real*)AllOnDev(m*n_box*2); Zero(Gg_b,m*n_box*2);
		g1=(Real*)AllOnDev(m*n_box); Zero(g1,m*n_box);
		rho=(Real*)AllOnDev(m*n_box*MolMonList.size()); Zero(rho,m*n_box*MolMonList.size());
		phi=(Real*)AllOnDev(M*MolMonList.size()); Zero(phi,M*MolMonList.size());
		if (save_memory) {Gs=(Real*)AllOnDev(m*n_box*2); Zero(Gs,m*n_box*2);}
	} else {
		Gg_f=(Real*)AllOnDev(M*N); //not yet for Markov =2
		Gg_b=(Real*)AllOnDev(M*2);
		phi=(Real*)AllOnDev(M*MolMonList.size());
		rho=phi;
		if (save_memory) Gs =(Real*)AllOnDev(M*2);
	}
	phitot=(Real*)AllOnDev(M);
	UNITY = (Real*)AllOnDev(M);
#else
	if (freedom=="clamped") {
		gn=H_gn;
		mask1=H_mask1; mask2=H_mask2;
		Bx=H_Bx; By=H_By; Bz=H_Bz;
		Px1=H_Px1; Py1=H_Py1; Pz1=H_Pz1;
		Px2=H_Px2; Py2=H_Py2; Pz2=H_Pz2;
		Gg_f = (Real*) malloc(m*N*n_box*sizeof(Real)); Zero(Gg_f,m*N*n_box);
		Gg_b = (Real*) malloc(m*2*n_box*sizeof(Real)); Zero(Gg_b,m*2*n_box);

		g1=(Real*) malloc(m*n_box*sizeof(Real)); Zero(g1,m*n_box);
		rho=(Real*)malloc(m*n_box*MolMonList.size()*sizeof(Real)); Zero(rho,m*n_box*MolMonList.size());
		if (save_memory) {Gs=(Real*) malloc(m*n_box*2*sizeof(Real));Zero(Gs,m*n_box*2);}
		phi=H_phi;
	} else {
		Gg_f = (Real*) malloc(M*N*sizeof(Real)*size);
		Gg_b = (Real*) malloc(M*2*sizeof(Real)*size);
		Zero(Gg_f,M*N*size);
		Zero(Gg_b,2*M*size);
		phi=H_phi;
		rho=phi;
		if (save_memory) {Gs=(Real*) malloc(2*M*sizeof(Real)*size); Zero(Gs,2*M*size);}
	}
	phitot = H_phitot;
	UNITY = (Real*) malloc(M*sizeof(Real)*size); Zero(UNITY,M*size);
#endif

	int length =MolAlList.size();
	if (freedom!="clamped") Seg[mon_nr[0]]->clamp_nr=0;
	for (int i=0; i<length; i++) Al[i]->AllocateMemory(Seg[mon_nr[0]]->clamp_nr,n_box);
	all_molecule=true;
}

bool Molecule:: PrepareForCalculations(int *KSAM) {
if (debug) cout <<"PrepareForCalculations in Mol " + name << endl;
int m=0;
if (freedom=="clamped") m=Lat[0]->m[Seg[mon_nr[0]]->clamp_nr];
int M=Lat[0]->M;
	if (freedom=="clamped") {
		std::fill(H_mask1,H_mask1+n_box*m,0);
		std::fill(H_mask2,H_mask2+n_box*m,0);
		int jx=Lat[0]->jx[Seg[mon_nr[0]]->clamp_nr];
		int jy=Lat[0]->jy[Seg[mon_nr[0]]->clamp_nr];
		int m=Lat[0]->m[Seg[mon_nr[0]]->clamp_nr];
		for (int i=0; i<n_box; i++) {
			H_Bx[i]=Seg[mon_nr[0]]->bx[i];
			H_By[i]=Seg[mon_nr[0]]->by[i];
			H_Bz[i]=Seg[mon_nr[0]]->bz[i];
			H_Px1[i]=Seg[mon_nr[0]]->px1[i];
			H_Py1[i]=Seg[mon_nr[0]]->py1[i];
			H_Pz1[i]=Seg[mon_nr[0]]->pz1[i];
			H_Px2[i]=Seg[mon_nr[0]]->px2[i];
			H_Py2[i]=Seg[mon_nr[0]]->py2[i];
			H_Pz2[i]=Seg[mon_nr[0]]->pz2[i];
			H_mask1[i*m + jx*(H_Px1[i]-H_Bx[i])+jy*(H_Py1[i]-H_By[i])+(H_Pz1[i]-H_Bz[i])]=1;
			H_mask2[i*m + jx*(H_Px2[i]-H_Bx[i])+jy*(H_Py2[i]-H_By[i])+(H_Pz2[i]-H_Bz[i])]=1;
		}
#ifdef CUDA
		TransferDataToDevice(H_mask1,mask1,m*n_box);
		TransferDataToDevice(H_mask2,mask2,m*n_box);
		TransferDataToDevice(H_Bx,Bx,n_box);
		TransferDataToDevice(H_By,By,n_box);
		TransferDataToDevice(H_Bz,Bz,n_box);
		TransferDataToDevice(H_Px1,Px1,n_box);
		TransferDataToDevice(H_Py1,Py1,n_box);
		TransferDataToDevice(H_Px2,Px2,n_box);
		TransferDataToDevice(H_Pz1,Pz1,n_box);
		TransferDataToDevice(H_Py2,Py2,n_box);
		TransferDataToDevice(H_Pz2,Pz2,n_box);
	//	TransferDataToDevice(H_u,u, (int)MolMonList.size()*M);
#endif
	}
	Cp(UNITY,KSAM,M);
	int length_al=MolAlList.size();
	for (int i=0; i<length_al; i++) Al[i]->PrepareForCalculations();
	bool success=true;
	Zero(phitot,M);

	//int length = MolMonList.size();
	//int i=0;
	//while (i<length) {

		//if (Seg[MolMonList[i]]->freedom=="tagged" || Seg[MolMonList[i]]->freedom=="clamp" ) Zero(u+i*M,M);
		//Lat[0]->set_bounds(u+i*M);
		//Boltzmann(G1+i*M , u+i*M, M); */

		//if (Seg[MolMonList[i]]->freedom=="tagged" || Seg[MolMonList[i]]->freedom=="clamp" ) Zero(u+i*M,M);
		//Lat[0]->set_bounds(u+i*M);
		//Cp(G1+i*M,Seg[MolMonList[i]]->G1,M);
		//Boltzmann(G1+i*M , u+i*M, M);


		//if (Seg[MolMonList[i]]->freedom=="pinned") Times(G1+i*M,G1+i*M,Seg[MolMonList[i]]->MASK,M);
		//if (Seg[MolMonList[i]]->freedom=="tagged") Cp(G1+i*M,Seg[MolMonList[i]]->MASK,M);
		//Lat[0]->set_bounds(G1+i*M);
		//if (!(Seg[MolMonList[i]]->freedom ==" frozen" || Seg[MolMonList[i]]->freedom =="tagged")) Times(G1+i*M,G1+i*M,KSAM,M);
		//i++;
	//}
	Zero(phi,M*MolMonList.size());
	compute_phi_alias=false;
	if (freedom=="clamped") {
		int chainlength_even=chainlength%2;
		int pathlength_even;
		for (int i=0; i<n_box; i++) {
			pathlength_even=0;
			pathlength_even=(H_Px2[i]-H_Px1[i]+H_Py2[i]-H_Py1[i]+H_Pz2[i]-H_Pz1[i])%2;
			if (chainlength_even == pathlength_even)
			cout <<" Warning, for chain part " << i << " the paths between clamps is not commensurate with the length of the chain fragment. Consider moving one of the calmp point by one (more) site!" << endl;
		}
		n = n_box;
		theta=n_box*chainlength;
	}


	return success;
}

bool Molecule::CheckInput(int start_, bool checking) {
if (debug) cout <<"Molecule:: CheckInput for mol " << name << endl;
start=start_;
phibulk=0;
n=0;
theta=0;
norm=0;
var_al_nr=-1;
if (debug) cout <<"CheckInput for Mol " + name << endl;
	bool success=true;
	if (!In[0]->CheckParameters("mol",name,start,KEYS,PARAMETERS,VALUES)) {
		success=false;
	} else {
		save_memory=false;
		if (GetValue("save_memory").size()>0) {
			save_memory=In[0]->Get_bool(GetValue("save_memory"),false);
		}
		if (GetValue("composition").size()==0) {cout << "For mol '" + name + "' the definition of 'composition' is required" << endl; success = false;
		} else {
			try {
				// The chain of success bools seems to be broken and handled incorrectly..
				// This is the result of calling the same function (interpret(), probably) multiple times with different results,
				// causing only the last success to transfer correctly.
				// Throw - catch for safety. Thrown by Interpret() and Decomposition() itself;
			if (!Decomposition(GetValue("composition")))
				// This error message is not even displayed.
				{cout << "For mol '" + name + "' the composition is rejected. " << endl; success=false;}
			} catch (const char* error) {
				cerr << error << endl;
				success = false;
			}
		}
		if (GetValue("restricted_range").size()>0) {
			if (GetValue("freedom")!="range_restricted") cout <<"For mol '" + name + "' freedom is not set to 'range_restricted' and therefore  the value of 'restricted_range' is ignored" << endl;
		}

		if (checking) {
			if (GetValue("freedom").size() > 0) freedom = GetValue("freedom"); else {
				cout <<"For molecule " << name << " no value for 'freedom' was found " << endl;
				success=false;
			}
			return success; //because moltype and freedom are known; as start <0 the checkinput can be terminated.
		}
		if (IsPinned()) {
			if (GetValue("freedom").size()==0) {
					cout <<"For mol " + name + " the setting for 'freedom' was not set" << endl; return false;
			} else {
				vector<string> free_list;
				free_list.push_back("restricted");
				free_list.push_back("fill_range");
				if (!In[0]->Get_string(GetValue("freedom"),freedom,free_list,"In mol " + name + " the value for 'freedom' is not recognised ")) return false;
				if (freedom=="restricted") {
					if (GetValue("theta").size() ==0 && GetValue("n").size()==0) {
							cout <<"In mol " + name + ", the setting 'freedom = restricted' or 'freedom = range_restricted',should be combined with a value for 'theta' or 'n'; do not use both settings! "<<endl; success=false;
					} else {
							if (GetValue("theta").size() >0 && GetValue("n").size()>0) {
							cout <<"In mol " + name + ", the setting 'freedom = restricted' of 'freedom = range_restricted' do not specify both 'n' and 'theta' "<<endl; success=false;
					} else {
							if (GetValue("n").size()>0) {n=In[0]->Get_Real(GetValue("n"),10*Lat[0]->volume);theta=n*chainlength;}
							if (GetValue("theta").size()>0) {theta = In[0]->Get_Real(GetValue("theta"),10*Lat[0]->volume);n=theta/chainlength;}
							if (theta < 0 || theta > Lat[0]->volume) {
								cout << "In mol " + name + ", the value of 'n' or 'theta' is out of range 0 .. 'volume', cq 'volume'/N." << endl; success=false;
							}
						}
					}
				}
				if (freedom=="fill_range") {
					freedom="restricted";
					Filling=true;
					if (GetValue("theta").size() > 0 || GetValue("n").size()>0 || GetValue("phibulk").size() > 0) {
						if (start==1) cout <<"For mol " + name + " the freedom is set to 'fill-range-of{-mon_name}' and therefore the value of 'theta', the value of 'n', or the value of 'phibulk' is ignored. " << endl;
					}
				}
			}
		} else
		if ( IsClamped() ){
			freedom="clamped";
			if (GetValue("freedom").size() > 0) freedom = In[0]->Get_string(GetValue("freedom"),"clamped");
			if (freedom !="clamped") {
				cout <<"For mol " + name + " the setting for 'freedom' was not equal to 'clamped'; This is not consistent with composition. " << endl;
				success=false;
			} else {
				n_box=Seg[mon_nr[0]]->n_box;
				//int m=Lat[0]->m[Seg[mon_nr[0]]->clamp_nr];
				if (GetValue("theta").size() >0) {
					theta=In[0]->Get_Real(GetValue("theta"),n_box*chainlength);
					if (theta!=n_box*chainlength)
					cout <<"Input value for 'theta' is ignored for clamped molecule '" + name + "'"  << endl;
				} theta= n_box*chainlength;
				if (GetValue("n").size() >0) { n=In[0]->Get_Real(GetValue("n"),n);
					if (n!=n_box)
					cout <<"Input value for 'n' is ignored for clamped molecule '" + name + "'"  << endl;
				} n= n_box;
				if (GetValue("phibulk").size() >0) {
					success=false;
					cout <<"For mol '" + name + "' the value of freedom is 'clamped' and therfore you can not set 'phibulk'. Problem terminated. " << endl;
				}
			}
		} else
		if (GetValue("freedom").size()==0 && !IsTagged() ) {
			cout <<"For mol " + name + " the setting 'freedom' is expected: options: 'free' 'restricted' 'solvent' 'neutralizer' 'range_restricted' 'clamped' 'tagged' . Problem terminated " << endl; success = false;
			} else {

				if (!IsTagged()) {
				vector<string> free_list;
				if (!IsPinned()) {
					free_list.push_back("free");
					free_list.push_back("solvent");
					free_list.push_back("neutralizer");
					free_list.push_back("range_restricted");
					free_list.push_back("gradient");
				}
				free_list.push_back("restricted");
				if (!In[0]->Get_string(GetValue("freedom"),freedom,free_list,"In mol " + name + " the value for 'freedom' is not recognised ")) success=false;
				if (freedom == "solvent") {
					if (IsPinned()) {success=false; cout << "Mol '" + name + "' is 'pinned' and therefore this molecule can not be the solvent" << endl; }
				}
				if (MolType == water && freedom!="solvent" ) {
					cout <<"MolType 'water' can only be used for the component with freedom 'solvent'. Job terminated." << endl; success=false;
				}
				if (MolType == water && freedom=="solvent" ) {
					Kw=50;
					if (GetValue("Kw").size()>0) {
						Kw=In[0]->Get_Real(GetValue("Kw"),-1);
						if (Kw < 0) {
								cout << "Value for assosication constant Kw (used in the water model) should be positive." << endl;
								cout << "In mol " + name + ", the value of 'Kw' serves in the 'association' water model. " << endl;
								cout << "A value of  K~50 is advised " << endl;
								cout << "X + X --Kw--> X_2; X_2 + X --Kw--> X_3; etc." << endl;
								cout << "See Phys REv E 67, 011910 (2003) for more information. " << endl;
								success=false;
							}
					} else {
						cout << "A default value for the value of 'Kw = 50' is used, because I failed to find an input for this quantity " << endl;
					}
				}
				if (freedom == "neutralizer") {
					if (IsPinned()) {success=false; cout << "Mol '" + name + "' is 'pinned' and therefore this molecule can not be the neutralizer" << endl; }
					if (!IsCharged()) {success=false; cout << "Mol '" + name + "' is not 'charged' and therefore this molecule can not be the neutralizer" << endl; }
					if (IsClamped()) {success=false; cout <<"Mol '" + name + "' is 'clamped' and therefore this molecule can not be the neutralizer" << endl;}
				}
				if (freedom == "free") {
					if (GetValue("phibulk").size() ==0) {
						cout <<"In mol " + name + ", the setting 'freedom = free' should be combined with a value for 'phibulk'. "<<endl; return false;
					} else {
						phibulk=In[0]->Get_Real(GetValue("phibulk"),-1);
						if (phibulk < 0 || phibulk >1) {
							cout << "In mol " + name + ", the value of 'phibulk' is out of range 0 .. 1." << endl; return false;
						}
					}
				}

				//if (freedom == "gradient") { //for the time being only in one-gradient systems; this is tested in system.
				B=1;
				if (GetValue("B").size()>0){
					B=In[0]->Get_Real(GetValue("B"),B);
					if (B<1e-9) {
						cout <<"for Mol" + name + " mobility B should have a posititve value. Default value B=1 is chosen. " << endl;
						B=1;
					}
				}

				if (freedom == "gradient") { //for the time being only in one-gradient systems; this is tested in system.
					//if (GetValue("phibulk").size() >0) {
					//	cout <<"In mol " + name + ", the setting 'freedom : gradient' should not be combined with a value for 'phibulk' but with values for 'phi_LB_x' and 'phi_UP_x' "<<endl; return false;
					//} else {
						if (GetValue("phi_LB_x").size()==0 || GetValue("phi_UB_x").size()==0) {
							cout <<"in mol " + name + "the setting 'freedom : gradient' should be combined with values of 'phi_LB_x' and 'phi_UB_x'=phibulk " << endl; return false;
						} else {
							phi_UB_X=In[0]->Get_Real(GetValue("phi_UB_x"),-1); phibulk=phi_UB_X;
							phi_LB_X=In[0]->Get_Real(GetValue("phi_LB_x"),-1);
							if (phi_UB_X < 0 || phi_UB_X >1 || phi_LB_X <0 || phi_UB_X > 1 ) {
								cout << "In mol " + name + ", the value of 'phi_UB_x' or 'phi_LB_x' is out of range 0 .. 1." << endl; return false;
							}
						}
					//}
				}

				if (freedom == "restricted" || freedom=="range_restricted") {
					//if (GetValue("phibulk").size()>0) {
					//	cout << "In mol " + name + ", the setting of 'freedom = restricted' or 'freedom = range_restricted', can not not be combined with 'phibulk'  use 'theta' or 'n'  instead." << endl; success=false;
					//} else {
						if (GetValue("theta").size() ==0 && GetValue("n").size()==0) {
							cout <<"In mol " + name + ", the setting 'freedom = restricted' or 'freedom = range_restricted',should be combined with a value for 'theta' or 'n'; do not use both settings! "<<endl; success=false;
						} else {
							if (GetValue("theta").size() >0 && GetValue("n").size()>0) {
							cout <<"In mol " + name + ", the setting 'freedom = restricted' of 'freedom = range_restricted' do not specify both 'n' and 'theta' "<<endl; success=false;
							} else {

								if (GetValue("n").size()>0) {n=In[0]->Get_Real(GetValue("n"),10*Lat[0]->volume);theta=n*chainlength;}
								if (GetValue("theta").size()>0) {theta = In[0]->Get_Real(GetValue("theta"),10*Lat[0]->volume);n=theta/chainlength;}
								if (theta < 0 || theta > Lat[0]->volume) {
									cout << "In mol " + name + ", the value of 'n' or 'theta' is out of range 0 .. 'volume', cq 'volume'/N." << endl; success=false;

								}
							}
						}
					//}
				}
				if (freedom =="range_restricted" ) {
					if (GetValue("restricted_range").size() ==0) {
						success=false;
						cout<<"In mol '" + name + "', freedom is set to 'range_restricted'. In this case we expect the setting for 'restricted_range'. This setting was not found. Problem terminated. " << endl;
					} else { //read range;
						int *HP=NULL;
						int M=Lat[0]->M;
						int npos=0;
						bool block;
						R_mask=(int*)malloc(M*sizeof(int));
						string s="restricted_range";
						int *r=(int*) malloc(6*sizeof(int));
						success=Lat[0]->ReadRange(r,HP,npos,block,GetValue("restricted_range"),0,name,s);
						Lat[0]->CreateMASK(R_mask,r,HP,npos,block);
						theta_range = theta;
						n_range = theta_range/chainlength;
						free(r);
					}

					//cout <<"restricted_range not yet implemented " << endl;
					//success=false;
				}

			} else {
				if (GetValue("theta").size() >0 || GetValue("n").size() > 0 || GetValue("phibulk").size() >0 || GetValue("freedom").size() > 0) cout <<"Warning. In mol " + name + " tagged segment(s) were detected. In this case no value for 'freedom' is needed, and also 'theta', 'n' and 'phibulk' values are ignored. " << endl;
			}
		}


		ring=false;
		if (GetValue("ring").size() > 0) {
			In[0]->Get_bool(GetValue("ring"),ring,"Input for ring is either 'true' or 'false'. Moreover, first and last segments of the backbone will be put on top of each other (chain length gets shorter by one). ");
			if (ring) {
				int length;
				switch (MolType) {
					case dendrimer:
					case asym_dendrimer:
						cout << "Do not know how to connect ends for dendrimer or asym_dendrimers; ring option is ignored " << endl;
						ring=false;
					break;
					case monomer:
						cout << "Can not make a ring from a molecule type monomer; ring option is ignored" << endl;
						ring=false;
					break;
					case linear:
					case branched:
					case comb:
						length=mon_nr.size();
						if (mon_nr[0]==mon_nr[length-1]) {
							if (n_mon[0] !=1 || n_mon[length-1] !=1) {
								cout <<"Fist and last block of main chain should consist of just one segment: e.g. (A)1(B)99(A)1 is a ring of 100 segments (one A only)."  << endl;
								ring =false;
							} else  {
								cout <<"ring can be implemented " << endl;
								cout <<"chain length reduced by one" << endl;
								chainlength--;
							}
						} else {
							ring =false;
							cout <<"first and last segment of the main chain should be of the same type, as these two are 'merged'. Chain length is reduced by one." << endl;
							cout <<"make sure to start the main chain and end the main chain by a block with length 1: for example '(B)1(B)99(B)1' will be a ring of 100 segments." << endl;
						}
					break;
					default:
						ring=false;
					break;

				}
			}
		}

	}
	Markov=1;
	if (GetValue("Markov").size()>0) Markov=In[0]->Get_int(GetValue("Markov"),1);
	if (Markov<1 || Markov>2) {
		cout <<" Integer value for 'Markov' is by default 1 and may be set to 2 for some mol_types and fjc-choices only. Markov value out of bounds. Proceed with caution. " << endl; success = false;
	}
	if (Markov==2) Lat[0]->Markov=2;
	k_stiff=Lat[0]->k_stiff; //pick up 'default' value from lattice.
	if (GetValue("k_stiff").size()>0) {
		k_stiff=In[0]->Get_Real(GetValue("k_stiff"),k_stiff);
		if (k_stiff<0 || k_stiff>10) {
			success =false;
			cout <<" Real value for 'k_stiff' out of bounds (0 < k_stiff < 10). " << endl;
			cout <<" For Markov == 2: u_bend (theta) = 0.5 k_stiff theta^2, where 'theta' is angle for bond direction deviating from the straight direction. " <<endl;
			cout <<" You may interpret 'k_stiff' as the molecular 'persistence length' " << endl;
			cout <<" k_stiff is a 'default value'. Use molecular specific values to overrule the default when appropriate (future implementation....) " << endl;
		}
		if (Lat[0]->fjc>1 && Lat[0]->gradients>1) {
			success=false;
			cout <<" Work in progress.... Currently, Markov == 2 is implemented in gradients>1 for FJC_choices = 3 " << endl;
		}
	}
	size=0;
	if (Markov ==2) {
		if (Lat[0]->gradients==1) size = Lat[0]->FJC;
		if (Lat[0]->gradients==2) { //assume fjc=1...
			if (Lat[0]->lattice_type==hexagonal && !Lat[0]->stencil_full) size = 7;
			if (Lat[0]->lattice_type==hexagonal && Lat[0]->stencil_full) size = 12; //2*FJC-1
			if (Lat[0]->lattice_type==simple_cubic && Lat[0]->stencil_full) size = 2*Lat[0]->FJC-1;
		}
		if (Lat[0]->gradients==3) {
			if (Lat[0]->lattice_type==hexagonal && !Lat[0]->stencil_full) size = 12;
			if (Lat[0]->lattice_type==simple_cubic && Lat[0]->stencil_full) size = 6;
		}
	    cout <<"size = " << size << endl;
	} else size = 1;
	if (size==0) {success=false; cout <<"Attention: size in molecule is not set; combination gradients (1,2,3),  Markov=2, stencil_full (true,false), lattice_type (hexagonal, simple_cubic) not implemented" << endl;}

	return success;
}

bool Molecule::PutVarInfo(string Var_type_,string Var_target_,Real Var_target_value_){
if (debug) cout <<"Molecule:: PutVarInfo in mol "+ name << endl;

	bool success=true;
	//Var_scan_value=-1;
	//Var_search_value=-1;
	//Var_target=-1;
	vector<string>sub;
	Var_type="";
	if (Var_type_=="scan") {
		var_al_nr=-1;
		Var_type="scan";
		In[0]->split(Var_target_,'-',sub);
		if (sub.size()==1) {
			if (Var_target_=="theta") {Var_scan_value=0; Var_start_value=theta; }
			if (Var_target_=="n") {Var_scan_value=1; Var_start_value = n;}
			if (Var_target_=="phibulk") {Var_scan_value=2; Var_start_value = phibulk;}
			if (Var_scan_value==-1) {
				cout <<"In var: scanning value for molecule is not recognized. Choose from {theta, n, phibulk}. You can also ask for 'aliasname-value'. " << endl;
				return false;
			}
		} else {
			Var_scan_value=3;
			int length=MolAlList.size();
			for (int k=0; k<length; k++) {
				if (Al[k]->name==sub[0] && "value" ==sub[1]) {
					if (Al[k]->value >0) {
						var_al_nr=k;
						Var_start_value=Al[k]->value;
					}
				}
			}
			if (var_al_nr==-1) {
				cout <<"In var: scanning value for molecule, the scanning value 'aliasname-value' was not recognised. The first element should contain a valid aliasname; the second element must contain the keyword 'value' " << endl;
			}
		}

	}
	if (Var_type_=="search") {
		Var_type="search";
		if (Var_target_=="theta") {Var_search_value=0; Var_start_search_value=theta;  }
		if (Var_target_=="n") {Var_search_value=1; Var_start_search_value=n; }
		if (Var_target_=="phibulk") {Var_search_value=2; Var_start_search_value=phibulk; }
		if (Var_target_=="equate_to_solvent") {
			Var_search_value=3; Var_start_search_value=theta;
			if (freedom=="solvent") {
				cout <<"in var: searching for theta to equate to solvent can not be done for a molecule with 'freedom' solvent " << endl;
				success=false;
			}
			if (freedom!="restricted") {
				success=false;
				cout <<"In var: searching for theta to equate to solvent can only be done for a molecule with 'freedom' restricted" << endl;
			}
		}
		if (Var_target_=="balance_membrane") {
			Var_search_value=4; Var_start_search_value=theta;
			if (freedom=="solvent") {
				cout <<"in var: searching for theta to balance membrane can not be done for a molecule with 'freedom' solvent " << endl;
				success=false;
			}
			if (freedom!="restricted") {
				success=false;
				cout <<"In var: searching for theta to balance membrane can only be done for a molecule with 'freedom' restricted" << endl;
			}
		}

		if (Var_search_value == -1) {
			cout <<"In var: searching value for molecule was not recognized. Choose from {theta, n, phibulk, equate_to_solvent}. " << endl;
			return false;
		}
	}
	if (Var_type_=="target") {
		Var_type="target";
		Var_target_value=Var_target_value_;
		if (Var_target_=="theta") {
			Var_target=0;
			if (Var_target_value <0 || Var_target_value>Lat[0]->volume){
				cout <<"In var: target value 'theta' out of range" << endl;
				return false;
			}
		}
		if (Var_target_=="n") {
			Var_target=1;
			if (Var_target_value <0 || Var_target_value*chainlength>Lat[0]->volume){
				cout <<"In var: target value 'n' out of range" << endl;
				return false;
			}

		}
		if (Var_target_=="phibulk") {
			Var_target=2;
			if (Var_target_value <0 || Var_target_value>1) {
				cout <<"in var: target value 'phibulk' out of range" << endl;
				return false;
			}
		}
		if (Var_target_=="mu") {Var_target=3; if (Var_target_value==12345.0) Var_target_value=0;}
		if (Var_target==-1) {
			cout <<"In var: molecule target should be selected from {theta, n, phibulk, mu}. " << endl;
			success = false;
		}
	}
	return success;

}

int Molecule::PutVarScan(Real step, Real end_value, int steps, string scale_) {
if (debug) cout <<"Molecule:: PutVarScan" << endl;
	num_of_steps = -1;
	scale=scale_;
	Var_end_value=end_value;
	if (scale=="exponential") {
		Var_steps=steps; Var_step=0;
		if (steps==0) {
			cout <<"In var scan: the value of 'steps' is zero, this is not allowed" << endl; return -1;
		}
		if (Var_end_value*Var_start_value<0) {
			cout <<"In var scan: the product end_value*start_value <0. This is not allowed. " << endl; return -1;
		}
		if (Var_end_value > Var_start_value)
			num_of_steps=steps*log10(Var_end_value/Var_start_value);
		else
			num_of_steps=steps*log10(Var_start_value/Var_end_value);
	} else {
		Var_steps=0; Var_step=step;
		if (step==0) {
			cout <<"In var san: of molecule variable, the value of step can not be zero" << endl; return -1;
		}
		num_of_steps=(Var_end_value-Var_start_value)/step;

		if (num_of_steps<0) {
			cout <<"In var scan : (end_value-start_value)/step is negative. This is not allowed. Try changing the sign of 'step'. " << endl;
			return -1;
		}

	}
	return num_of_steps;
}

bool Molecule::ResetInitValue() {
if (debug) cout <<"Molecule:: ResetInitValue" << endl;
	bool success=true;
	cout <<"reset: ";
	switch (Var_scan_value) {
		case 0:
			theta=Var_start_value;
			n=theta/chainlength;
			cout <<"mol : " + name + " : theta : " << theta << endl;
			break;
		case 1:
			n=Var_start_value;
			theta=n*chainlength;
			cout <<"mol : " + name + " : n : " << n << endl;
			break;
		case 2:
			phibulk=Var_start_value;
			cout <<"mol : " + name + " : phibulk " << phibulk << endl;
			break;
		case 3:
			Gnr.clear();
			first_s.clear();
			last_s.clear();
			first_b.clear();
			last_b.clear();
			mon_nr.clear();
			n_mon.clear();
			molmon_nr.clear();
			memory.clear();
			last_stored.clear();
			//MolAlList.clear();
			MolMonList.clear();
			Al[var_al_nr]->value=Var_start_value;
			Decomposition(GetValue("composition"));
			cout <<"alias : " + Al[var_al_nr]->name + " : value : " << Al[var_al_nr]->value << endl;
			break;
		default:
			break;
	}

	switch (Var_search_value) {
		case 0:
			theta=Var_start_search_value;
			n=theta/chainlength;
			cout <<"mol : " + name + " : theta : " << theta << endl;
			break;
		case 1:
			n=Var_start_search_value;
			theta=n*chainlength;
			cout <<"mol : " + name + " : n : " << n << endl;
			break;
		case 2:
			phibulk=Var_start_search_value;
			cout <<"mol : " + name + " : phibulk " << phibulk << endl;
			break;
		case 3: theta=Var_start_search_value;
			cout <<"mol : " + name + " : theta " << theta << endl;
		case 4: theta=Var_start_search_value;
			cout <<"mol : " + name + " : theta " << theta << endl;
		default:
			break;
	}

	return success;
}

bool Molecule::UpdateVarInfo(int step_nr) {
if (debug) cout <<"Molecule:: UpdateVarInfo" << endl;
	bool success=true;
	switch(Var_scan_value) {
		case 0:
			if (scale=="exponential") {
				theta=pow(10,(1-1.0*step_nr/num_of_steps)*log10(Var_start_value)+(1.0*step_nr/num_of_steps)*log10(Var_end_value));
			} else {
				theta=Var_start_value+step_nr*Var_step;
			}
			n=theta/chainlength;
			cout <<"mol : " + name + " : theta : " << theta << endl;
			break;
		case 1:
			if (scale=="exponential") {
				n=pow(10,(1-1.0*step_nr/num_of_steps)*log10(Var_start_value)+(1.0*step_nr/num_of_steps)*log10(Var_end_value));
			} else {
				n=Var_start_value+step_nr*Var_step;
			}
			cout <<"mol : " + name + " : n : " << n << endl;
			theta=n*chainlength;
			break;
		case 2:
			if (scale=="exponential") {
				phibulk=pow(10,(1-1.0*step_nr/num_of_steps)*log10(Var_start_value)+(1.0*step_nr/num_of_steps)*log10(Var_end_value));
			} else {
				phibulk=Var_start_value+step_nr*Var_step;
			}
			cout <<"mol : " + name + " : phibulk " << phibulk << endl;
			break;
		case 3:
			if (scale=="exponential") {
				Al[var_al_nr]->value =(int)pow(10,(1-1.0*step_nr/num_of_steps)*log10(Var_start_value) + (1.0*step_nr/num_of_steps)*log10(Var_end_value));
			} else {
				Gnr.clear();
				first_s.clear();
				last_s.clear();
				first_b.clear();
				last_b.clear();
				mon_nr.clear();
				n_mon.clear();
				molmon_nr.clear();
				memory.clear();
				last_stored.clear();
				//MolAlList.clear();
				MolMonList.clear();
				Al[var_al_nr]->value =(int) (Var_start_value+step_nr*Var_step);
				Decomposition(GetValue("composition"));
			}
			cout <<"alias : " + Al[var_al_nr]->name + " : value : " << Al[var_al_nr]->value << endl;
			break;
		default:
			cout <<"program error in Molecule::UpdateInfo " << endl;
			break;
	}
	return success;
}

Real Molecule::GetError() {
if (debug) cout <<"Molecule:: GetError" << endl;
	Real Error=0;
	switch (Var_target) {
		case 0:
			if (Var_target_value!=0) Error=theta/Var_target_value-1.0; else Error = theta;
			break;
		case 1:
			if (Var_target_value!=0) Error=n/Var_target_value-1.0; else Error = n;
			break;
		case 2:
			if (Var_target_value!=0) Error=phibulk/Var_target_value-1.0; else Error = phibulk;
			Error *=-1;
			break;
		case 3:
			if (Var_target_value!=0) Error=Mu/Var_target_value-1.0; else Error = Mu;
			break;
		default:
			cout <<"program error in Molecule::GetError" << endl;
			break;
	}

	if (Var_search_value==3||Var_search_value==4) {
		Error=theta;
		cout<<"Get Error in mol " + name + "unexpectedly called. Shouldn't sys[0]->GetError be called? Returning the value for theta " << endl;
	}
	return Error;
}

Real Molecule::GetValue(){
if (debug) cout <<"Molecule:: GetValue" << endl;
	Real X=0;
	switch (Var_search_value) {
		case 0:
			X=theta;
			break;
		case 1:
			X=n;
			break;
		case 2:
			X=phibulk;
			break;
		case 3:
			X=theta;
			break;
		case 4:
			X=theta;
			break;
		default:
			cout <<"program error in Molecule::GetValue" << endl;
	}
	return X;
}
void Molecule::PutValue(Real X){
if (debug) cout <<"Molecule:: PutValue" << endl;
	switch (Var_search_value) {
		case 0:
			theta=X; n=theta/chainlength;
			break;
		case 1:
			n=X; theta=n*chainlength;
			break;
		case 2:
			phibulk=X;
			break;
		case 3:
			theta=X;
			break;
		case 4:
			theta=X;
			break;
		default:
			cout <<"program error in Molecule::GetValue" << endl;
	}
}

int Molecule::GetAlNr(string s){
if (debug) cout <<"GetAlNr for Mol " + name << endl;
	int n_als=MolAlList.size();
	int found=-1;
	int i=0;
	while(i<n_als) {
		if (Al[i]->name ==s) found=i;
		i++;
	}
	return found;
}

int Molecule::GetMonNr(string s){
if (debug) cout <<"GetMonNr for Mon " + name << endl;
	int n_segments=In[0]->MonList.size();
	int found=-1;
	int i=0;
	while(i<n_segments) {
		if (Seg[i]->name ==s) found=i;
		i++;
	}
	return found;
}

bool Molecule::ExpandAlias(vector<string> sub, string &s) {
if (debug) cout <<"Molecule:: ExpandAlias" << endl;
	bool success=true;
	vector<int> open;
	vector<int> close;
	int length_al=sub.size();
	int i=0;
	while (i<length_al-1) {
		string sA;
		sA=sub[i+1];
		if (!In[0]->InSet(In[0]->AliasList,sA)) {
			cout <<"In composition of mol '" + name + "' Alias '" + sA + "' was not found"<<endl; success=false; return success;
		} else {
			int Alnr =GetAlNr(sA);
			if (Alnr<0) {
				Al.push_back(new Alias(In,Lat,sA));
				Alnr=Al.size();
				if (!Al[Alnr-1]->CheckInput(start)) {cout <<"Alias '" + sA + "' in composition not recognised " << endl; return false;}
				MolAlList.push_back(Alnr);
			}
			Alnr =GetAlNr(sA);
			int iv = Al[Alnr]->value;
			string al_comp=Al[Alnr]->composition;
			if (iv < 0) {
				string si;
				stringstream sstm;
				sstm << Alnr;
				si = sstm.str();
				string ssub="";
				sub[i+1]=":"+si+":"+al_comp+":"+si+":";

			} else {
				string sss;
				stringstream sstm;
				sstm << iv;
				sss = sstm.str();
				sub[i+1]=sss;
			}
		}
		i+=2;
	}

	string ss;
	for (int i=0; i<length_al; i++) {
		ss=ss.append(sub[i]);
	}
	s=ss;
	//cout <<"expanded to " << s << endl;
	return success;
}

bool Molecule::ExpandBrackets(string &s) {
if (debug) cout <<"Molecule:: ExpandBrackets" << endl;
	bool success=true;
	if (s[0] != '(') {cout <<"illegal composition. Expects composition to start with a '(' in: " << s << endl; return false;}
	vector<int> open;
	vector<int> close;
	bool done=false; //now interpreted the (expanded) composition
	while (!done) { done = true;
		open.clear(); close.clear();
		if (!In[0]->EvenBrackets(s,open,close)) {
			cout << "s : " << s << endl;
			cout << "In composition of mol '" + name + "' the backets are not balanced."<<endl; success=false; return success;
		 }
		int length=open.size();
		int pos_open;
		int pos_close;
		int pos_low=0;
		int i_open=0; {pos_open=open[0]; pos_low=open[0];}
		int i_close=0; pos_close=close[0];
		if (pos_open > pos_close) {cout << "Brackets open in composition not correct" << endl; return false;}
		while (i_open < length-1 && done) {
			i_open++;
			pos_open=open[i_open];
			if (pos_open < pos_close && done) {
				i_close++; if (i_close<length) pos_close=close[i_close];
			} else {
				if (pos_low==open[i_open-1] ) {
					pos_low=pos_open;
					i_close++; if (i_close<length) pos_close=close[i_close];
					if (pos_open > pos_close) {cout << "Brackets open in composition not correct" << endl; return false;}
				} else {
					done=false;

					int x=In[0]->Get_int(s.substr(pos_close+1),-1);
					if (x<1) {
							cout <<"Number of 'repeats' smaller or equal to zero (or '# repeats' is missing) in composition at pos : " << pos_close+1 << " for: " << s << endl; return false;
					}
					string sA,sB,sC;
					if (s.substr(pos_open-1,1)=="]") {pos_open --;  }
					sA=s.substr(0,pos_low);
					sB=s.substr(pos_low+1,pos_close-pos_low-1);
					sC=s.substr(pos_open,s.size()-pos_open+1);
					s=sA;for (int k=0; k<x; k++) s.append(sB); s.append(sC);
				}
			}
		}
		if (pos_low < open[length-1]&& done) {
			done=false;
			pos_close=close[length-1];
			int x=In[0]->Get_int(s.substr(pos_close+1),0);
			string sA,sB,sC;
			sA=s.substr(0,pos_low);
			sB=s.substr(pos_low+1,pos_close-pos_low-1);
			sC="";
			s=sA;for (int k=0; k<x; k++) s.append(sB); s.append(sC);
		}

	}
	if (s[s.size()-1]==']') {cout <<"illegal composition. Composition can not end with a ']' in: " << s << endl; return false;}
//cout <<"expand brackets to " << s << endl;
	return success;
}

bool Molecule::Interpret(string s,int generation){
if (debug) cout <<"Molecule:: Interpret" << endl;
	if (s=="[") return true;
	bool success=true;
	vector<string>sub;
	vector<int>open;
	vector<int>close;
	In[0]->split(s,':',sub);
	int length_sub =sub.size();
	int AlListLength=MolAlList.size();
	int i=0;

	while (i<length_sub) {
		open.clear(); close.clear();
		In[0]->EvenBrackets(sub[i],open,close);
		if (open.size()==0) {
			int a=In[0]->Get_int(sub[i],0);
			if (Al[a]->active) Al[a]->active=false; else Al[a]->active=true;
		} else {
			int k=0;
			int length=open.size();

			while (k<length) {
				string segname=sub[i].substr(open[k]+1,close[k]-open[k]-1);
				int mnr=GetMonNr(segname);
				if (mnr <0)  {cerr <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised"  << endl; success = false;
				// The entire chain of successes seems to be broken, as they are not transferred to or handled correctly by functions down the stack.
				// This function is called multiple times with different results, causing the success handling to break by only transferring the last success result.
				// Throwing to prevent segfaults and other undefined behavior. Caught by CheckInput.
				throw "Composition Error";
				} else {

					int length=Gnr.size();
					if (length>0) {//fragments at branchpoint need to be just 1 segment long.
						if (Gnr[length-1]<generation) {
							if (n_mon[length-1]>1) {
								n_mon[length-1]--;
								n_mon.push_back(1);
								mon_nr.push_back(mon_nr[length-1]);
								Gnr.push_back(Gnr[length-1]);
							 	last_b[Gnr[length-1]]++;
							}
						}
					}
					mon_nr.push_back(mnr);
					Gnr.push_back(generation);
					if (first_s[generation] < 0) first_s[generation]=chainlength;
					if (first_b[generation] < 0) first_b[generation]=mon_nr.size()-1;
					last_b[generation]=mon_nr.size()-1;
					for (int i=0; i<AlListLength; i++) {if (Al[i]->active) Al[i]->frag.push_back(1); else Al[i]->frag.push_back(0);}
				}
				int nn = In[0]->Get_int(sub[i].substr(close[k]+1,s.size()-close[k]-1),0);
				if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity " << endl; success=false; return success;
				//throw "Composition error";
				} else {
					n_mon.push_back(nn);
				}
				chainlength +=nn; last_s[generation]=chainlength;
				k++;
			}
		}
		i++;
	}
	return success;
}

bool Molecule::GenerateTree(string s,int generation,int &pos, vector<int> open,vector<int> close) {
if (debug) cout <<"Molecule:: GenerateTree" << endl;
	bool success=true;
	string ss;
	int i=0;
	int newgeneration=0;
	int new_generation=0;
	int length=open.size();
	int pos_open=0;
	int pos_close=s.length();
	bool openfound,closedfound;
	while  (pos_open<pos_close && success) {
		pos_open =s.length()+1;
		pos_close=s.length();
		openfound=closedfound=false;
		i=0;
		while (i<length && !(openfound && closedfound) ){

			if (close[i]>pos && !closedfound) {closedfound=true; pos_close=close[i]+1; new_generation=i+1;}
			if (open[i]>=pos && !openfound) {openfound=true; pos_open=open[i]+1; newgeneration=i+1;}
			i++;
		}

		if (pos_close<pos_open) {
			ss=s.substr(pos,pos_close-pos);
			if (ss.substr(0,1)=="[") {
				pos=pos+1;
				first_s.push_back(-1);
				last_s.push_back(-1);
				first_b.push_back(-1);
				last_b.push_back(-1);
				success=GenerateTree(s,new_generation,pos,open,close);
				if (!success) {cout <<"error in generate tree ." << endl; return success; }
				pos_close=pos_open+1;
			} else {
				pos=pos_close;
				success =Interpret(ss,generation);
				if (!success)  {cout <<"error in interpret ." << endl; return success; }
			}
		} else {
			ss=s.substr(pos,pos_open-pos);
			pos=pos_open;
			success =Interpret(ss,generation);
			if (!success) {cout <<"error in Interpret" << endl;  return success;}
			first_s.push_back(-1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			success=GenerateTree(s,newgeneration,pos,open,close);
			if (!success) {cout <<"error in generate tree .." << endl; return success;}
		}
	}
	return success;
}

//bool Molecule::GenerateDend(string s, int generation) {
//if (debug) cout <<"Molecule:: GenerateDend" << endl;
//return true;
//}

bool Molecule::Decomposition(string s){
if (debug) cout <<"Decomposition for Mol " + name << endl;
	bool success = true;
	bool aliases = true;
	MolType=linear;//default;
	chainlength=0;
	int loopnr=0;
	vector<int> open;
	vector<int> close;
	vector<string>sub;

	sub.clear();
	In[0]->split(s,'@',sub);

	if (s!=sub[0]) {
		bool keyfound=false;
		vector<string>SUB;
		In[0]->split(sub[1],'(',SUB);
		if (SUB[0]=="water"){
			MolType=water; keyfound=true;
			s=s.substr(7,s.length()-8);
			int mnr=GetMonNr(s);
			if (mnr<0) {
				cout <<"Language for MolType 'water' is as follows: @water(mon_name) wherein mon_name is a valid monomer name with freedom free" << endl;
				success=false; return success;
			} else {
				mon_nr.push_back(mnr);
				n_mon.push_back(1);
			}
		}

		if (SUB[0]=="dend") {
			MolType=dendrimer; keyfound=true;
			s=s.substr(6,s.length()-7);
			if (save_memory) {success=false; cout <<"In dendrimer no save_memory implemented yet. contact frans.leermakers@wur.nl....." << endl; return success;}
			//if (aliases) {success=false; cout <<"In dendrimer no aliases are not implemented yet. Can be done relatively easily....." << endl; return success;}
		}
		if (SUB[0]=="comb") {
			if (save_memory) {success=false; cout <<"In comb no save_memory implemented yet. contact frans.leermakers@wur.nl ....." << endl; return success;}
			MolType=comb; keyfound=true;
			s=s.substr(6,s.length()-7);
		}
		if (!keyfound) { success=false; cout << "Keyword specifying Moltype not recognised: select from @dend, @comb. @water Problem terminated "<< endl ;
			return success;
		 }
	}

	while (aliases && success) {
		loopnr++;
		sub.clear();

		In[0]->split(s,'#',sub);
		if (sub.size()%2!=1) {
			if (s.back()!='#') {
				cout << " Alias in composition should be bracketed on both sides by '#'. For example, (A)10#alias_name#(B)3, when e.g., 'alias : alias_name : value : (X)5' is defined, so that you oubtain (A)10(X)5(B)3 " << endl; success=false;
			}
		}
		aliases=(s!=sub[0]);
		if (aliases) {if (!ExpandAlias(sub,s)) {cout << "expand alias failed. " << endl; success=false; return success; }}
		if (loopnr == 20) {
			cout << "Nesting nr 20 reached in aliases for mol " + name + " -composition. It is decided that this is too deep to continue; Possible, you have defined an alias-A inside alias-B which itself refers to alias-A. This is not allowed. Problem terminated. " << endl;
			success=false; return success;
		}
	}
	sub.clear();
	In[0]->split(s,',',sub);
	int length=sub.size();
	int length_open;
	int i,j,k,a,f,dd;
	string ss;
	switch(MolType) {
		case water:

			break;
		case dendrimer:
			for (i=0; i<length; i++) {
				open.clear(); close.clear();
				In[0]->EvenBrackets(sub[i],open,close);
				length_open=open.size();
				if (length_open>0) {
					if (!ExpandBrackets(sub[i])) {
						cout <<"brackets '(' ')' not well positioned in "+sub[i] << endl;
          				success=false; return success;
					}
				}
			}
			for (i=0; i<length; i++) {
				ss.append(sub[i]);
				if (i<length-1) ss.append(",");
			}
			s=ss;
		break;
		case comb:
			//success=false;
		break;
		default:
			if (!ExpandBrackets(s)) {success=false; return success;}
		break;
	}
	if (!In[0]->EvenSquareBrackets(s,open,close)) {
		cout << "Error in composition of mol '" + name + "'; the square brackets are not balanced in " << s << endl;
		success=false; return success;
	}
	if (open.size()>0) {
		if (MolType==linear) { //overwrite default.
			MolType=branched;
		} else {
			success=false;
			cout <<" In 'composition' you can not combine special keywords such as 'dend, or comb' with branched compositions. This implies that square brackets are not allowed. " <<endl ;
			return success;
		}
	}
	int generation=0;
	int pos=0;
	MolMonList.clear();
	int AlListLength=MolAlList.size();
	for (int i=0;i<AlListLength; i++) Al[i]->active=false;
	vector<string>sub_gen;
	vector<string>sub_dd;
	int length_g,length_dd,mnr,nn,arm=-1,degeneracy=1,arms=0;
	string segname;
	int N=0;
	int lopen;
	int chainlength_backbone,chainlength_arm;
	switch(MolType) {
		case water:
			break;
		case linear:
			first_s.clear();
			last_s.clear();
			first_b.clear();
			last_b.clear();
			first_s.push_back(-1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			success = GenerateTree(s,generation,pos,open,close);
			if (!success) {cout << " GenerateTree failed " << endl; return success; }
			break;
		case branched:
			first_s.clear();
			last_s.clear();
			first_b.clear();
			last_b.clear();
			first_s.push_back(-1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			lopen=open.size();
			for (int i=1; i<lopen; i++) {
				if ((open[i]-open[i-1])==1 || (close[i]-close[i-1])==1) {cout <<"In molecule " + name + " in 'composition', two similar square brackets in a row '[[' or ']]' is not allowed" << endl; success=false; return success;  }
			}
			success = GenerateTree(s,generation,pos,open,close);
			if (!success) {cout <<"GenerateTree failed" << endl; return success;}
			break;
		case dendrimer:
			first_a.clear();
			last_a.clear();
			first_b.clear();
			last_b.clear();
			first_s.clear();
			last_s.clear();
			n_arm.clear();
			mon_nr.clear();
			n_mon.clear();
			In[0]->split(s,';',sub_gen);
			n_generations=sub_gen.size();
			sym_dend=true; //default.
			if (save_memory) {
				success = false;
				cout <<"In dendrimer the use of 'save_memory' is not allowed (yet). " << endl; return success;
			}


			//cout <<"n_generations " << n_generations << endl;
			chainlength=0; N=-1;
			for (i=0; i<n_generations; i++) {
				sub.clear();	arms=0;
				In[0]->split(sub_gen[i],',',sub);
				int sublength=sub.size();
				if ((sublength-1)%2 >0) {success=false; cout << "In composition of dend for generation " << i << " that is, in " + sub_gen[i] + " the number of arguments is not a multiple of (1+ multiple of 2); use @dend(?) for help. " << endl; }
				if (sublength>3) MolType=asym_dendrimer;
				if (sublength<2) { success=false;

					cout<<" ------Dendrimer language:------ " << endl;
					cout<<" example 1: consider segment name 'A' for the branch points and spacers (B)2 with functionality 3 and 4 generations: " << endl;
					cout<<" @dend(A,(B)2,3; A,(B)2,3; A,(B)2,3; A,(B)2,3)  " << endl;
					cout<<" hence generations are seperated by ';', branch points are 'monomer_names' without '(', ')'. Each generation can have unique numbers and structures. " << endl;
					cout<<" expample 2: asymmetric dendrimers, with asymmetry in branches of first generation:" << endl;
					cout <<" @dend(A,(B)2,3,(B)4,1; A,(B)2,3,(B)2,2; A,(B)2,3; A,(B)2,3)  " << endl;
					cout <<" example 2: star with 10 arms " << endl;
					cout <<" @dend(A,(B)100,10) " << endl;
					cout <<" Note that in the 'standard' dendrimer the branches are symmetric -all are equal- " << endl;
					cout <<" example 4: 10 arm star end-grafted by one arm to a surface by way of segment 'X' " << endl;
					cout <<" @dend(A,(B)100,9,(B)99(X)1,1)" << endl;
					cout<<" ------Dendrimer language:------ " << endl;
					return success;
				}

				length_g = sub.size();
				sub_dd.clear();
				In[0]->split(sub[0],':',sub_dd);
				length_dd =sub_dd.size();
				if (sub_dd.size()>1) {
					cout <<"In dendrimer alias not yet implemented" << endl; success=false; return success;
				}


				if (length_dd==4) {
					mnr=GetMonNr(sub_dd[2]);
					a=In[0]->Get_int(sub_dd[2],0);
					if (Al[a]->active) Al[a]->active=false; else Al[a]->active=true;

				} else mnr=GetMonNr(sub[0]);
				if (mnr <0)  {success=false; cout <<"In composition of mol '" + name + "', segment name '" + sub_dd[0] + "' is not recognised"  << endl; }
				for (a=0; a<AlListLength; a++) {if (Al[a]->active) Al[a]->frag.push_back(1); else Al[a]->frag.push_back(0);}

				n_mon.push_back(1);
				mon_nr.push_back(mnr);
				d_mon.push_back(degeneracy);
				first_a.push_back(-1);
				last_a.push_back(-1);
				k=1; chainlength+=degeneracy; N++;

				while (k<length_g-1) {
					arm++;
					first_s.push_back(N+1);
          last_s.push_back(-1);
          first_b.push_back(-1);
          last_b.push_back(-1);
					if (first_a[first_a.size()-1]==-1) first_a[first_a.size()-1]=arm;
					last_a[last_a.size()-1]=arm;

					f=In[0]->Get_int(sub[k+1],0); //should not contain double dots....
					if (f<1) {
						success=false; cout <<"In dendrimer-composition, in generation "<<i << " an integer number is expected at argument " << k+1 << " problem terminated" << endl;
					}
					n_arm.push_back(f); arms+=f;

					sub_dd.clear();
					In[0]->split(sub[k],':',sub_dd);
					dd=0; length_dd=sub_dd.size();
					while (dd<length_dd) {
						open.clear(); close.clear();
            In[0]->EvenBrackets(sub_dd[dd],open,close);
						if (open.size()==0) {
							a=In[0]->Get_int(sub_dd[dd],-1);
							if (a==-1) {
								cout <<"No integer found. Possibly you have a segment name in composition that is not surrounded by brackets " << endl; success=false;
							} else {if (Al[a]->active) Al[a]->active=false; else Al[a]->active=true;}
						} else {
							j=0; length=open.size();
							while (j<length) {
								segname=sub_dd[dd].substr(open[j]+1,close[j]-open[j]-1);
								mnr=GetMonNr(segname);
								if (mnr<0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised; this occurs at generation " <<i << " arm " << k << "."  << endl; success=false;}
								mon_nr.push_back(mnr);
								nn=In[0]->Get_int(sub_dd[dd].substr(close[j]+1,s.size()-close[j]-1),0);
								if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity; this occurs at generation " <<i << " arm " <<k <<"."<< endl; success=false;}
								n_mon.push_back(nn); N+=nn;
								d_mon.push_back(degeneracy*f);
								chainlength +=degeneracy*nn*f;
                if (first_b[first_b.size()-1] < 0) first_b[first_b.size()-1]=mon_nr.size()-1;
                last_b[first_b.size()-1]=mon_nr.size()-1;
								last_s[last_s.size()-1]=N;
								for (a=0; a<AlListLength; a++) {if (Al[a]->active) Al[a]->frag.push_back(1); else Al[a]->frag.push_back(0);}
								j++;
							}
						}
						dd++;
					}
					k=k+2;
				}
				degeneracy*=arms;
			}
   //anticipated that asymmetric dendrimers will be of interest in the future
//for (int i=0; i<n_generations; i++) cout<<"generation " << i << " first_a " << first_a[i] << " last_a " << last_a[i] << endl;
//length=n_arm.size();
//for (int i=0; i<length; i++) cout <<"arm " << i << " first_b " << first_b[i] << " last_b " << last_b[i] << endl;
//for (int i=0; i<length; i++) cout <<"arm " << i << " first_s " << first_s[i] << " last_s " << last_s[i] << endl;
//for (int i=0; i<length; i++) cout <<"arm " << i << " n_arm " << n_arm[i] << endl;
//length=n_mon.size();
//for (int i=0; i<length; i++) cout <<"block " << i << " n_mon " << n_mon[i] << " mon_nr " << mon_nr[i] << " d_mon " << d_mon[i] << endl;
//cout <<"Chain length =" << chainlength << endl;
//cout <<" N = " << N << endl;
//cout <<"Composition = " <<  s << endl;

			break;
		case comb:
			first_a.clear();
			last_a.clear();
			first_b.clear();
			last_b.clear();
			first_s.clear();
			last_s.clear();
			n_arm.clear();
			mon_nr.clear();
			n_mon.clear();
			In[0]->split(s,';',sub_gen);
			n_generations=sub_gen.size();
			int j;
			int mnr;
			int nn;

			sub_dd.clear();
			In[0]->split(s,':',sub_dd);
			length_dd=sub_dd.size();
			if (length_dd>1) {
				success = false;
				cout <<"In comb the use of 'aliases' is not allowed (yet). " << endl; return success;
			}
			if (save_memory) {
				success = false;
				cout <<"In comb the use of 'save_memory' is not allowed (yet). " << endl; return success;
			}

			if (n_generations != 3) {
				success=false;
				cout<<" ------Comb language:------ " << endl;
				cout<<" generic example: @comb((A)5; A,(B)3,(A)6,4; (A)7 )" << endl;
				cout<<" Backbnone has a first part of 5 A-segments, then spacers with 6 A units and a trailing of 7 A-segments. " << endl;
				cout<<" Teeth are composed of 3 B-segments " << endl;
				cout<<" Number of repeats (inbetween the two ';') is eqaul to 4. " << endl;
				cout<<" Number of segments in molelecule is 5+(1+3+6)x4+7 = 52 segments" << endl;
				cout<<" Backbone length is 5+(1+6)x4+7 = 40 segments " << endl;
				cout<<" ------Comb language:------ " << endl;
				return success;
			}
			chainlength=0; N=-1;
			i=0;
			//success=Interpret(sub_gen[0],i);

			first_s.push_back(N+1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			open.clear(); close.clear();
			In[0]->EvenBrackets(sub_gen[0],open,close);
			j=0; length=open.size();
			while (j<length) {
				segname=sub_gen[0].substr(open[j]+1,close[j]-open[j]-1);
				mnr=GetMonNr(segname);
				if (mnr<0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised; this occurs at generation " << endl; success=false;}
				mon_nr.push_back(mnr); d_mon.push_back(1);
				nn=In[0]->Get_int(sub_gen[0].substr(close[j]+1,s.size()-close[j]-1),0);
				if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity; this occurs at generation " << endl; success=false;}
				n_mon.push_back(nn); N+=nn;
				chainlength +=nn;
				if (first_b[first_b.size()-1] < 0) first_b[first_b.size()-1]=mon_nr.size()-1;
				last_b[first_b.size()-1]=mon_nr.size()-1;
				last_s[last_s.size()-1]=N;
				j++;
			}

			chainlength_backbone=chainlength;
			i=1;
			sub.clear();
			In[0]->split(sub_gen[1],',',sub);
			//cout <<"sub_gen[1] " << sub_gen[1] << "size of sub : " << sub.size() << endl;
			if (sub.size() != 4) {
				cout << "Central part in comb definition should contain four arguments and three ',' to separare them. Use comb(? ) for details."  << endl;
				success=false; return success;
			}

			n_arm.push_back(In[0]->Get_int(sub[3],0));
			first_s.push_back(N+1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			open.clear(); close.clear();
			In[0]->EvenBrackets(sub[1],open,close);
			j=0; length=open.size();
			while (j<length) {
				segname=sub[1].substr(open[j]+1,close[j]-open[j]-1);
				mnr=GetMonNr(segname);
				if (mnr<0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised; this occurs at generation " << endl; success=false; return success; }
				mon_nr.push_back(mnr); d_mon.push_back(n_arm[0]);
				nn=In[0]->Get_int(sub[1].substr(close[j]+1,s.size()-close[j]-1),0);
				if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity; this occurs at generation " << endl; success=false; return success; }
				n_mon.push_back(nn); N+=nn;
				chainlength +=nn;
				if (first_b[first_b.size()-1] < 0) first_b[first_b.size()-1]=mon_nr.size()-1;
				last_b[first_b.size()-1]=mon_nr.size()-1;
				last_s[last_s.size()-1]=N;
				j++;
			}

			chainlength_arm=chainlength-chainlength_backbone;
			chainlength=chainlength_backbone;

			if (n_arm[0] <1) {
				success=false; cout <<" Error in composition of mol "+ name + " number of arms is less than unity. Use comb(? ) for details. " << endl; return success;
			}
			n_generations=3;
			for (int a=1; a<=n_arm[0]; a++) {
				mnr=GetMonNr(sub[0]);
				if (mnr <0)  {success=false; cout <<"In composition of mol '" + name + "', segment name '" + sub_dd[0] + "' is not recognised. For the branching point the ( ) are not needed."  << endl; return success;}


				n_mon.push_back(1); d_mon.push_back(1);
				mon_nr.push_back(mnr); N++;
				chainlength +=chainlength_arm+1;
				n_generations++;
				//success=Interpret(sub[2],n_generations-1);

				first_s.push_back(N+1);
				last_s.push_back(-1);
				first_b.push_back(-1);
				last_b.push_back(-1);
				open.clear(); close.clear();
				In[0]->EvenBrackets(sub[2],open,close);
				j=0; length=open.size();
				while (j<length) {
					segname=sub[2].substr(open[j]+1,close[j]-open[j]-1);
					mnr=GetMonNr(segname);
					if (mnr<0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised; Use ring(?) for details." << endl; success=false; return success; }
					mon_nr.push_back(mnr); d_mon.push_back(1);
					nn=In[0]->Get_int(sub[2].substr(close[j]+1,s.size()-close[j]-1),0);
					if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity; Use ring(? ) for details. "<< endl; success=false; return success; }
					n_mon.push_back(nn); N+=nn;
					chainlength +=nn;
					if (first_b[first_b.size()-1] < 0) first_b[first_b.size()-1]=mon_nr.size()-1;
					last_b[first_b.size()-1]=mon_nr.size()-1;
					last_s[last_s.size()-1]=N;
					j++;
				}

			}
			n_generations++;

			//success=Interpret(sub_gen[2],n_generations-1);
			first_s.push_back(N+1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			open.clear(); close.clear();
			In[0]->EvenBrackets(sub_gen[2],open,close);
			j=0; length=open.size();
			while (j<length) {
				segname=sub_gen[2].substr(open[j]+1,close[j]-open[j]-1);
				mnr=GetMonNr(segname);
				if (mnr<0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised; " << endl; success=false; return success; }
				mon_nr.push_back(mnr); d_mon.push_back(1);
				nn=In[0]->Get_int(sub_gen[2].substr(close[j]+1,s.size()-close[j]-1),0);
				if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity;  " << endl; success=false; return success; }
				n_mon.push_back(nn); N+=nn;
				chainlength +=nn;
				if (first_b[first_b.size()-1] < 0) first_b[first_b.size()-1]=mon_nr.size()-1;
				last_b[first_b.size()-1]=mon_nr.size()-1;
				last_s[last_s.size()-1]=N;
				j++;
			}

			//length=first_b.size();
			//for (int i=0; i<length; i++) cout <<"first_b " <<  first_b[i] << " last_b " << last_b[i] << endl;
			//for (int i=0; i<length; i++) cout <<"first_s " <<  first_s[i] << " last_s " << last_s[i] << endl;
			break;
					default:
			break;
	}
	if (MolType==branched) { //invert numbers;
		int g_length=first_s.size();
		int length = n_mon.size();
		int xxx;
		int ChainLength=last_s[0];

		for (int i=0; i<g_length; i++) {
			first_s[i] = ChainLength-first_s[i]-1;
			last_s[i] = ChainLength-last_s[i]-1;
			xxx=first_s[i]; first_s[i]=last_s[i]+1; last_s[i]=xxx;
			first_b[i]=length-first_b[i]-1;
			last_b[i]=length-last_b[i]-1;
			xxx=first_b[i]; first_b[i]=last_b[i]; last_b[i]=xxx;
		}
		//AlListLength=MolAlList.size();

		for (int i=0; i<length/2; i++) {
			xxx=Gnr[i]; Gnr[i]=Gnr[length-1-i]; Gnr[length-1-i]=xxx;
			xxx=n_mon[i]; n_mon[i]=n_mon[length-1-i]; n_mon[length-1-i]=xxx;
			xxx=mon_nr[i]; mon_nr[i]=mon_nr[length-1-i]; mon_nr[length-1-i]=xxx;
			//for (int j=0; j<AlListLength; j++) {
			//	xxx=Al[j]->frag[i]; Al[j]->frag[i]=Al[j]->frag[length-1-i]; Al[j]->frag[length-1-i]=xxx;
			//} //dit geeft problemen in valgrind. Als alias verbetert is dan moet dit ook weer aangepakt.
		}
		//for (int i=1; i<length-1; i++) { //needs a fix
		//		if ((Gnr[i-1]!=Gnr[i]) && (Gnr[i]!=Gnr[i+1])&& Gnr[i-1]!=Gnr[i+1]) {
		//			cout <<"Violation of the architecture in branched molecule '" + name + "': spacers between branches must be more than one segment in length.  " << endl;
	        //			return false;
		//		}
		//}

	}

	success=MakeMonList();
	if (chainlength==1 && MolType!=water) MolType=monomer;
	return success;
}

int Molecule::GetChainlength(void){
if (debug) cout <<"GetChainlength for Mol " + name << endl;
	return chainlength;
}

bool Molecule:: MakeMonList(void) {
if (debug) cout <<"Molecule:: MakeMonList" << endl;
	MolMonList.clear();
	bool success=true;
	int length = mon_nr.size();
	int i=0;
	while (i<length) {
		if (!In[0]->InSet(MolMonList,mon_nr[i])) {
			if (Seg[mon_nr[i]]->GetFreedom()=="frozen") {
				success = false;
				cout << "In 'composition of mol " + name + ", a segment was found with freedom 'frozen'. This is not permitted. " << endl;

			}
			MolMonList.push_back(mon_nr[i]);
		}
		i++;
	}
	i=0;
	int pos;
	while (i<length) {
		if (In[0]->InSet(MolMonList,pos,mon_nr[i])) {molmon_nr.push_back(pos);
		//cout << "in frag i " << i << " there is segment nr " <<  mon_nr[i] << " and it is on molmon  pos " << pos << endl;
		} else {cout <<"program error in mol PrepareForCalcualations" << endl; }
		i++;
	}
	return success;
}

bool Molecule::IsClamped() {
if (debug) cout <<"IsClamped for Mol " + name << endl;
	bool success=false;
	int length=mon_nr.size();
	// In case of failed composition reading, mon_nr will be empty and mon_nr[0] will segfault, hence the empty check.
	if (not mon_nr.empty() and mon_nr[0] == mon_nr[length-1]) {
		if (Seg[mon_nr[0]]->freedom == "clamp") success=true;
	}
	if (success && MolType!=linear) {
		success=false;
		cout <<"Sorry, currently clamped molecules should be linear" << endl;
	}
	return success;
}

bool Molecule::IsPinned() {
if (debug) cout <<"IsPinned for Mol " + name << endl;
	bool success=false;
	int length=MolMonList.size();
	int i=0;
	while (i<length) {
        	if (Seg[MolMonList[i]]->GetFreedom()=="pinned") success=true;
		i++;
	}
	return success;
}

int Molecule::GetPinnedSeg() {
if (debug) cout <<"GetPinnedSeg for Mol " + name << endl;
	int segnr=-1;
	int length=MolMonList.size();
	int i=0;
	while (i<length) {
        	if (Seg[MolMonList[i]]->GetFreedom()=="pinned") {
				if (segnr > -1) cout <<"multiple segnrs in GetPinnedSeg() " << endl; else segnr=MolMonList[i];
			}
		i++;
	}
	return segnr;
}

bool Molecule::IsTagged() {
if (debug) cout <<"IsTagged for Mol " + name << endl;
	bool success=false;
	int length=MolMonList.size();
	int i=0;
	while (i<length) {
		if (Seg[MolMonList[i]]->freedom=="tagged") {success = true; tag_segment=MolMonList[i]; }
		i++;
	}
	return success;
}

Real Molecule::Charge() {
if (debug) cout <<"Molecule:: Charge" << endl;
	Real charge=0;
	int length=mon_nr.size();
	int length_states;
	for (int i=0; i<length; i++) {
		if (Seg[mon_nr[i]]->state_name.size() >1) {
			length_states=Seg[mon_nr[i]]->state_name.size();
			for (int j=0; j<length_states; j++) charge +=Seg[mon_nr[i]]->state_alphabulk[j]*Seg[mon_nr[i]]->state_valence[j]*n_mon[i];
		} else	charge +=Seg[mon_nr[i]]->valence*n_mon[i];
	}
//cout <<"mol " << name << "charge " << charge << endl;
	return charge/chainlength;
}

bool Molecule::IsCharged() {
if (debug) cout <<"IsCharged for Mol " + name << endl;
	Real charge =0;
	bool ischarged=false;
	int length = n_mon.size();
	int length_states;
	int i=0;
	while (i<length) {
		if (Seg[mon_nr[i]]->state_name.size()>0) {
			length_states=Seg[mon_nr[i]]->state_name.size();
			for (int j=0; j<length_states; j++) {if (Seg[mon_nr[i]]->state_valence[j]!=0) ischarged=true; }
		} else
		charge +=n_mon[i]*Seg[mon_nr[i]]->valence;
		i++;
	}
	if (charge !=0) ischarged=true;
	return ischarged;
}

void Molecule::PutParameter(string new_param) {
if (debug) cout <<"PutParameter for Mol " + name << endl;
	KEYS.push_back(new_param);
}

string Molecule::GetValue(string parameter) {
if (debug) cout <<"GetValue " + parameter + " for Mol " + name << endl;
	int length = PARAMETERS.size();
	int i=0;
	while (i<length) {
		if (PARAMETERS[i]==parameter) { return VALUES[i];}
		i++;
	}
	return "";
}

void Molecule::push(string s, Real X) {
if (debug) cout <<"push (Real) for Mol " + name << endl;
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void Molecule::push(string s, int X) {
if (debug) cout <<"push (int) for Mol " + name << endl;
	ints.push_back(s);
	ints_value.push_back(X);
}
void Molecule::push(string s, bool X) {
if (debug) cout <<"push (bool) for Mol " + name << endl;
	bools.push_back(s);
	bools_value.push_back(X);
}
void Molecule::push(string s, string X) {
if (debug) cout <<"push (string) for Mol " + name << endl;
	strings.push_back(s);
	strings_value.push_back(X);
}



void Molecule::PushOutput() {
if (debug) cout <<"PushOutput for Mol " + name << endl;
	int length_al=MolAlList.size();
	for (int i=0; i<length_al; i++) Al[i]->PushOutput();
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	push("composition",GetValue("composition"));
	if (IsTagged()) {string s="tagged"; push("freedom",s);} else {push("freedom",freedom);}
	if (freedom=="free") theta = Lat[0]->WeightedSum(phitot);
	push("Markov",Markov);
	push("k_stiff",k_stiff);
	if (Lat[0]->gradients==3) {
		int MZ=Lat[0]->MZ;
		int MY=Lat[0]->MY;
		int MX=Lat[0]->MX;
		int JX=Lat[0]->JX;
		int JY=Lat[0]->JY;
		for (int z=1; z<MZ+1; z++) {
			Real phiz=0;
			for (int x=1; x<MX+1; x++) for (int y=1;y<MY+1;y++) {
				phiz +=phitot[x*JX+y*JY+z];
			}
			phiz /= MX*MY;
			if (z==1) push("phiz[1]",phiz);
			if (z==2) push("phiz[2]",phiz);
			if (z==3) push("phiz[3]",phiz);
			if (z==4) push("phiz[4]",phiz);
			if (z==5) push("phiz[5]",phiz);
			if (z==6) push("phiz[6]",phiz);
			if (z==7) push("phiz[7]",phiz);
			if (z==8) push("phiz[8]",phiz);
			if (z==9) push("phiz[9]",phiz);
			if (z==10) push("phiz[10]",phiz);
			if (z==11) push("phiz[11]",phiz);
			if (z==12) push("phiz[12]",phiz);
			if (z==13) push("phiz[13]",phiz);
			if (z==14) push("phiz[14]",phiz);
			if (z==15) push("phiz[15]",phiz);
			if (z==16) push("phiz[16]",phiz);
			if (z==17) push("phiz[17]",phiz);
			if (z==18) push("phiz[18]",phiz);
			if (z==19) push("phiz[19]",phiz);
			if (z==20) push("phiz[20]",phiz);
		}
	}
	if (Markov==2) {
		for (int k=0; k<size; k++){
			if (k==0) push("P[0]",P[0]);
			if (k==1) push("P[1]",P[1]);
			if (k==2) push("P[2]",P[2]);
			if (k==3) push("P[3]",P[3]);
			if (k==4) push("P[4]",P[4]);
		}
	}
	//if (theta==0) {
        //  if (Lat[0]->gradients==3) {
	//	    Lat[0]->remove_bounds(phitot);
	//	    Sum(theta,phitot,Lat[0]->M);
	//  }
	//}

	push("Rg",pow((Lat[0]->Moment(phitot,0.0,2)/chainlength),0.5));
	Lat[0]->remove_bounds(phitot);
	theta=Lat[0]->WeightedSum(phitot);
	push("theta",theta);
        //cout <<"theta " << name << " = " << theta << endl;
	Real thetaexc=theta-Lat[0]->volume*phibulk;
	push("theta_exc",thetaexc);
	push("thetaexc",thetaexc);
	push("theta_Gibbs",theta_Gibbs);
	if (R_Gibbs>0) push("R_Gibbs",R_Gibbs);
	push("n",n);
	push("chainlength",chainlength);
	push("phibulk",phibulk);
	if (MolType == water) {
		push("phib1",phib1);
		push("Kw",Kw);
	}
	push("Mu",Mu);
	if (Lat[0]->gradients==3) {
		Real TrueVolume=Lat[0]->MX*Lat[0]->MY*Lat[0]->MZ;
		Real Volume_particles=0;
		int num_of_seg=In[0]->MonList.size();
		for (int i=0; i<num_of_seg; i++) {
			if (Seg[i]->freedom=="frozen") Volume_particles += Seg[i]->Volume_particles();
		}
		push("Gamma",theta-(TrueVolume-Volume_particles)*phibulk);
	}
	if (GetValue("compute_width_interface").size()>0){
		if (!ComputeWidth()) {
			cout <<"Computation of width of interface is rejected" <<endl;
		}
	}
	push("width",width);
	push("phi1",phitot[Lat[0]->fjc]);
	push("phiM",phitot[Lat[0]->M-2*Lat[0]->fjc]);
	push("Dphi",phi1-phiM);
	push("pos_interface",pos_interface);
	push("phi_average",phi_av);
	if (chainlength==1) {
		int seg=MolMonList[0];
		if (Seg[seg]->ns >1) {
			if (mu_state.size() ==0) for (int i=0; i<Seg[seg]->ns; i++) mu_state.push_back(Mu);
			for (int i=0; i<Seg[seg]->ns; i++) {
				mu_state[i]+=log(Seg[seg]->state_alphabulk[i]);
				push("mu-"+Seg[seg]->state_name[i],mu_state[i]);
			}
		}
	}
	int M=Lat[0]->M;
	Real phimax=phitot[M/2];
	bool maxfound=false;
	int i=M/2;
	while (!maxfound) {
		i++;
		if (phitot[i]> phimax ) phimax =phitot[i]; else maxfound=true;
	}
	maxfound=false;
	i=M/2;
	while (!maxfound) {
		i--;
		if (phitot[i]> phimax ) phimax =phitot[i]; else maxfound=true;
	}

	push("phiMax",phimax);

	push("GN",GN);
	push("norm",norm);
	push("phi_LB_x",phi_LB_X);
	push("phi_UB_x",phi_UB_X);
	J=0;

	int molmonlength=MolMonList.size();
	for (int i=0; i<molmonlength; i++) {
		J+=Seg[MolMonList[i]]->J;
	}
	push("J",J/chainlength);
	push("DeltaMu",Delta_MU);
	string s="profile;0"; push("phi",s);
	int length = MolMonList.size();
	for (int i=0; i<length; i++) {
		stringstream ss; ss<<i+1; string str=ss.str();
		s= "profile;"+str; push("phi_"+Seg[MolMonList[i]]->name,s);
	}
	for (int i=0; i<length_al; i++) {
		push(Al[i]->name+"_value",Al[i]->value);
		push(Al[i]->name+"_composition",Al[i]->composition);
		stringstream ss; ss<<i+length; string str=ss.str();
		s="profile;"+str; push(Al[i]->name+"phi",s);
	}
	s="vector;0"; push("gn",s);
#ifdef CUDA
int M = Lat[0]->M;
	TransferDataToHost(H_phitot,phitot,M);
	TransferDataToHost(H_phi,phi,M*MolMonList.size());
#endif
}

Real* Molecule::GetPointer(string s, int &SIZE) {
if (debug) cout <<"GetPointer for Mol " + name << endl;
	vector<string> sub;
	int M= Lat[0]->M;
	In[0]->split(s,';',sub);
	if (sub[0]=="profile") {
		SIZE=M;

		if (sub[1]=="0") {
			Lat[0]->set_bounds(phitot);
			return H_phitot;
		}

		int length=MolMonList.size();
		int i=0;
		while (i<length) {
			stringstream ss; ss<<i+1; string str=ss.str();
			if (sub[1]==str) {
				Lat[0]->set_bounds(phi+i*M);
				return H_phi+i*M;
			}
			i++;
		}
		int length_al=MolAlList.size();
		i=0;
		while (i<length_al) {
			stringstream ss; ss<<i+length; string str=ss.str();
			if (sub[i]==str)  {
				Lat[0]->set_bounds(Al[i]->H_phi);
				return Al[i]->H_phi;
			}
		}
	} else { //sub[0]=="vector";
		if (sub[1]=="0") {SIZE=n_box; return gn;}
	}
	return NULL;
}
int* Molecule::GetPointerInt(string s, int &SIZE) {
if (debug) cout <<"GetPointerInt for Mol " + name << endl;
	vector<string> sub;
	In[0]->split(s,';',sub);
	if (sub[0]=="array") { //set SIZE and return array pointer.
	}
	return NULL;
}
int Molecule::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
if (debug) cout <<"GetValue (long) for Mol " + name << endl;
	int i=0;
	int length = ints.size();
	while (i<length) {
		if (prop==ints[i]) {
			int_result=ints_value[i];
			return 1;
		}
		i++;
	}
	i=0;
	length = Reals.size();
	while (i<length) {
		if (prop==Reals[i]) {
			Real_result=Reals_value[i];
			return 2;
		}
		i++;
	}
	i=0;
	length = bools.size();
	while (i<length) {
		if (prop==bools[i]) {
			if (bools_value[i]) string_result="true"; else string_result="false";
			return 3;
		}
		i++;
	}
	i=0;
	length = strings.size();
	while (i<length) {
		if (prop==strings[i]) {
			string_result=strings_value[i];
			return 3;
		}
		i++;
	}
	return 0;
}

Real Molecule::ComputeGibbs(Real R_gibbs) {
if (debug) cout <<"ComputeGibbs for Mol " + name << endl;
	int fjc=Lat[0]->fjc;
	int M=Lat[0]->M;
	int gradients=Lat[0]->gradients;
	if (gradients>1) {cout <<"Error in ComputeGibbs; gadients not equal to 1 " << endl; return 0.0;}
	Lat[0]->remove_bounds(phitot);
	Real phi_low =phitot[fjc+1];
	Real phi_high=phitot[M-2*fjc-1];
	theta=Lat[0]->WeightedSum(phitot);
	Real theta_exc=theta-Lat[0]->volume*phibulk;
	//cout <<"for mol : " << name << " theta " << theta << " theta_exc " << theta_exc << endl;
	//cout <<"volume " <<Lat[0]->volume << endl;
	//cout <<"phibulk " << phibulk << endl;
	if (freedom=="solvent") {
		theta_Gibbs=0;
		R_Gibbs=theta_exc/(phi_low-phi_high);
		//cout <<"R_gibbs " << R_Gibbs << endl;
		return  R_Gibbs;
	} else {
		R_Gibbs=R_gibbs;
		theta_Gibbs=theta_exc-(phi_low-phi_high)*R_gibbs;
		//cout <<" philow = " << phi_low << endl;
		//cout <<" phihigh = " << phi_high << endl;
		//cout <<" theta_Gibbs " << theta_Gibbs << endl;
		return theta_Gibbs;
	}
}

bool Molecule::ComputeWidth() {
if (debug) cout <<"ComputeWidth for Mol " + name << endl;
	bool success=true;
	int M=Lat[0]->M;
	if (Lat[0]->gradients>1) {success=false; cout <<" Compute width of interface only in system with 'one-gradient'" << endl; return success; }
	if (GetValue("compute_width_interface")!="true") {
		cout <<"Interfacial width not computed because value for 'compute_width_interface' was not set to 'true'. " << endl; success=false; return success;
	} else {
		int fjc=Lat[0]->fjc;
		width=0;
		Dphi=0;
		phi_av=0;
		phi1=phitot[fjc]; phiM=phitot[M-2*fjc]; Dphi=phi1-phiM;
		if (Dphi*Dphi<1e-20) {cout << "There is no interface present. Width evaluation failed;" << endl; return success; }
		for (int x=fjc; x<M-fjc-1; x++) {
			if ((phitot[x]-phitot[x+1])/Dphi > width) {width = (phitot[x]-phitot[x+1])/Dphi; pos_interface=x+0.5; phi_av=(phitot[x]+phitot[x+1])/2;}
			}
	}
	if (width >0) width = 1.0/width/Lat[0]->fjc;
	pos_interface = (pos_interface)/Lat[0]->fjc;
	return success;
}

void Molecule::NormPerBlock(int split) {
	int MX=Lat[0]->MX/split;
	int MY=Lat[0]->MY/split;
	int MZ=Lat[0]->MZ/split;
	int JX=Lat[0]->JX;
	int JY=Lat[0]->JY;
	int M = Lat[0]->M;
	Real theta_block;
	int blocknr=-1;
	for(int i=0; i<split; i++)
	for(int j=0; j<split; j++)
	for(int k=0; k<split; k++){
		theta_block=0;
		blocknr++;
		for (int x=1; x<=MX; x++)
		for (int y=1; y<=MY; y++)
		for (int z=1; z<=MZ; z++) {
			theta_block+=phitot[(i*MX+x)*JX+(j*MY+y)*JY+(k*MZ+z)];
		}
		for (int x=1; x<=MX; x++)
		for (int y=1; y<=MY; y++)
		for (int z=1; z<=MZ; z++) {
			phitot[(i*MX+x)*JX+(j*MY+y)*JY+(k*MZ+z)]*=block[blocknr]/theta_block;
		        int length=MolMonList.size();
			for (int kk=0; kk<length; kk++)
			     phi[kk*M+(i*MX+x)*JX+(j*MY+y)*JY+(k*MZ+z)]*=block[blocknr]/theta_block;
		}
	}
	//cout <<"theta : " << theta << " theta_tot : " << theta_tot << endl;
	//int length=block.size();
	//for (int p=0; p<length; p++) cout << "block("<<p<<")= " << block[p] << endl;
}

void Molecule::SetThetaBlocks(int split) {
	int MX=Lat[0]->MX/split;
	int MY=Lat[0]->MY/split;
	int MZ=Lat[0]->MZ/split;
	int JX=Lat[0]->JX;
	int JY=Lat[0]->JY;
	Real theta_block;
	Real theta_tot=0;
	block.clear();
	for(int i=0; i<split; i++)
	for(int j=0; j<split; j++)
	for(int k=0; k<split; k++){
		theta_block=0;
		for (int x=1; x<=MX; x++)
		for (int y=1; y<=MY; y++)
		for (int z=1; z<=MZ; z++) {
			theta_block+=phitot[(i*MX+x)*JX+(j*MY+y)*JY+(k*MZ+z)];
		}
		theta_tot+=theta_block;
		block.push_back(theta_block);
	}
	//cout <<"theta : " << theta << " theta_tot : " << theta_tot << endl;
	//int length=block.size();
	//for (int p=0; p<length; p++) cout << "block("<<p<<")= " << block[p] << endl;
}

Real* Molecule::propagate_forward(Real* G1, int &s, int block, int generation, int M) {
if (debug) cout <<"1. propagate_forward for Mol " + name << endl;

	int N= n_mon[block];

	if (save_memory) {
		int k,k0,t0,v0,t;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];

		t=1;
		v0=t0=k0=0;

		if (s==first_s[generation]) {

			//Cp(Gs+M,G1,M); Cp(Gs,G1,M);
			Lat[0]->Initiate(Gs+M,G1,Markov,M);
			Lat[0]->Initiate(Gs,G1,Markov,M); //not sure why this is done....
		} else {
			Lat[0] ->propagate(Gs,G1,0,1,M); //assuming Gs contains previous end-point distribution on pos zero;
		}
		Cp(Gg_f+n0*M,Gs+M,M); last_stored[block]=n0;
		s++;
		t=1;
		v0=t0=k0=0;
		for (k=2; k<=N; k++) {
			t++; s++;
			Lat[0]->propagate(Gs,G1,(k-1)%2,k%2,M);
			if (t>n) {
				t0++;
				if (t0 == n) t0 = ++v0;
				t = t0 + 1;
				k0 = k - t0 - 1;
			}

			if ((t == t0+1 && t0 == v0)
		  	 || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(k0+t0)))
		  	 || (2*(n-t+k) >= N-1)) {
				Cp(Gg_f+(n0+t-1)*M,Gs+(k%2)*M,M);
				last_stored[block]=n0+t-1;
			}
		}
		if ((N)%2!=0) {
			Cp(Gs,Gs+M,M);
		}
	} else {
		for (int k=0; k<N; k++) {
			if (s>first_s[generation]) {

				Lat[0] ->propagate(Gg_f,G1,s-1,s,M);
			} else {
				Lat[0]->Initiate(Gg_f+first_s[generation]*M,G1,Markov,M);
			}
			 s++;
		}
	}
	if (save_memory) {
		return Gg_f+last_stored[block]*M;
	} else {
		 return Gg_f+(s-1)*M;
	}

}

void Molecule::propagate_backward(Real* G1, int &s, int block, int unity, int M) {
if (debug) cout <<"propagate_backward for Mol " + name << endl;

	int N= n_mon[block];
	if (save_memory) {
		int k,k0,t0,v0,t,rk1;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];

		t=1;
		v0=t0=k0=0;
		for (k=2; k<=N; k++) {t++; if (t>n) { t0++; if (t0 == n) t0 = ++v0; t = t0 + 1; k0 = k - t0 - 1;}}
		for (k=N; k>=1; k--) {
			if (k==N) {
				if (s==chainlength-1) {
					Cp(Gg_b+(k%2)*M,G1,M);
				} else {
					Lat[0]->propagate(Gg_b,G1,(k+1)%2,k%2,M);
				}
			} else {
				Lat[0]->propagate(Gg_b,G1,(k+1)%2,k%2,M);
			}
			t = k - k0;

			if (t == t0) {
				k0 += - n + t0;
				if (t0 == v0 ) {
					k0 -= ((n - t0)*(n - t0 + 1))/2;
				}
				t0 --;
				if (t0 < v0) {
					v0 = t0;
				}
				Cp(Gs+(t%2)*M,Gg_f+(n0+t-1)*M,M);
				for (rk1=k0+t0+2; rk1<=k; rk1++) {
					t++;
					Lat[0]->propagate(Gs,G1,(t-1)%2,t%2,M);
					if (t == t0+1 || k0+n == k) {
						Cp(Gg_f+(n0+t-1)*M,Gs+(t%2)*M,M);
					}
					if (t == n && k0+n < k) {
						t  = ++t0;
						k0 += n - t0;
					}
				}
				t = n;
			}
			Lat[0]->AddPhiS(rho+molmon_nr[block]*M,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,Markov,M);

			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) {
						Lat[0]->AddPhiS(Al[i]->rho,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,G1,norm,Markov,M);
					}
				}
			}
			s--;
		}
		Cp(Gg_b,Gg_b+M,M);
	} else {
		for (int k=0; k<N; k++) {
			if (s<chainlength-1) {
				Lat[0]->propagate(Gg_b,G1,(s+1)%2,s%2,M);
			} else {
				Lat[0]->Initiate(Gg_b+(s%2)*M,G1,Markov,M);
			}

			Lat[0]->AddPhiS(rho+molmon_nr[block]*M, Gg_f+(s*M), Gg_b+(s%2)*M,Markov, M);

			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) {
						//  Composition(Al[i]->rho,Gg_f+s*M,Gg_b+(s%2)*M,G1,norm,M);
						Lat[0]->AddPhiS(Al[i]->rho,Gg_f+s*M,Gg_b+(s%2)*M,G1,norm,Markov,M);
					}
				}
			}
			s--;
		}
	}
}



Real* Molecule::propagate_forward(Real* G1, int &s, int block, Real* P, int generation, int M) {
if (debug) cout <<"1. propagate_forward for Mol " + name << endl;

	int N= n_mon[block];
	if (save_memory) {
		int k,k0,t0,v0,t;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];
		if (s==first_s[generation]) {

			Lat[0]->Initiate(Gs+size*M,G1,Markov,M);
			//Lat[0]->Initiate(Gs,G1,Markov,M); //not necessary.
		} else {
			Lat[0] ->propagateF(Gs,G1,P,0,1,M); //assuming Gs contains previous end-point distribution on pos zero;

		}
		s++;
		Cp(Gg_f+n0*M*size,Gs+M*size,M*size);
		last_stored[block]=n0;

		t=1;
		v0=t0=k0=0;
		for (k=2; k<=N; k++) {
			t++; s++;
			Lat[0]->propagateF(Gs,G1,P,(k-1)%2,k%2,M);
			if (t>n) {
				t0++;
				if (t0 == n) { t0 = ++v0; cout <<"v0 >0 .... save memory may fail!" << endl; }
				t = t0 + 1;
				k0 = k - t0 - 1;
			}
			if ((t == t0+1 && t0 == v0)
		  	 || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(k0+t0)))
		  	 || (2*(n-t+k) >= N-1)) {
				Cp(Gg_f+(n0+t-1)*M*size,Gs+(k%2)*M*size,M*size);
				last_stored[block]=n0+t-1;
			}
		}
		if ((N)%2!=0) {
			Cp(Gs,Gs+M*size,M*size); //FL
		}
	} else {
		for (int k=0; k<N; k++) {
			if (s>first_s[generation]) {
				Lat[0] ->propagateF(Gg_f,G1,P,s-1,s,M);
			} else {
				//Cp(Gg_f+first_s[generation]*M,G1,M);
				Lat[0]->Initiate(Gg_f+first_s[generation]*M*size,G1,Markov,M);
			}
			 s++;
		}
	}
	if (save_memory) {
		return Gg_f+last_stored[block]*M*size;
	} else {
		 return Gg_f+(s-1)*M*size;
	}

}

void Molecule::propagate_backward(Real* G1, int &s, int block, Real* P, int& unity, int M) {
if (debug) cout <<"propagate_backward for Mol " + name << endl;
	int N= n_mon[block];
	if (save_memory) {
		int k,k0,t0,v0,t,rk1;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];

		t=1;
		v0=t0=k0=0;
		for (k=2; k<=N; k++) {t++; if (t>n) { t0++; if (t0 == n) t0 = ++v0; t = t0 + 1; k0 = k - t0 - 1;}}
		for (k=N; k>=1; k--) {
			if (k==N) {
				if (s==chainlength-1) {
					Lat[0]->Initiate(Gg_b+(k%2)*M*size,G1,Markov,M);
					//Cp(Gg_b+(k%2)*M,G1,M);
				} else {
					if (unity==-1) {
						unity=0;
						//cout <<"SM unity " << s << endl;
						Real* GB= (Real*) malloc(2*M*sizeof(Real)); //must be adjusted for cuda
						Lat[0]->Terminate(GB,Gg_b+((k+1)%2)*M*size,Markov,M);
						Lat[0]->propagate(GB,G1,0,1,M); //first step is freely joined
						Lat[0]->Initiate(Gg_b+(k%2)*M*size,GB+M,Markov,M);
						free(GB);
					} else {
						Lat[0]->propagateB(Gg_b,G1,P,(k+1)%2,k%2,M); //FL
					}
				}
			} else {
				Lat[0]->propagateB(Gg_b,G1,P,(k+1)%2,k%2,M); //FL
			}
			t = k - k0;
			if (t == t0) {
				k0 += - n + t0;
				if (t0 == v0 ) {
					k0 -= ((n - t0)*(n - t0 + 1))/2;
				}
				t0 --;
				if (t0 < v0) {
					v0 = t0;
				}
				Cp(Gs+(t%2)*M*size,Gg_f+(n0+t-1)*M*size,M*size); //FL
				for (rk1=k0+t0+2; rk1<=k; rk1++) {
					t++;
					Lat[0]->propagateF(Gs,G1,P,(t-1)%2,t%2,M);
					if (t == t0+1 || k0+n == k) {
						Cp(Gg_f+(n0+t-1)*M*size,Gs+(t%2)*M*size,M*size); //FL
					}
					if (t == n && k0+n < k) {
						t  = ++t0;
						k0 += n - t0;
					}
				}
				t = n;
			}

			//AddTimes(rho+molmon_nr[block]*M,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,M);
			Lat[0]->AddPhiS(rho+molmon_nr[block]*M,Gg_f+(n0+t-1)*M*size,Gg_b+(k%2)*M*size,Markov,M);

			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) {
						//Composition(Al[i]->rho,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,G1,norm,M);
						Lat[0]->AddPhiS(Al[i]->rho,Gg_f+(n0+t-1)*M*size,Gg_b+(k%2)*M*size,G1,norm,Markov,M);
					}
				}
			}
			s--;
		}
		Cp(Gg_b,Gg_b+M*size,size*M);//FL
	} else {
		for (int k=0; k<N; k++) {
			if (s<chainlength-1) {
				if (unity==-1) {
					unity=0;
					//cout <<"unity " << s << endl;
					Real* GB= (Real*) malloc(2*M*sizeof(Real)); //must be adjusted for cuda
					Lat[0]->Terminate(GB,Gg_b+((s+1)%2)*M*size,Markov,M);
					Lat[0]->propagate(GB,G1,0,1,M); //first step is freely joined
					Lat[0]->Initiate(Gg_b+(s%2)*M*size,GB+M,Markov,M);
					free(GB);
				} else {
					Lat[0]->propagateB(Gg_b,G1,P,(s+1)%2,s%2,M); //FL
				}
			} else {
				Lat[0]->Initiate(Gg_b+(s%2)*M*size,G1,Markov,M); //FL
			}

			Lat[0]->AddPhiS(rho+molmon_nr[block]*M, Gg_f+s*M*size, Gg_b+(s%2)*M*size, Markov, M);


			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) {
						//Composition(Al[i]->rho,Gg_f+s*M,Gg_b+(s%2)*M,G1,norm,M);
						Lat[0]->AddPhiS(Al[i]->rho,Gg_f+s*M*size,Gg_b+(s%2)*M*size,G1,norm,Markov, M);
					}
				}
			}
			s--;
		}
	}
}

bool Molecule::ComputePhi(Real* BETA,int id){
if (debug) cout <<"ComputePhi for Mol " + name << endl;
	bool success=true;
	int M=Lat[0]->M;
	Lat[0]->sub_box_on=0;//selecting 'standard' boundary condition
	if (id !=0) {
		int molmonlistlength= MolMonList.size();
		for (int i=0; i<molmonlistlength; i++)
		if (id==1) {
			Times(Seg[MolMonList[i]]->G1,Seg[MolMonList[i]]->G1,BETA,M);
		} else {
			Div(Seg[MolMonList[i]]->G1,BETA,M);
		}
	}
	if (MolType==water) phib1=0;
	if (freedom == "clamped") {
		Lat[0]->sub_box_on=Seg[mon_nr[0]]->clamp_nr; //slecting sub_box boundary conditions.
				//success=ComputePhiLin();
		success=ComputePhi();
		Lat[0]->sub_box_on=0;//selecting 'standard' boundary condition
	} else {
		success=ComputePhi();
	}


	if (id !=0) {
		int molmonlistlength=MolMonList.size();

		for (int i=0; i<molmonlistlength; i++) {
			if (id==1) {
				Div(phi+i*M,BETA,M);
				Div(Seg[MolMonList[i]]->G1,BETA,M);
			}  else {
				Times(phi+i*M,phi+i*M,BETA,M);
				Times(Seg[MolMonList[i]]->G1,Seg[MolMonList[i]]->G1,BETA,M);
			}
		}
	}
	return success;
}

Real Molecule::GetPhib1() {
	return 0;
}

void Molecule::AddToGP(Real*) {
}

void Molecule::AddToF(Real*) {
}


bool Molecule::ComputePhi(){
if (debug) cout <<"ComputePhi for Molecule " + name << endl; //default computation for monomer only....
	int M=Lat[0]->M;
	bool success=true;
	Cp(phi,Seg[mon_nr[0]]->G1,M);
	GN=Lat[0]->WeightedSum(phi);
	if (compute_phi_alias)
		for (auto& alias : Al) //For every alias in the Al vector (same as Al[i])
			if (alias->frag[0]==1) {
				Cp(alias->phi,phi,M);
				Norm(alias->phi,norm,M);
			}
	Times(phi,phi,Seg[mon_nr[0]]->G1,M);
	return success;
}


Real Molecule::fraction(int segnr){
if (debug) cout <<"fraction for mol_test " + name << endl; //default for monomer.
	int Nseg=0;
	int length = mon_nr.size();
	int i=0;
	if (ring) i++; //first segment is not counted in fraction;
	while (i<length) {
		if (segnr==mon_nr[i]) {Nseg+=n_mon[i];}
		i++;
	}
	return 1.0*Nseg/chainlength;
}


