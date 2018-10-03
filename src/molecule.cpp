#include "molecule.h"


Molecule::Molecule(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) {
	In=In_; Seg=Seg_; name=name_;  Lat=Lat_;
if (debug) cout <<"Constructor for Mol " + name << endl;
	KEYS.push_back("freedom");
	KEYS.push_back("composition");
	KEYS.push_back("theta");
	KEYS.push_back("phibulk");
	KEYS.push_back("n");
	KEYS.push_back("save_memory");
	KEYS.push_back("restricted_range");
}

Molecule::~Molecule() {
	DeAllocateMemory();
}

void Molecule :: DeAllocateMemory(){
if (debug) cout <<"Destructor for Mol " + name << endl;
	if (H_phi!=NULL) free(H_phi);
	free(H_phitot);
	//free(H_u);
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
		//cudaFree(u);
		//cudaFree(G1);
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
	//free(G1);
	free(UNITY);
	free(Gg_f);
	free(Gg_b);
	if (save_memory) free(Gs);
#endif
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
	int M=Lat[0]->M;
	int m=0;
	if (freedom=="clamped") {
		m=Lat[0]->m[Seg[mon_nr[0]]->clamp_nr];
		n_box = Seg[mon_nr[0]]->n_box;
	}

	if (save_memory) {
		int length_ = mon_nr.size();
		for (int i=0; i<length_; i++) last_stored.push_back(0);
		for (int i=0; i<length_; i++) {
			int n=int(pow(n_mon[i]*6,1.0/3)+0.5);
			if (n_mon[i]<120) n++;
			if (n_mon[i]<60) n++;
			if (n>n_mon[i]) n=n_mon[i]; //This seems to me to be enough...needs a check though..
			if (i==0) memory.push_back(n); else memory.push_back(n+memory[i-1]);
		}
	}
	N=0;
	if (save_memory) {
		N=memory[n_mon.size()-1];
	} else {
		int length_ = mon_nr.size();
		for (int i=0; i<length_; i++) {N+=n_mon[i];}
		//cout <<"allocate memory for " << N << " EPD" << endl;
		//N=chainlength; //in case of dendrimers this is not correct. Way fewer EDF needed in this case.
	}

	H_phi = (Real*) malloc(M*MolMonList.size()*sizeof(Real)); Zero(H_phi,M*MolMonList.size());
//cout <<"molmonlist.size in mol" << MolMonList.size() << endl;
	//H_u = (Real*) malloc(M*MolMonList.size()*sizeof(Real));
	//Zero(H_u, M*MolMonList.size());
	H_phitot = (Real*) malloc(M*sizeof(Real)); Zero(H_phitot,M);
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
		rho=(Real*)AllOnDev(m*n_box*MolMonList.size()); Zero(rho,m*n_box,MolMonList.size());
		phi=(Real*)AllOnDev(M*MolMonList.size()); Zero(phi,M*MolMonList.size());
		if (save_memory) {Gs=(Real*)AllOnDev(m*n_box*2); Zero(Gs,m*n_box*2);}
	} else {
		Gg_f=(Real*)AllOnDev(M*N);
		Gg_b=(Real*)AllOnDev(M*2);
		phi=(Real*)AllOnDev(M*MolMonList.size());
		rho=phi;
		if (save_memory) Gs =(Real*)AllOnDev(M*2);
	}
	phitot=(Real*)AllOnDev(M);
	UNITY = (Real*)AllonDev(M);
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
		Gg_f = (Real*) malloc(M*N*sizeof(Real));
		Gg_b = (Real*) malloc(M*2*sizeof(Real));
		Zero(Gg_f,M*N);
		Zero(Gg_b,2*M);
		phi=H_phi;
		rho=phi;
		if (save_memory) Gs=(Real*) malloc(2*M*sizeof(Real));
	}
	//u = H_u;
	//G1 = (Real*)malloc(M*MolMonList.size()*sizeof(Real));
	phitot = H_phitot;
	UNITY = (Real*) malloc(M*sizeof(Real));
#endif

	if (save_memory) Zero(Gs,2*M);
	int length =MolAlList.size();
	if (freedom!="clamped") Seg[mon_nr[0]]->clamp_nr=0;
	for (int i=0; i<length; i++) Al[i]->AllocateMemory(Seg[mon_nr[0]]->clamp_nr,n_box);
}

bool Molecule:: PrepareForCalculations(int *KSAM) {
if (debug) cout <<"PrepareForCalculations in Mol " + name << endl;
int m=0;
if (freedom=="clamped") m=Lat[0]->m[Seg[mon_nr[0]]->clamp_nr];
int M=Lat[0]->M;
	if (freedom=="clamped") {
		Zero(H_mask1,n_box*m);
		Zero(H_mask2,n_box*m);
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
			H_mask1[i*m + jx*(Px1[i]-Bx[i])+jy*(Py1[i]-By[i])+(Pz1[i]-Bz[i])]=1;
			H_mask2[i*m + jx*(Px2[i]-Bx[i])+jy*(Py2[i]-By[i])+(Pz2[i]-Bz[i])]=1;
		}
#ifdef CUDA
		TransferDataToDevice(H_mask1,mask1,m*n_box);
		TransferDataToDevice(H_mask2,mask2,m*n_box);
		TransferDataToDevice(H_Bx,Bx,n_box);
		TransferDataToDevice(H_By,By,n_box);
		TransferDataToDevice(H_Bz,Bz,n_box);
		TransferDataToDevice(H_Px1,Px1,n_box);
		TransferDataToDevice(H_Py1,Py1,n_box);
		TransferDataToDevice(H_Pz1,Pz1,n_box);
		TransferDataToDevice(H_Px2,Px2,n_box);
		TransferDataToDevice(H_Py2,Py2,n_box);
		TransferDataToDevice(H_Pz2,Pz2,n_box);
		TransferDataToDevice(H_u, u, MolMonList.size()*M);
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
			pathlength_even=(Px2[i]-Px1[i]+Py2[i]-Py1[i]+Pz2[i]-Pz1[i])%2;
			if (chainlength_even == pathlength_even)
			cout <<" Warning, for chain part " << i << " the paths between clamps is not commensurate with the length of the chain fragment. Consider moving one of the calmp point by one (more) site!" << endl;
		}
		n = n_box;
		theta=n_box*chainlength;
	}
	return success;
}

bool Molecule::CheckInput(int start_) {
if (debug) cout <<"Molecule:: CheckInput" << endl;
start=start_;
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
			if (!Decomposition(GetValue("composition"))) {cout << "For mol '" + name + "' the composition is rejected. " << endl; success=false;}

		}
		if (GetValue("restricted_range").size()>0) {
			if (GetValue("freedom")!="range_restricted") cout <<"For mol '" + name + "' freedom is not set to 'range_restricted' and therefore  the value of 'restricted_range' is ignored" << endl;
		}
		if ( IsClamped()) {
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
		if (GetValue("freedom").size()==0 && !IsTagged()) {
			cout <<"For mol " + name + " the setting 'freedom' is expected: options: 'free' 'restricted' 'solvent' 'neutralizer' 'range_restricted' 'clamped' 'tagged' . Problem terminated " << endl; success = false;
			} else { if (!IsTagged()) {
				vector<string> free_list;
				if (!IsPinned()) {
					free_list.push_back("free");
					free_list.push_back("solvent");
					free_list.push_back("neutralizer");
					free_list.push_back("range_restricted");
				}
				free_list.push_back("restricted");
				if (!In[0]->Get_string(GetValue("freedom"),freedom,free_list,"In mol " + name + " the value for 'freedom' is not recognised ")) success=false;
				if (freedom == "solvent") {
					if (IsPinned()) {success=false; cout << "Mol '" + name + "' is 'pinned' and therefore this molecule can not be the solvent" << endl; }
				}
				if (freedom == "neutralizer") {
					if (IsPinned()) {success=false; cout << "Mol '" + name + "' is 'pinned' and therefore this molecule can not be the neutralizer" << endl; }
					if (!IsCharged()) {success=false; cout << "Mol '" + name + "' is not 'charged' and therefore this molecule can not be the neutralizer" << endl; }
					if (IsClamped()) {success=false; cout <<"Mol '" + name + "' is 'clamped' and therefore this molecule can not be the neutralizer" << endl;}
				}
				if (freedom == "free") {
					if (GetValue("theta").size()>0 || GetValue("n").size() > 0) {
						cout << "In mol " + name + ", the setting of 'freedom = free', can not not be combined with 'theta' or 'n': use 'phibulk' instead." << endl; success=false;
					} else {
						if (GetValue("phibulk").size() ==0) {
							cout <<"In mol " + name + ", the setting 'freedom = free' should be combined with a value for 'phibulk'. "<<endl; success=false;
						} else {

					phibulk=In[0]->Get_Real(GetValue("phibulk"),-1);
							if (phibulk < 0 || phibulk >1) {
								cout << "In mol " + name + ", the value of 'phibulk' is out of range 0 .. 1." << endl; success=false;
							}
						}
					}
				}
				if (freedom == "restricted" || freedom=="range_restricted") {
					if (GetValue("phibulk").size()>0) {
						cout << "In mol " + name + ", the setting of 'freedom = restricted' or 'freedom = range_restricted', can not not be combined with 'phibulk'  use 'theta' or 'n'  instead." << endl; success=false;
					} else {
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
						success=Lat[0]->ReadRange(r,HP,npos,block,GetValue("restricted_range"),name,s);
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
	}
	return success;
}

bool Molecule::PutVarInfo(string Var_type_,string Var_target_,Real Var_target_value_){
if (debug) cout <<"Molecule:: PutVarInfo" << endl;
	bool success=true;
	Var_scan_value=-1;
	Var_search_value=-1;
	vector<string>sub;
	Var_target=-1;
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
		if (Var_target_=="theta") {Var_search_value=0; Var_start_search_value=theta; }
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
		num_of_steps=(Var_end_value-Var_start_value)/step+1;

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
				theta=pow(10,(1-step_nr/num_of_steps)*log10(Var_start_value)+(step_nr/num_of_steps)*log10(Var_end_value));
			} else {
				theta=Var_start_value+step_nr*Var_step;
			}
			n=theta/chainlength;
			cout <<"mol : " + name + " : theta : " << theta << endl;
			break;
		case 1:
			if (scale=="exponential") {
				n=pow(10,(1-step_nr/num_of_steps)*log10(Var_start_value)+(step_nr/num_of_steps)*log10(Var_end_value));
			} else {
				n=Var_start_value+step_nr*Var_step;
			}
			cout <<"mol : " + name + " : n : " << n << endl;
			theta=n*chainlength;
			break;
		case 2:
			if (scale=="exponential") {
				phibulk=pow(10,(1-step_nr/num_of_steps)*log10(Var_start_value)+(step_nr/num_of_steps)*log10(Var_end_value));
			} else {
				phibulk=Var_start_value+step_nr*Var_step;
			}
			cout <<"mol : " + name + " : phibulk " << phibulk << endl;
			break;
		case 3:
			if (scale=="exponential") {
				Al[var_al_nr]->value =(int)pow(10,(1-step_nr/num_of_steps)*log10(Var_start_value) + (step_nr/num_of_steps)*log10(Var_end_value));
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
			break;
		case 3:
			if (Var_target_value!=0) Error=Mu/Var_target_value-1.0; else Error = Mu;
			break;
		default:
			//cout <<"program error in Molecule::GetError" << endl;
			break;
	}

	if (Var_search_value==3) {
		Error=theta;
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
			cout <<"In composition of mol '" + name + "' Alias '" + sA + "' was not found"<<endl; success=false;
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
	return success;
}



bool Molecule::ExpandBrackets(string &s) {
if (debug) cout <<"Molecule:: ExpandBrackets" << endl;
	bool success=true;
	vector<int> open;
	vector<int> close;
	bool done=false; //now interpreted the (expanded) composition
	while (!done) { done = true;
		open.clear(); close.clear();
		if (!In[0]->EvenBrackets(s,open,close)) {
			cout << "s : " << s << endl;
			cout << "In composition of mol '" + name + "' the backets are not balanced."<<endl; success=false;
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
					int x=In[0]->Get_int(s.substr(pos_close+1),0);
					string sA,sB,sC;
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
	return success;
}

bool Molecule::Interpret(string s,int generation){
if (debug) cout <<"Molecule:: Interpret" << endl;
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
				if (mnr <0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised"  << endl; success=false;
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
				if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity " << endl; success=false;
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
	while  (pos_open<pos_close) {
		pos_open =s.length()+1;
		pos_close=s.length();
		openfound=closedfound=false;
		i=0;
		while (i<length && !(openfound && closedfound) ){
			if (close[i]>pos && !closedfound) {closedfound=true; pos_close=close[i]+1; new_generation=i+1;}
			if (open[i]>pos && !openfound) {openfound=true; pos_open=open[i]+1; newgeneration=i+1;}
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
				GenerateTree(s,new_generation,pos,open,close);
				pos_close=pos_open+1;
			} else {
				pos=pos_close;
				Interpret(ss,generation);
			}
		} else {
			ss=s.substr(pos,pos_open-pos);
			pos=pos_open;
			Interpret(ss,generation);
			first_s.push_back(-1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			GenerateTree(s,newgeneration,pos,open,close);
		}
	}
	return success;
}

bool Molecule::GenerateDend(string s, int generation) {
if (debug) cout <<"Molecule:: GenerateDend" << endl;
return true;
}



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
		if (SUB[0]=="dend") {
			MolType=dendrimer; keyfound=true;
			s=s.substr(6,s.length()-7);
		}
		if (SUB[0]=="comb") {
			MolType=comb; keyfound=true;
			s=s.substr(6,s.length()-7);
		}
		if (SUB[0]=="ring") {
			MolType=ring; keyfound=true;
			s=s.substr(6,s.length()-7);
		}
		if (!keyfound) { success=false; cout << "Keyword specifying Moltype not recognised: select from @dend, @comb, @ring. Problem terminated "<< endl ;
			return false;
		 }
	}

	while (aliases && success) {
		loopnr++;
		sub.clear();
		In[0]->split(s,'#',sub);
		if (sub.size()%2!=1) {cout << " Alias in composition should be bracketed on both sides by '#'. For example, (A)10#alias_name#(B)3, when e.g., 'alias : alias_name : value : (X)5' is defined, so that you oubtain (A)10(X)5(B)3 " << endl; success=false; }
		aliases=(s!=sub[0]);
		if (aliases) ExpandAlias(sub,s);
		if (loopnr == 20) {
			cout << "Nesting nr 20 reached in aliases for mol " + name + " -composition. It is decided that this is too deep to continue; Possible, you have defined an alias-A inside alias-B which itself refers to alias-A. This is not allowed. Problem terminated. " << endl;
			success=false;
		}
	}
	sub.clear();
	In[0]->split(s,',',sub);
	int length=sub.size();
	int length_open;
	int i,j,k,a,f,dd;
	string ss;
	switch(MolType) {
		case dendrimer:
			for (i=0; i<length; i++) {
				open.clear(); close.clear();
				In[0]->EvenBrackets(sub[i],open,close);
				length_open=open.size();
				//if (sub[i].substr(0,1)=="(") {
				if (length_open>0) {
					if (!ExpandBrackets(sub[i])) cout <<"brackets '(' ')' not well positioned in "+sub[i] << endl;
          success=false;
				}
			}
			for (i=0; i<length; i++) {
				ss.append(sub[i]);
				if (i<length-1) ss.append(",");
			}
			s=ss;
		break;
		case comb:
			success=false;
		break;
		default:
			if (!ExpandBrackets(s)) success=false;
		break;
	}

	if (!In[0]->EvenSquareBrackets(s,open,close)) {

		cout << "Error in composition of mol '" + name + "'; the square brackets are not balanced. " << endl;
		success=false;
	}

	if (open.size()>0) {
		if (MolType==linear) { //overwrite default.
			MolType=branched;
		} else {
			success=false;
			cout <<" In 'composition' you can not combine special keywords such as 'dend, comb or ring' with branched compositions. This implies that square brackets are not allowed. " <<endl ;
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

	switch(MolType) {
		case linear:
			first_s.push_back(-1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			GenerateTree(s,generation,pos,open,close);
			break;
		case branched:
			first_s.push_back(-1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			GenerateTree(s,generation,pos,open,close);
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
			cout <<" s is " << s << endl;
			In[0]->split(s,';',sub_gen);
			n_generations=sub_gen.size();
			sym_dend=true; //default.
			cout <<"n_generations " << n_generations << endl;
			chainlength=0; N=-1;
			for (i=0; i<n_generations; i++) {
				sub.clear();	arms=0;
				In[0]->split(sub_gen[i],',',sub);
				if (sub.size()%2==0) {success=false; cout << "In composition of dend for generation " << i << " that is, in " + sub_gen[i] + " the number of arguments is not odd; use @dend(?) for help. " << endl; }
				if (sub.size()<2) { success=false;

					cout<<" ------Dendrimer language:------ " << endl;
					cout<<" example 1: consider segment name 'A' for the branch points and spacers (B)2 with functionality 3 and 4 generations: " << endl;
					cout<<" @dend(A,(B)2,3; A,(B)2,3; A,(B)2,3; A,(B)2,3)  " << endl;
					cout<<" hence generations are seperated by ';', branch points are 'monomer_names' without '(', ')'. Each generation can have unique numbers and structures. " << endl;
					cout<<" expample 2: asymmetric dendrimers, with asymmetry in branches of first generation:" << endl;
					cout <<" @dend(A,(B)2,3,(B)4,1; A,(B)2,3,(B)2,2; A,(B)2,3; A,(B)2,3)  " << endl;
					cout <<" example 3: star with 10 arms " << endl;
					cout <<" @dend(A,(B)100,10) " << endl;
					cout <<" example 4: 10 arm star end-grafted by one arm to a surface by way of segment 'X' " << endl;
					cout <<" @dend(A,(B)100,9,(B)99(X)1,1)" << endl;
					cout<<" ------Dendrimer language:------ " << endl;
				}


				length_g = sub.size(); //check nr generations
				if (length_g != 3) sym_dend=false;
				sub_dd.clear();
				In[0]->split(sub[0],':',sub_dd);
				length_dd =sub_dd.size();
//cout <<"length_dd: " << length_dd << endl;

				if (length_dd==4) {
					mnr=GetMonNr(sub_dd[2]);
					a=In[0]->Get_int(sub_dd[2],0);
					if (Al[a]->active) Al[a]->active=false; else Al[a]->active=true;

				} else mnr=GetMonNr(sub[0]);
				if (mnr <0)  {success=false; cout <<"In composition of mol '" + name + "', segment name '" + sub_dd[0] + "' is not recognised"  << endl; success=false; }
//cout <<"mon_nr is " << mnr << endl;
				for (a=0; a<AlListLength; a++) {if (Al[a]->active) Al[a]->frag.push_back(1); else Al[a]->frag.push_back(0);}

				n_mon.push_back(1);
				mon_nr.push_back(mnr);
				d_mon.push_back(degeneracy);
				first_a.push_back(-1);
				last_a.push_back(-1);
				k=1; chainlength+=degeneracy; N++;

				while (k<length_g-1) {
					arm++;
//cout <<"arm " << arm << endl;
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
//cout <<"length _dd = " <<length_dd << endl;
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
//cout <<"sub[k] " << sub[k] << endl;
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

//for (int i=0; i<n_generations; i++) cout<<"generation " << i << " first_a " << first_a[i] << " last_a " << last_a[i] << endl;
//length=n_arm.size();
//for (int i=0; i<length; i++) cout <<"arm " << i << " first_b " << first_b[i] << " last_b " << last_b[i] << endl;
//for (int i=0; i<length; i++) cout <<"arm " << i << " first_s " << first_s[i] << " last_s " << last_s[i] << endl;
//for (int i=0; i<length; i++) cout <<"arm " << i << " n_arm " << n_arm[i] << endl;
//length=n_mon.size();
//for (int i=0; i<length; i++) cout <<"block " << i << " n_mon " << n_mon[i] << " mon_nr " << mon_nr[i] << " d_mon " << d_mon[i] << endl;



			//success=false;
			break;
		case comb:
			cout <<"comb polymers not implemented" << endl; 	success=false;
			break;
		case ring:
			cout <<"ring polymers not implemented" << endl; 	success=false;
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

		for (int i=0; i<length/2; i++) {
			xxx=Gnr[i]; Gnr[i]=Gnr[length-1-i]; Gnr[length-1-i]=xxx;
			xxx=n_mon[i]; n_mon[i]=n_mon[length-1-i]; n_mon[length-1-i]=xxx;
			xxx=mon_nr[i]; mon_nr[i]=mon_nr[length-1-i]; mon_nr[length-1-i]=xxx;
			for (int j=0; j<AlListLength; j++) {
				xxx=Al[j]->frag[i]; Al[j]->frag[i]=Al[j]->frag[length-1-i]; Al[j]->frag[length-1-i]=xxx;
			}
		}

//	length=Gnr.size();
//	for (int i=0; i<length; i++)
//		cout << "segnr " << mon_nr[i] << " nn " << n_mon[i] << " Gnr " << Gnr[i] << endl;
//	length=last_s.size();
//	for (int i=0; i<length; i++)
//		cout << "Gen " << i  << " first_s " << first_s[i] << " last_s " << last_s[i] << " first_b " << first_b[i] << " last_b " <<  last_b[i] << endl;


	}

	if (MolType==dendrimer){
		int ChainLength=N;
		int length_a=first_b.size();
		int xxx;
		int length=mon_nr.size();

		for (int i=0; i<length_a; i++) {
			first_s[i] = ChainLength-first_s[i];
			last_s[i] = ChainLength-last_s[i];
			xxx=first_s[i]; first_s[i]=last_s[i]; last_s[i]=xxx;
			first_b[i]=length-first_b[i]-1;
			last_b[i]=length-last_b[i]-1;
			xxx=first_b[i]; first_b[i]=last_b[i]; last_b[i]=xxx;

		}
		for (int i=0; i<(length_a)/2; i++) {
			xxx=n_arm[i]; n_arm[i]=n_arm[length_a-i-1]; n_arm[length_a-i-1]=xxx;
			xxx=first_s[i]; first_s[i]=first_s[length_a-i-1]; first_s[length_a-i-1]=xxx;
			xxx=last_s[i]; last_s[i]=last_s[length_a-i-1]; last_s[length_a-i-1]=xxx;
			xxx=first_b[i]; first_b[i]=first_b[length_a-i-1]; first_b[length_a-i-1]=xxx;
			xxx=last_b[i]; last_b[i]=last_b[length_a-i-1]; last_b[length_a-i-1]=xxx;
		}
		length_a=first_a.size();
		length=first_b.size();
		for (int i=0; i<length_a; i++) {
			first_a[i]=length-first_a[i]-1;
			last_a[i]=length-last_a[i]-1;
			xxx=first_a[i]; first_a[i]=last_a[i]; last_a[i]=xxx;
		}
		for (int i=0; i<(length_a)/2; i++) {
			xxx=first_a[i]; first_a[i]=first_a[length_a-i-1]; first_a[length_a-i-1]=xxx;
			xxx=last_a[i]; last_a[i]=last_a[length_a-i-1]; last_a[length_a-i-1]=xxx;
		}
		length=n_mon.size();
		for (int i=0; i<length/2; i++) {
			xxx=n_mon[i]; n_mon[i]=n_mon[length-1-i]; n_mon[length-1-i]=xxx;
			xxx=mon_nr[i]; mon_nr[i]=mon_nr[length-1-i]; mon_nr[length-1-i]=xxx;
			xxx=d_mon[i]; d_mon[i]=d_mon[length-1-i]; d_mon[length-1-i]=xxx;
			for (int j=0; j<AlListLength; j++) {
				xxx=Al[j]->frag[i]; Al[j]->frag[i]=Al[j]->frag[length-1-i]; Al[j]->frag[length-1-i]=xxx;
			}
		}


//for (int i=0; i<n_generations; i++) cout<<"generation " << i << " first_a " << first_a[i] << " last_a " << last_a[i] << endl;
//length=n_arm.size();
//for (int i=0; i<length; i++) cout <<"arm " << i << " first_b " << first_b[i] << " last_b " << last_b[i] << endl;
//for (int i=0; i<length; i++) cout <<"arm " << i << " first_s " << first_s[i] << " last_s " << last_s[i] << endl;
//for (int i=0; i<length; i++) cout <<"arm " << i << " n_arm " << n_arm[i] << endl;
//length=n_mon.size();
//for (int i=0; i<length; i++) cout <<"block " << i << " n_mon " << n_mon[i] << " mon_nr " << mon_nr[i] << " d_mon " << d_mon[i] << endl;
	}
	success=MakeMonList();
	if (chainlength==1) MolType=monomer;


	return success;
}

int Molecule::GetChainlength(void){
if (debug) cout <<"GetChainlength for Mol " + name << endl;
	return chainlength;
}

bool Molecule:: MakeMonList(void) {
if (debug) cout <<"Molecule:: MakeMonList" << endl;
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
		//	cout << "in frag i " << i << " there is segment nr " <<  mon_nr[i] << " and it is on molmon  pos " << pos << endl;
		} else {cout <<"program error in mol PrepareForCalcualations" << endl; }
		i++;
	}
	return success;
}

bool Molecule::IsClamped() {
if (debug) cout <<"IsClamped for Mol " + name << endl;
	bool success=false;
	int length=mon_nr.size();
	if (mon_nr[0]==mon_nr[length-1]) {
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
	push("theta",theta);
	Real thetaexc=theta-phibulk*Lat[0]->Accesible_volume;
	push("theta_exc",thetaexc);
	push("n",n);
	push("chainlength",chainlength);
	push("phibulk",phibulk);
	push("mu",Mu);
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
	push("GN",GN);
	push("norm",norm);
	string s="profile;0"; push("phi",s);
	int length = MolMonList.size();
	for (int i=0; i<length; i++) {
		stringstream ss; ss<<i+1; string str=ss.str();
		s= "profile;"+str; push("phi-"+Seg[MolMonList[i]]->name,s);
	}
	for (int i=0; i<length_al; i++) {
		push(Al[i]->name+"value",Al[i]->value);
		push(Al[i]->name+"composition",Al[i]->composition);
		stringstream ss; ss<<i+length; string str=ss.str();
		s="profile;"+str; push(Al[i]->name+"-phi",s);
	}
	s="vector;0"; push("gn",s);
#ifdef CUDA
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
		if (sub[1]=="0") return H_phitot;
		int length=MolMonList.size();
		int i=0;
		while (i<length) {
			stringstream ss; ss<<i+1; string str=ss.str();
			if (sub[1]==str) return H_phi+i*M;
			i++;
		}
		int length_al=MolAlList.size();
		i=0;
		while (i<length_al) {
			stringstream ss; ss<<i+length; string str=ss.str();
			if (sub[i]==str) return Al[i]->H_phi;
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

Real Molecule::fraction(int segnr){
if (debug) cout <<"fraction for Mol " + name << endl;
	int Nseg=0;
	int length = mon_nr.size();
	int i=0;
	while (i<length) {
		if (segnr==mon_nr[i]) Nseg+=n_mon[i];
		i++;
	}
	return 1.0*Nseg/chainlength;
}

bool Molecule::ComputePhi(){
if (debug) cout <<"ComputePhi for Mol " + name << endl;
	bool success=true;
	Lat[0]->sub_box_on=0;//selecting 'standard' boundary condition
	switch (MolType) {
		case monomer:
			success=ComputePhiMon();
			break;
		case linear:
			if (freedom == "clamped") {
				Lat[0]->sub_box_on=Seg[mon_nr[0]]->clamp_nr; //slecting sub_box boundary conditions.
				//success=ComputePhiLin();
				success=ComputeClampLin();
				Lat[0]->sub_box_on=0;//selecting 'standard' boundary condition
			} else success=ComputePhiBra();
    	break;
		case branched:
			success=ComputePhiBra();
			break;
		case dendrimer:
			success=ComputePhiDendrimer();
			break;
		default:
			cout << "Programming error " << endl;
			break;
	}
	return success;
}

bool Molecule::ComputePhiMon(){
if (debug) cout <<"ComputePhiMon for Mol " + name << endl;
	int M=Lat[0]->M;
	bool success=true;
	Cp(phi,Seg[mon_nr[0]]->G1,M);
	//Cp(phi,G1,M);
	//Lat[0]->remove_bounds(phi);
	GN=Lat[0]->WeightedSum(phi);
	if (compute_phi_alias) {
		int length = MolAlList.size();
		for (int i=0; i<length; i++) {
			if (Al[i]->frag[0]==1) {
				Cp(Al[i]->phi,phi,M); Norm(Al[i]->phi,norm,M);
			}
		}
	}
	Times(phi,phi,Seg[mon_nr[0]]->G1,M);
	//Times(phi,phi,G1,M);
	return success;
}

Real* Molecule::propagate_forward(Real* G1, int &s, int block, int generation, int M) {
if (debug) cout <<"propagate_forward for Mol " + name << endl;

	int N= n_mon[block];
	if (save_memory) {
		int k,k0,t0,v0,t;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];
		if (s==first_s[generation]) { s++;
			Cp(Gs+M,G1,M); Cp(Gs,G1,M); Cp(Gg_f+n0*M,Gs+M,M); last_stored[block]=n0;
		} else {
			Lat[0] ->propagate(Gs,G1,0,1,M); //assuming Gs contains previous end-point distribution on pos zero;
			Cp(Gg_f+n0*M,Gs+M,M); s++;
			last_stored[block]=n0;
		}
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
			Cp(Gs,Gs+M,M); //make sure that Gs has last end-point distribution on spot zero.
			//return Gs;
		}
	} else {
		for (int k=0; k<N; k++) {
			if (s>first_s[generation]) {
				Lat[0] ->propagate(Gg_f,G1,s-1,s,M);
			} else {
				Cp(Gg_f+first_s[generation]*M,G1,M);
			}
			 s++;
		} 	}
	if (save_memory) {
		return Gg_f+last_stored[block]*M;
	} else return Gg_f+(s-1)*M;
}



void Molecule::propagate_backward(Real* G1, int &s, int block, int generation, int M) {
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
			AddTimes(rho+molmon_nr[block]*M,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,M);
			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) {
						Composition(Al[i]->rho,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,G1,norm,M);
					}
				}
			}
			s--;
		}
		Cp(Gg_b,Gg_b+M,M);  //make sure that on both spots the same end-point distribution is stored
	} else {
		for (int k=0; k<N; k++) {
			if (s<chainlength-1) Lat[0]->propagate(Gg_b,G1,(s+1)%2,s%2,M); else Cp(Gg_b+(s%2)*M,G1,M);

			AddTimes(rho+molmon_nr[block]*M,Gg_f+(s)*M,Gg_b+(s%2)*M,M);
			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) {
						Composition(Al[i]->rho,Gg_f+s*M,Gg_b+(s%2)*M,G1,norm,M);
					}
				}
			}
			s--;
		}
	}
}



bool Molecule::ComputeClampLin(){
	if (debug) cout <<"ComputeClampLin for Mol " + name << endl;
	bool success=true;
	int M=Lat[0]->M;
	int m=0;
	if (freedom=="clamped") m=Lat[0]->m[Seg[mon_nr[0]]->clamp_nr];
	int blocks=mon_nr.size();
	Zero(rho,m*n_box*MolMonList.size());
	int s=1;
	if (save_memory) {
		Cp(Gs,mask1,m*n_box);
	} else {
		Cp(Gg_f,mask1,m*n_box); //first block just contains the clamp
	}
	for (int i=1; i<blocks-1; i++) {
		Lat[0]->DistributeG1(Seg[mon_nr[i]]->G1,g1,Bx,By,Bz,n_box);
		//Lat[0]->DistributeG1(G1+molmon_nr[i]*M,g1,Bx,By,Bz,n_box);
		//propagate_forward(Gg_f,g1,s,n_mon[i],i,m*n_box);
		propagate_forward(g1,s,i,0,m*n_box);
	}
	if (save_memory) {
		int k=last_stored[blocks-2];
		int N=memory[n_mon.size()-1];
		Lat[0]->propagate(Gg_f,mask2,k,N-1,m*n_box);
		Lat[0]->ComputeGN(gn,Gg_f,H_Bx,H_By,H_Bz,H_Px2,H_Py2,H_Pz2,N-1,n_box);
	} else {
		Lat[0]->propagate(Gg_f,mask2,s-1,s,m*n_box);
		Lat[0]->ComputeGN(gn,Gg_f,H_Bx,H_By,H_Bz,H_Px2,H_Py2,H_Pz2,chainlength-1,n_box);
	}
	s=chainlength-1; //for last segment (clamp) no densities computed
	Cp(Gg_b+(s%2)*m*n_box,mask2,m*n_box); //last block just contains the clamp;
	if (save_memory) Cp(Gg_b+((s-1)%2)*m*n_box,Gg_b+(s%2)*m*n_box,m*n_box);
	s--;
	for (int i=blocks-2; i>0; i--) {
		Lat[0]->DistributeG1(Seg[mon_nr[i]]->G1,g1,Bx,By,Bz,n_box);
		//Lat[0]->DistributeG1(G1+molmon_nr[i]*M,g1,Bx,By,Bz,n_box);
		//propagate_backward(Gg_f,Gg_b,g1,s,n_mon[i],i,m*n_box);
		propagate_backward(g1,s,i,0,m*n_box);
	} //for first segment (clamp) no densities computed.
	int length=MolMonList.size();
	for (int i=1; i<length; i++) {
		Lat[0]->CollectPhi(phi+M*i,gn,rho+m*n_box*i,Bx,By,Bz,n_box);
	}

	return success;
}



void Molecule::BackwardBra(Real* G_start, int generation, int &s){//not yet robust for GPU computations: GS and GX need to be available on GPU
	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=Lat[0]->M;
	//Real* GS = new Real[4*M];
	Real* GS= (Real*) malloc(4*M*sizeof(Real));
	int k=bN;
	int ss=0;
	while (k>=b0){
		if (k>b0 && k<bN) {
			if (Gnr[k]!=generation) {
				Br.clear(); Gb.clear();
				while (Gnr[k] != generation){
					Br.push_back(Gnr[k]);
					if (save_memory) Gb.push_back(Gg_f+last_stored[k]*M); else Gb.push_back(Gg_f+last_s[Gnr[k]]*M);
					ss=first_s[Gnr[k]];
					k-=(last_b[Gnr[k]]-first_b[Gnr[k]]+1) ;
				}
				Br.push_back(generation); ss--;
				if (save_memory) Gb.push_back(Gg_f+last_stored[k]*M); else Gb.push_back(Gg_f+ss*M);
				int length = Br.size();
				//Real *GX = new Real[length*M];
				Real* GX= (Real*) malloc(length*M*sizeof(Real));
				for (int i=0; i<length; i++) Cp(GX+i*M,Gb[i],M);
				Cp(GS+3*M,Gg_b,M);
				for (int i=0; i<length; i++) {
					Cp(GS+2*M,GS+3*M,M);
					for (int j=0; j<length; j++) {
						if (i !=j) {
							Cp(GS,GX+j*M,M);
							Lat[0]->propagate(GS,UNITY,0,1,M);
							Times(GS+2*M,GS+2*M,GS+M,M);
						}
					}
					Cp(Gg_b,GS+2*M,M); Cp(Gg_b+M,GS+2*M,M);
					if (i<length-1) BackwardBra(Gg_b,Br[i],s);
				}
				free(GX);
				k++;
			} else propagate_backward(Seg[mon_nr[k]]->G1,s,k,generation,M);
				//propagate_backward(G1+molmon_nr[k]*M,s,k,generation,M);

		} else propagate_backward(Seg[mon_nr[k]]->G1,s,k,generation,M);
			//propagate_backward(G1+molmon_nr[k]*M,s,k,generation,M);

		k--;
	}
	free(GS);
}

Real* Molecule::ForwardBra(int generation, int &s) {
	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=Lat[0]->M;
	//Real* GS = new Real[3*M];
	Real* GS= (Real*) malloc(3*M*sizeof(Real));

	Real* Glast=NULL;
	int k=b0;
	while (k<=bN) {
		if (b0<k && k<bN) {
			if (Gnr[k]==generation ){
				Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,generation,M);
				//Glast=propagate_forward(G1+molmon_nr[k]*M,s,k,generation,M);
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
				//Lat[0]->propagate(GS,G1+molmon_nr[k]*M,0,2,M);

				for (int i=0; i<length; i++) {
					Cp(GS,Gb[i],M);
					Lat[0]->propagate(GS,UNITY,0,1,M);
					Times(GS+2*M,GS+2*M,GS+M,M);
				}
				if (save_memory) {
					Cp(Gs,GS+2*M,M); Cp(Gs+M,GS+2*M,M);
					Cp(Gg_f+(memory[k]-1)*M,GS+2*M,M); //correct because in this block there is just one segment.
				} else Cp(Gg_f+s*M,GS+2*M,M);
				s++;
			}
		} else {
			Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,generation,M);
			//Glast=propagate_forward(G1+molmon_nr[k]*M,s,k,generation,M);
		}
		k++;
	}
	free(GS);
	return Glast;
}

bool Molecule::ComputePhiBra() {
	if (debug) cout <<"ComputePhiBra for Mol " + name << endl;
	int M=Lat[0]->M;
	bool success=true;
	int generation=0;
	int s=0;
	Real* G=ForwardBra(generation,s);
	//Lat[0]->remove_bounds(G);
	GN=Lat[0]->WeightedSum(G);
//cout << "GN " << GN << endl;
	s--;
	if (save_memory) {Cp(Gg_b,Seg[mon_nr[last_b[0]]]->G1,M); Cp(Gg_b+M,Seg[mon_nr[last_b[0]]]->G1,M);} //toggle; initialize on both spots the same G1, so that we always get proper start.
	//if (save_memory) {Cp(Gg_b,G1+molmon_nr[last_b[0]]*M,M); Cp(Gg_b+M,G1+molmon_nr[last_b[0]]*M,M);} //toggle; initialize on both spots the same G1, so that we always get proper start.
	BackwardBra(Seg[mon_nr[last_b[0]]]->G1,generation,s);
	//BackwardBra(G1+molmon_nr[last_b[0]]*M,generation,s);
	return success;
}

Real* Molecule::propagate_forward(Real* G1, int &s, int block, int generation, int arm, int M) { //for dendrimer
if (debug) cout <<"propagate_forward for Molecule " + name << endl;
	int N= n_mon[block];
	if (save_memory) {
		int k,k0,t0,v0,t;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];
		if (s==first_s[arm] && generation ==0) { s++;
			Cp(Gs+M,G1,M); Cp(Gs,G1,M); Cp(Gg_f+n0*M,Gs+M,M); last_stored[block]=n0;
		} else {
			if (s==first_s[arm]) {
				//Nog niet af: Cp(Gs,last_stored[last_b[last_a[generation-1]]+1]*M,M); //Cp(Gs+M,Gs,M);
			}
			Lat[0] ->propagate(Gs,G1,0,1,M); //assuming Gs contains previous end-point distribution on pos zero;
			Cp(Gg_f+n0*M,Gs+M,M); s++;
			last_stored[block]=n0;
		}
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
			Cp(Gs,Gs+M,M); //make sure that Gs has last end-point distribution on spot zero.
			//return Gs;
		}
	} else {
		for (int k=0; k<N; k++) {

			if (s==first_s[arm] && generation==0) {
				Cp(Gg_f+first_s[generation]*M,G1,M);
			} else {
				if (s==first_s[arm])
					Lat[0] ->propagate(Gg_f,G1,last_s[last_a[generation-1]]+1,s,M);
				else Lat[0] ->propagate(Gg_f,G1,s-1,s,M);
			}
			 s++;
		}
	}
	if (save_memory) {
		return Gg_f+last_stored[block]*M;
	} else return Gg_f+(s-1)*M;
}


void Molecule::propagate_backward(Real* G1, int &s, int block, int generation, int arm, int M) { //this one is for the dendrimer molecules.
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
			AddTimes(rho+molmon_nr[block]*M,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,M);
			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) {
						Composition(Al[i]->rho,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,G1,norm,M);
					}
				}
			}
			s--;
		}
		Cp(Gg_b,Gg_b+M,M);  //make sure that on both spots the same end-point distribution is stored
	} else {
		for (int k=0; k<N; k++) {
			if (s<chainlength-1) Lat[0]->propagate(Gg_b,G1,(s+1)%2,s%2,M); else Cp(Gg_b+(s%2)*M,G1,M);

			AddTimes(rho+molmon_nr[block]*M,Gg_f+(s)*M,Gg_b+(s%2)*M,M);
			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) {
						Composition(Al[i]->rho,Gg_f+s*M,Gg_b+(s%2)*M,G1,norm,M);
					}
				}
			}
			s--;
		}
	}
}

void Molecule::BackwardDen(Real* GBr, int generation, int &s,int deg){//not yet robust for GPU computations: GS and GX need to be available on GPU
if (debug) cout << "BackwardDen in Molecule " << endl;
	int a0 = first_a[generation];
	int aN = last_a[generation];

	vector<int> Br;
	vector<Real*> Gb;
	int M=Lat[0]->M;
	Real* GS = new Real[3*M];
	int n_arms;
	propagate_backward(GBr+(generation-1)*M,s,last_b[aN]+1,generation,deg,M); //this should compute phi for first central segment in dendrimer s=N
	for (int a=aN; a>=a0; a--) {
		Cp(GS+2*M,GBr+(generation-1)*M,M);
		for (int aa=aN; aa>=a0; a--) {
			n_arms=n_arm[aa]; if (aa==a) n_arms--;
			Cp(Gs,Gg_f+last_stored[last_b[aa]],M);
			Lat[0]->propagate(GS,UNITY,0,1,M);
			for (int k=0; k<n_arms; k++) {Times(GS+2*M,GS+2*M,GS+M,M);}
		}
		//cp to GS
		int b0=first_b[a];
		int bN=last_b[a];
		for (int k=bN; k>=b0; k--) 	propagate_backward(Seg[mon_nr[k]]->G1,s,k,generation,a,M);
		//propagate_backward(G1+molmon_nr[k]*M,s,k,generation,a,M);
		//goto branch
		//bookkeep GBr
		//call BackwardDen for previous generation when generation still larger than 0
	}

	delete [] GS;
}


bool Molecule::ComputePhiDendrimer() {
	if (debug) cout <<"ComputePhiDendrimer for Mol " + name << endl;
	int M=Lat[0]->M;
	Real* GS = new Real[3*M];
	int s=0;
	Real* Glast=NULL;
	bool success=true;
	int n_g=first_a.size();
	for (int g=0; g<n_g; g++) {
		//int n_a=first_b.size();
		int a0=first_a[g],aN=last_a[g];
		for (int a=a0; a<=aN; a++) {
			Cp(GS+2*M,UNITY,M);
			int b0=first_b[a], bN=last_b[a];
			for (int b=b0; b<=bN; b++) {
				Glast=propagate_forward(Seg[mon_nr[b]]->G1,s,b,g,a,M);
				//Glast=propagate_forward(G1+molmon_nr[b]*M,s,b,g,a,M);
			}
			Cp(GS,Glast,M);
			Lat[0] ->propagate(GS,UNITY,0,1,M);
			for (int k=0; k<n_arm[a]; a++) Times(Gs+2*M,Gs+2*M,Gs+M,M);
		}
		Times(Gs+2*M,Gs+2*M,Seg[mon_nr[last_b[aN]+1]]->G1,M);
		//Times(Gs+2*M,Gs+2*M,G1+molmon_nr[last_b[aN]+1]*M,M);

		Glast=Gg_f+memory[last_b[aN]]; Cp(Glast,Gs+2*M,M);
		s++;
	}

	//Lat[0]->remove_bounds(Glast);
	GN=Lat[0]->WeightedSum(Glast);
	s=N;
	if (save_memory) {Cp(Gg_b,Seg[mon_nr[mon_nr.size()-1]]->G1,M); Cp(Gg_b+M,Seg[mon_nr[mon_nr.size()-1]]->G1,M);}
	//if (save_memory) {Cp(Gg_b,G1+molmon_nr[mon_nr.size()-1]*M,M); Cp(Gg_b+M,G1+molmon_nr[mon_nr.size()-1]*M,M);}
	Real* GBr=new Real[n_g*M];
	Cp(GBr+(n_g-1)*M,Seg[mon_nr[mon_nr.size()-1]]->G1,M);
	//Cp(GBr+(n_g-1)*M,G1+molmon_nr[mon_nr.size()-1]*M,M);
	BackwardDen(GBr,n_g,s,1);
	delete [] GBr;
	delete [] GS;
	return success;
}


/*   --------------------------------trash-----------------------------------------------------


bool Molecule::ComputePhiLin(){
	if (debug) cout <<"ComputePhiLin for Mol " + name << endl;
	int M=Lat[0]->M;
	bool success=true;
	int blocks=mon_nr.size();
	int s=0;
	Cp(Gg_f,Seg[mon_nr[0]]->G1,M);
	//Cp(Gg_f,G1+molmon_nr[0]]*M,M);
	//for (int i=0; i<blocks; i++) propagate_forward(Gg_f,Seg[mon_nr[i]]->G1,s,n_mon[i],i,M);
	for (int i=0; i<blocks; i++) propagate_forward(Gg_f,G1+molmon_nr[i]*M,s,n_mon[i],i,M);

	s=chainlength-1;
	Cp(Gg_b+(s%2)*M,Seg[mon_nr[blocks-1]]->G1,M);
	//Cp(Gg_b+(s%2)*M,G1+molmon_nr[blocks-1]*M,M);
	for (int i=blocks-1; i>-1; i--) propagate_backward(Gg_f,Gg_b,Seg[mon_nr[i]]->G1,s,n_mon[i],i,M);
	//for (int i=blocks-1; i>-1; i--) propagate_backward(Gg_f,Gg_b,G1+molmon_nr[i]*M,s,n_mon[i],i,M);
	//Lat[0]->remove_bounds(Gg_b);
	GN=Lat[0]->WeightedSum(Gg_b);
	return success;
}

void Molecule::propagate_backward(Real* Gg_f, Real* Gg_b, Real* G1, int &s, int N, int block, int M) {
if (debug) cout <<"propagate_backward for Mol " + name << endl;
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
			AddTimes(rho+molmon_nr[block]*M,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,M);
			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) { //alias is not yet read for clamp
						Composition(Al[i]->rho,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,G1,norm,M);
					}
				}
			}
			s--;
		}
		Cp(Gg_b,Gg_b+M,M);  //make sure that on both spots the same end-point distribution is stored
	} else {
		for (int k=0; k<N; k++) {
			if (s<chainlength-1) Lat[0]->propagate(Gg_b,G1,(s+1)%2,s%2,M);
			AddTimes(rho+molmon_nr[block]*M,Gg_f+(s)*M,Gg_b+(s%2)*M,M);
			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) { //alias is not yet ready for clamp
						Composition(Al[i]->rho,Gg_f+s*M,Gg_b+(s%2)*M,G1,norm,M);
					}
				}
			}
			s--;
		}
	}
}


void Molecule::propagate_forward(Real* Gg_f, Real* G1,int &s, int N, int block, int M) {
if (debug) cout <<"propagate_forward for Mol " + name << endl;
	//int M=Lat[0]->M;
	if (save_memory) {
		int k,k0,t0,v0,t;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];
		if (s==0) {
			Cp(Gs+M,G1,M);Cp(Gg_f,Gs+M,M);
		} else {
			Lat[0] ->propagate(Gs,G1,0,1,M); //assuming Gs contains previous end-point distribution on pos zero;
			Cp(Gg_f+n0*M,Gs+M,M); s++;
		}
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
		  	 || (2*(n-t+k) >= N-1))
				Cp(Gg_f+(n0+t-1)*M,Gs+(k%2)*M,M);

		}
		if ((N)%2!=0) {
			Cp(Gs,Gs+M,M); //make sure that Gs has last end-point distribution on spot zero.
		}
	} else
	for (int k=0; k<N; k++) {
		if (s>0) Lat[0] ->propagate(Gg_f,G1,s-1,s,M); s++;
	}
}



void ComputePhi(){
	Lat[0]->DisG1(G1,g1,Bx,By,Bz,n_box);
	Cp(Gg_f,mask1,M*n_box);
	for (s=1; s<=N; s++) Propagate(Gg_f,g1,s-1,s); Propagate(Gg_f,mask2,N,N+1);
	Cp(Gg_b+((N+1)%2)*M*n_box,mask2,M*n_box);
	Zero(rho,M*n_box);
	for (int s=N; s>0; s--) {
		Propagate(Gg_b,g1,((s+1)%2),(s%2));
		AddTimes(rho,Gg_f+s*M*n_box,Gg_b+(s%2)*M*n_box,M*n_box);
	}
	ComputeGN(GN,Gg_f,H_Bx,H_By,H_Bz,H_Px2,H_Py2,H_Pz2,N,M,n_box);
	Lat[0]->ColPhi(phi,GN,rho,Bx,By,Bz,n_box);
	Div(phi+MM,G1,MM);
	Add(phi+MM,MASK,MM);
}

void ComputeGN(Real* GN, Real* Gg_f, int* H_Bx, int* H_By, int* H_Bz, int* H_Px2, int* H_Py2, int* H_Pz2, int N, int M, int n_box) {
	for (int p=0; p<n_box; p++) Cp(GN+p,Gg_f+n_box*M*(N+1) +p*M+ jx*(H_Px2[p]-H_Bx[p])+jy*(H_Py2[p]-H_By[p])+(H_Pz2[p]-H_Bz[p]),1);

#ifdef CUDA //this transfer can go away when all is on GPU.
	TransferDataToHost(H_GN,GN,n_box);
#endif

}
*/

/*
bool Molecule::ExpandAlias(vector<string> sub, string &s) {
if (debug) cout <<"Molecule:: ExpandAlias" << endl;
	bool success=true;
	int correction=0;
	vector<int> open;
	vector<int> close;
	int length_al=sub.size();
	for (int i=0; i<length_al; i++) {
		open.clear(); close.clear();
		string sA;
		In[0]->EvenBrackets(sub[i],open,close);
		if (open.size() ==0) sA=sub[i]; else sA=sub[i].substr(0,open[0]);
		if (i==0 && sA=="") { //the 'composition' does not start with an alias. This should not be a problem.
		} else {

			if (sA[sA.length()-1]==',') {sA=sA.substr(0,sA.length()-1); correction=1;}
			if (sA[sA.length()-1]==';') {sA=sA.substr(0,sA.length()-1); correction=2;}

			if (!In[0]->InSet(In[0]->AliasList,sA)) {
				cout <<"In composition of mol '" + name + "' Alias '" + sA + "' was not found"<<endl; success=false;
			} else {
				int Alnr =GetAlNr(sA);
				if (Alnr<0) {
					Al.push_back(new Alias(In,Lat,sA));
					Alnr=Al.size();
					if (!Al[Alnr-1]->CheckInput(start)) {return false;}
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
					if (open[0]!=0) ssub=sub[i].substr(open[0]);
					switch(correction) {
						case 0:
							sub[i]=":"+si+":"+al_comp+":"+si+":"+ssub;
							break;
						case 1:
							sub[i]=":"+si+":"+al_comp+",:"+si+":"+ssub;
							break;
						case 2:
							sub[i]=":"+si+":"+al_comp+";:"+si+":"+ssub;
							break;
						default:
							break;

					}
				} else {
					string sss;
					stringstream sstm;
					sstm << iv;
					sss = sstm.str();
					sub[i]=sss+sub[i].substr(open[0]);
				}
			}
		}
	}
	string ss;
	for (int i=0; i<length_al; i++) {
		ss=ss.append(sub[i]);
	}
	s=ss;
	return success;
}
*/
