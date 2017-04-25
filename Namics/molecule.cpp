#include "molecule.h"
#include <sstream> 

Molecule::Molecule(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, string name_) {
	In=In_; Seg=Seg_; name=name_;  Lat=Lat_; 
if (debug) cout <<"Constructor for Mol " + name << endl;
	KEYS.push_back("freedom"); 
	KEYS.push_back("composition"); 
	KEYS.push_back("theta");
	KEYS.push_back("phibulk");
	KEYS.push_back("n");
	KEYS.push_back("save_memory"); 
}

Molecule::~Molecule() {
	DeAllocateMemory();
}

void Molecule :: DeAllocateMemory(){
if (debug) cout <<"Destructor for Mol " + name << endl;
	free(H_phi);
	free(H_phitot);
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
	free(UNITY);
	free(Gg_f);
	free(Gg_b);
	if (save_memory) free(Gs);
#endif
}

void Molecule:: AllocateMemory() {
if (debug) cout <<"AllocateMemory in Mol " + name << endl;
	int M=Lat[0]->M;
	int m=0;
	if (freedom=="clamped") m=Lat[0]->m[Seg[mon_nr[0]]->clamp_nr];

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
	int N; 
	if (save_memory) {
		N=memory[n_mon.size()-1]; 
	} else N=chainlength;

	H_phi = (Real*) malloc(M*MolMonList.size()*sizeof(Real)); 
	H_phitot = (Real*) malloc(M*sizeof(Real)); 
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
		Gg_f=(Real*)AllOnDev(m*n_box*N);
		Gg_b=(Real*)AllOnDev(m*n_box*2);
		g1=(Real*)AllOnDev(m*n_box);
		rho=(Real*)AllOnDev(m*n_box*MolMonList.size());
		phi=(Real*)AllOnDev(M*MolMonList.size());
		if (save_memory) Gs=(Real*)AllOnDev(m*n_box*2);
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
		Px1=H_Px1; Py1=H_Py1; Pz1=H_Pz2;
		Px2=H_Px2; Py2=H_Py2; Pz2=H_Pz2;
		Gg_f = (Real*) malloc(m*N*n_box*sizeof(Real));
		Gg_b = (Real*) malloc(m*2*n_box*sizeof(Real));
		g1=(Real*) malloc(m*n_box*sizeof(Real));
		rho=(Real*)malloc(m*n_box*MolMonList.size()*sizeof(Real));
		if (save_memory) Gs=(Real*) malloc(m*n_box*2*sizeof(Real));
		phi=H_phi;
	} else {
		Gg_f = (Real*) malloc(M*N*sizeof(Real));
		Gg_b = (Real*) malloc(M*2*sizeof(Real));
		phi=H_phi;
		rho=phi;
		if (save_memory) Gs=(Real*) malloc(2*M*sizeof(Real));
	}		
	phitot = H_phitot;
	UNITY = (Real*) malloc(M*sizeof(Real));
#endif
	Zero(Gg_f,M*N);
	Zero(Gg_b,2*M);
	if (save_memory) Zero(Gs,2*M);
	int length =MolAlList.size();
	for (int i=0; i<length; i++) Al[i]->AllocateMemory(Seg[mon_nr[0]]->clamp_nr,n_box); 
}

bool Molecule:: PrepareForCalculations(int *KSAM) {
if (debug) cout <<"PrepareForCalculations in Mol " + name << endl; 
	int M=Lat[0]->M;
	if (freedom=="clamped") {
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
#endif
	}
	Cp(UNITY,KSAM,M);
	int length_al=MolAlList.size();
	for (int i=0; i<length_al; i++) Al[i]->PrepareForCalculations();
	bool success=true;
	Zero(phitot,M);
	Zero(phi,M*MolMonList.size()); 
	compute_phi_alias=false;
	return success;
}

bool Molecule::CheckInput(int start_) {
start=start_;
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
			cout <<"For mol " + name + " the setting 'freedom' is expected: options: 'free' 'restricted' 'solvent' 'neutralizer' . Problem terminated " << endl; success = false;
			} else { if (!IsTagged()) {
				vector<string> free_list; 
				if (!IsPinned()) free_list.push_back("free"); free_list.push_back("restricted"); free_list.push_back("solvent"); 
				free_list.push_back("neutralizer");
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
				if (freedom == "restricted") {
					if (GetValue("phibulk").size()>0) {
						cout << "In mol " + name + ", the setting of 'freedom = restricted', can not not be combined with 'phibulk'  use 'theta' or 'n'  instead." << endl; success=false;  
					} else {
						if (GetValue("theta").size() ==0 && GetValue("n").size()==0) {
							cout <<"In mol " + name + ", the setting 'freedom = restricted' should be combined with a value for 'theta' or 'n'; do not use both settings! "<<endl; success=false;
						} else {
							if (GetValue("theta").size() >0 && GetValue("n").size()>0) {
							cout <<"In mol " + name + ", the setting 'freedom = restricted' do not specify both 'n' and 'theta' "<<endl; success=false;
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
 
			} else {
				if (GetValue("theta").size() >0 || GetValue("n").size() > 0 || GetValue("phibulk").size() >0 || GetValue("freedom").size() > 0) cout <<"Warning. In mol " + name + " tagged segment(s) were detected. In this case no value for 'freedom' is needed, and also 'theta', 'n' and 'phibulk' values are ignored. " << endl;  
			}
		} 	 	
	}	
	return success; 
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
if (debug) cout <<"GetMonNr for Mol " + name << endl;
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
	bool success=true;
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
					sub[i]=":"+si+":"+al_comp+":"+si+":"+ssub;
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

bool Molecule::ExpandBrackets(string &s) {
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
		int i_open=0; pos_open=open[0]; 
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
					sC=s.substr(pos_open,s.size()-pos_open);  
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
				int nn = In[0]->Get_int(sub[i].substr(close[k]+1,s.size()-close[k]),0); 
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
//cout <<"string ss IN ][ " << ss << "and new_generation " << new_generation << endl; 
				pos=pos+1; 
				first_s.push_back(-1);
				last_s.push_back(-1);
				first_b.push_back(-1);
				last_b.push_back(-1);
				GenerateTree(s,new_generation,pos,open,close);
				pos_close=pos_open+1;
			} else {
				pos=pos_close;
//cout <<"go to interpret ] " << ss << " and generation " << generation << " pos " << pos <<  endl;
				Interpret(ss,generation);
			}
		} else {
			ss=s.substr(pos,pos_open-pos);
			pos=pos_open;
//cout <<"go to interpret [ " << ss << " and generation " << generation << " pos " << pos << " new generation " << newgeneration << endl; 
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


bool Molecule::Decomposition(string s){
if (debug) cout <<"Decomposition for Mol " + name << endl;
	bool success = true;
	bool aliases = true;
	MolType=linear;//default; 
	int loopnr=0;
	vector<int> open;
	vector<int> close;
	vector<string>sub;
	while (aliases && success) {
		loopnr++;
		sub.clear();
		In[0]->split(s,'#',sub);
		aliases=(s!=sub[0]);
		if (aliases) ExpandAlias(sub,s);
		if (loopnr == 20) {
			cout << "Nesting nr 20 reached in aliases for mol " + name + " -composition. It is decided that this is too deep to continue; Possible, you have defined an alias-A inside alias-B which itself refers to alias-A. This is not allowed. Problem terminated. " << endl;  
			success=false; 
		}
	}
	if (!ExpandBrackets(s)) success=false;

	//test dendrimer
	//test asymmetric dendrimer
	//test ring
	//test comb
	//test star

	if (!In[0]->EvenSquareBrackets(s,open,close)) {cout << "Error in composition of mol '" + name + "'; the square brackets are not balanced. " << endl; success=false;  }
	if (open.size()>0) MolType=branched;

	int generation=0;
	int pos=0;
	chainlength=0;
	MolMonList.clear();
	int AlListLength=MolAlList.size();
	for (int i=0;i<AlListLength; i++) Al[i]->active=false; 
	first_s.push_back(-1);
	last_s.push_back(-1);
	first_b.push_back(-1);
	last_b.push_back(-1);
	GenerateTree(s,generation,pos,open,close); 

	//invert nubers;
	if (MolType==branched) {
		int g_length=first_s.size();
		int length = n_mon.size(); 
		int ss;
		chainlength=last_s[0]; 	

		for (int i=0; i<g_length; i++) { 
			first_s[i] = chainlength-first_s[i]-1;
			last_s[i] = chainlength-last_s[i]-1;
			ss=first_s[i]; first_s[i]=last_s[i]+1; last_s[i]=ss;
			first_b[i]=length-first_b[i]-1;
			last_b[i]=length-last_b[i]-1;
			ss=first_b[i]; first_b[i]=last_b[i]; last_b[i]=ss; 
		}
		
		int xxx;
		for (int i=0; i<length/2; i++) {
			xxx=Gnr[i]; Gnr[i]=Gnr[length-1-i]; Gnr[length-1-i]=xxx;
			xxx=n_mon[i]; n_mon[i]=n_mon[length-1-i]; n_mon[length-1-i]=xxx;
			xxx=mon_nr[i]; mon_nr[i]=mon_nr[length-1-i]; mon_nr[length-1-i]=xxx; 
			for (int j=0; j<AlListLength; j++) {
				xxx=Al[j]->frag[i]; Al[j]->frag[i]=Al[j]->frag[length-1-i]; Al[j]->frag[length-1-i]=xxx;	
			} 
		}
	}

	success=MakeMonList();
	if (chainlength==1) MolType=monomer;	

//	if (MolType==branched) {
//		int length=Gnr.size();
//		for (int i=0; i<length; i++) 
//			cout << "segnr " << mon_nr[i] << " nn " << n_mon[i] << " Gnr " << Gnr[i] << endl;  
//		length=last_s.size();
//		for (int i=0; i<length; i++) 
//			cout << "Gen " << i  << " first_s " << first_s[i] << " last_s " << last_s[i] << " first_b " << first_b[i] << " last_b " <<  last_b[i] << endl; 
//	}
	return success; 
}

int Molecule::GetChainlength(void){
if (debug) cout <<"GetChainlength for Mol " + name << endl;
	return chainlength; 
} 

bool Molecule:: MakeMonList(void) {
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
	Real charge=0;
	int length=mon_nr.size();
	for (int i=0; i<length; i++) {
		charge +=Seg[mon_nr[i]]->valence*n_mon[i];
	} 
	return charge;
}

bool Molecule::IsCharged() {
if (debug) cout <<"IsCharged for Mol " + name << endl;
	Real charge =0;
	int length = n_mon.size(); 
	int i=0;
	while (i<length) {
		charge +=n_mon[i]*Seg[mon_nr[i]]->valence; 
		i++;
	}
	return charge<-1e-6 || charge > 1e-6; 
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
	push("theta exc",theta-phibulk*Lat[0]->volume);
	push("n",n);
	push("chainlength",chainlength);
	push("phibulk",phibulk);
	push("Mu",Mu);
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
#ifdef CUDA
	TransferDataToHost(H_phitot,phitot,M);
	TransferDataToHost(H_phi,phi,M*MolMonList.size());
#endif	
}

Real* Molecule::GetPointer(string s) {
if (debug) cout <<"GetPointer for Mol " + name << endl;	vector<string> sub;
	int M= Lat[0]->M;
	In[0]->split(s,';',sub);
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
	Lat[0]->remove_bounds(phi);
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
	return success;
}

Real* Molecule::propagate_forward(Real* G1, int &s, int block, int generation, int M) {
if (debug) cout <<"propagate_forward for Mol " + name << endl;
	
	// Real* G1 = Seg[mon_nr[block]]->G1;
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
		} 
	}
	if (save_memory) {
		return Gg_f+last_stored[block]*M;
	} else return Gg_f+(s-1)*M;	
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

void Molecule::propagate_backward(Real* G1, int &s, int block, int generation, int M) {
if (debug) cout <<"propagate_backward for Mol " + name << endl;
	//Real *G1=Seg[mon_nr[block]]->G1;
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

bool Molecule::ComputeClampLin(){
	if (debug) cout <<"ComputeClampLin for Mol " + name << endl;
	bool success=true;
	int M=Lat[0]->M;
	int m=Lat[0]->m[Seg[mon_nr[0]]->clamp_nr];	
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
		//propagate_backward(Gg_f,Gg_b,g1,s,n_mon[i],i,m*n_box);	
		propagate_backward(g1,s,i,0,m*n_box);
	} //for first segment (clamp) no densities computed. 
	int length=MolMonList.size();
	for (int i=1; i<length; i++) {
		Lat[0]->CollectPhi(phi+M*i,gn,rho+m*n_box*i,Bx,By,Bz,n_box);
	}

	return success;
}

bool Molecule::ComputePhiLin(){
	if (debug) cout <<"ComputePhiLin for Mol " + name << endl;
	int M=Lat[0]->M;
	bool success=true;
	int blocks=mon_nr.size(); 
	int s=0; Cp(Gg_f,Seg[mon_nr[0]]->G1,M); 
	for (int i=0; i<blocks; i++) propagate_forward(Gg_f,Seg[mon_nr[i]]->G1,s,n_mon[i],i,M);
	s=chainlength-1; Cp(Gg_b+(s%2)*M,Seg[mon_nr[blocks-1]]->G1,M);
	for (int i=blocks-1; i>-1; i--) propagate_backward(Gg_f,Gg_b,Seg[mon_nr[i]]->G1,s,n_mon[i],i,M);
	Lat[0]->remove_bounds(Gg_b);
	GN=Lat[0]->WeightedSum(Gg_b);
	return success;
}

void Molecule::Backward(Real* G_start, int generation, int &s){//not yet robust for GPU computations: GS and GX need to be available on GPU 
	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;	
	int M=Lat[0]->M;
	Real* GS = new Real[4*M];
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
				Real *GX = new Real[length*M];
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
					if (i<length-1) Backward(Gg_b,Br[i],s);
				}
				delete GX;
				k++;
			} else 	propagate_backward(Seg[mon_nr[k]]->G1,s,k,generation,M);
		} else propagate_backward(Seg[mon_nr[k]]->G1,s,k,generation,M);
		k--; 
	}
	delete GS;
}

Real* Molecule::Forward(int generation, int &s) { 
	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=Lat[0]->M;
	Real* GS = new Real[3*M]; 
	Real* Glast=NULL;   
	int k=b0; 
	while (k<=bN) {
		if (b0<k && k<bN) { 
			if (Gnr[k]==generation ){
				Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,generation,M);	 
			} else {
				Br.clear(); Gb.clear();
				Cp(GS,Glast,M);
				while (Gnr[k] !=generation) {
					Br.push_back(Gnr[k]);
					Gb.push_back(Forward(Gnr[k],s));
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
					Cp(Gg_f+(memory[k]-1)*M,GS+2*M,M); //correct becuse in this block there is just one segment. 
				} else Cp(Gg_f+s*M,GS+2*M,M);			
				s++; 
			}
		} else {
			Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,generation,M);
		}
		k++;
	}
	delete GS; 
	return Glast;
}


bool Molecule::ComputePhiBra() {
	int M=Lat[0]->M;
	if (debug) cout <<"ComputePhiBra for Mol " + name << endl; 
	bool success=true;
	int generation=0;
	int s=0; 
	Real* G=Forward(generation,s);
	Lat[0]->remove_bounds(G);
	GN=Lat[0]->WeightedSum(G);
	s--;
	if (save_memory) {Cp(Gg_b,Seg[mon_nr[last_b[0]]]->G1,M); Cp(Gg_b+M,Seg[mon_nr[last_b[0]]]->G1,M);}
	Backward(Seg[mon_nr[last_b[0]]]->G1,generation,s);	
	return success;
}

/*   --------------------------------trash-----------------------------------------------------
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

