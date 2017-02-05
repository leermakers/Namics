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
if (debug) cout <<"Destructor for Mol " + name << endl;
	free(H_phi);
	free(H_phitot);
#ifdef CUDA
	cudaFree(phi);
	cudaFree(phitot);
	cudaFree(Gg_f);
	cudaFree(Gg_b);
	if (save_memory) cudaFree(Gs); 
#else
	free(Gg_f);
	free(Gg_b);
	if (save_memory) free(Gs); 
#endif
}

void Molecule:: AllocateMemory() {
if (debug) cout <<"AllocateMemory in Mol " + name << endl;
	int M=Lat[0]->M;
	if (save_memory) {
		int length_ = mon_nr.size();
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

	H_phi = (double*) malloc(M*MolMonList.size()*sizeof(double)); 
	H_phitot = (double*) malloc(M*sizeof(double)); 
#ifdef CUDA
	phi=(double*)AllOnDev(M*MolMonList.size());
	phitot=(double*)AllOnDev(M);
	Gg_f=(double*)AllOnDev(M*N); 
	Gg_b=(double*)AllOnDev(M*2);
	if (save_memory) {
		Gs=(double*)AllOnDev(M*2);
	}	
#else
	phi = H_phi;
	phitot = H_phitot;
	Gg_f = (double*) malloc(M*N*sizeof(double)); 
	Gg_b = (double*) malloc(M*2*sizeof(double)); 
	if (save_memory) {
		Gs = (double*) malloc(M*2*sizeof(double));
	}
#endif
	Zero(Gg_f,M*N);
	Zero(Gg_b,2*M);
	if (save_memory) Zero(Gs,2*M);
	int length =MolAlList.size();
	for (int i=0; i<length; i++) Al[i]->AllocateMemory(); 
}

bool Molecule:: PrepareForCalculations() {
if (debug) cout <<"PrepareForCalculations in Mol " + name << endl; 
	int M=Lat[0]->M;
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
		if (GetValue("save_memory").size()>0) {
			save_memory=In[0]->Get_bool(GetValue("save_memory"),false); 
		}
		if (GetValue("composition").size()==0) {cout << "For mol '" + name + "' the definition of 'composition' is required" << endl; success = false;
		} else {
			if (!Decomposition(GetValue("composition"))) {cout << "For mol '" + name + "' the composition is rejected. " << endl; success=false;}

		}
		if (GetValue("freedom").size()==0 && !IsTagged()) {
			cout <<"For mol " + name + " the setting 'freedom' is expected. Problem terminated " << endl; success = false;
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
				}
					if (IsCharged()) {success=false; cout << "Mol '" + name + "' is not 'charged' and therefore this molecule can not be the neutralizer" << endl; } 
				if (freedom == "free") {
					if (GetValue("theta").size()>0 || GetValue("n").size() > 0) {
						cout << "In mol " + name + ", the setting of 'freedom = free', can not not be combined with 'theta' or 'n': use 'phibulk' instead." << endl; success=false;  
					} else {
						if (GetValue("phibulk").size() ==0) {
							cout <<"In mol " + name + ", the setting 'freedom = free' should be combined with a value for 'phibulk'. "<<endl; success=false;
						} else {
							phibulk=In[0]->Get_double(GetValue("phibulk"),-1); 
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

								if (GetValue("n").size()>0) {n=In[0]->Get_double(GetValue("n"),10*Lat[0]->volume);theta=n*chainlength;}
								if (GetValue("theta").size()>0) {theta = In[0]->Get_double(GetValue("theta"),10*Lat[0]->volume);n=theta/chainlength;} 
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

bool Molecule::Decomposition(string s){
if (debug) cout <<"Decomposition for Mol " + name << endl;
	bool success = true;
	vector<int> open;
	vector<int> close;
	vector<string>sub;
	In[0]->split(s,'@',sub);
	aliases=(s!=sub[0]);
	if (aliases) {//first do the work for aliases....
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
	}
	bool done=false; //now interpreted the (expended) composition
	while (!done) { done = true; 
		open.clear(); close.clear();
		if (!In[0]->EvenBrackets(s,open,close)) {cout << "In composition of mol '" + name + "' the backets are not balanced."<<endl; success=false; }
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

	sub.clear();
	In[0]->split(s,':',sub);
	int length_sub =sub.size();
	int AlListLength=MolAlList.size();
	for (int i=0;i<AlListLength; i++) Al[i]->active=false; 
	int i=0; chainlength=0; MolMonList.clear();
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
					mon_nr.push_back(mnr); 
					for (int i=0; i<AlListLength; i++) {if (Al[i]->active) Al[i]->frag.push_back(1); else Al[i]->frag.push_back(0);}
				}
				int nn = In[0]->Get_int(sub[i].substr(close[k]+1,s.size()-close[k]),0); 
				if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity " << endl; success=false;
				} else {n_mon.push_back(nn); }
				chainlength +=nn;
				k++;
			}
			
		}
		
		i++;
	}
	success=MakeMonList();
	if (chainlength==1) MolType=monomer; else MolType=linear;	
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

bool Molecule::IsCharged() {
if (debug) cout <<"IsCharged for Mol " + name << endl;
	double charge =0;
	int length = n_mon.size(); 
	int i=0;
	while (i<length) {
		charge +=n_mon[i]*Seg[mon_nr[i]]->valence; 
		i++;
	}
	return charge<-1e-5 || charge > 1e-5; 
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

void Molecule::push(string s, double X) {
if (debug) cout <<"push (double) for Mol " + name << endl;
	doubles.push_back(s);
	doubles_value.push_back(X); 
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
	doubles.clear();
	doubles_value.clear();
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

double* Molecule::GetPointer(string s) {
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

int Molecule::GetValue(string prop,int &int_result,double &double_result,string &string_result){
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
	length = doubles.size();
	while (i<length) {
		if (prop==doubles[i]) { 
			double_result=doubles_value[i];
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

double Molecule::fraction(int segnr){
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
	switch (MolType) {
		case monomer:
			success=ComputePhiMon();
		break;
		case linear:			
			success=ComputePhiLin();
		break;
		default:
		cout << "Programming error " << endl; 
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

void Molecule::propagate_forward(double* Gg_f,double* G1,int &s, int N, int block) {
if (debug) cout <<"propagate_forward for Mol " + name << endl;
	int M=Lat[0]->M;
	if (save_memory) {
		int k,k0,t0,v0,t;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];
		if (s==0) {
			Cp(Gs+M,G1,M);
		} else {
			Lat[0] ->propagate(Gs,G1,0,1); //assuming Gs contains previous end-point distribution on pos 1; 
			Cp(Gg_f+n0*M,Gs+M,M);
		} 
		t=1;
		v0=t0=k0=0;
		for (k=2; k<=N; k++) {
			t++; s++;
			Lat[0]->propagate(Gs,G1,(k-1)%2,k%2);
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
			Cp(Gs,Gs+M,M); //make sure that Gs has last end-point distribution on both spots.
		} 
	} else
	for (int k=0; k<N; k++) {if (s>0) Lat[0] ->propagate(Gg_f,G1,s-1,s); s++;} 	
}

void Molecule::propagate_backward(double* Gg_f, double* Gg_b, double* G1, int &s, int N, int block) {
if (debug) cout <<"propagate_backward for Mol " + name << endl;
	int M = Lat[0]->M; 
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
					Lat[0]->propagate(Gg_b,G1,(k+1)%2,k%2);
				}
			} else {
				Lat[0]->propagate(Gg_b,G1,(k+1)%2,k%2);
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
					Lat[0]->propagate(Gs,G1,(t-1)%2,t%2);
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
			AddTimes(phi+molmon_nr[block]*M,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,M); 
			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) {
						Composition(Al[i]->phi,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,G1,norm,M);
					}
				}
			}
			s--;
		}
		Cp(Gg_b,Gg_b+M,M);  //make sure that on both spots the same end-point distribution is stored
	} else {
		for (int k=0; k<N; k++) {
			if (s<chainlength-1) Lat[0]->propagate(Gg_b,G1,(s+1)%2,s%2);
			AddTimes(phi+molmon_nr[block]*M,Gg_f+(s)*M,Gg_b+(s%2)*M,M); 
			if (compute_phi_alias) {
				int length = MolAlList.size();
				for (int i=0; i<length; i++) {
					if (Al[i]->frag[k]==1) {
						Composition(Al[i]->phi,Gg_f+s*M,Gg_b+(s%2)*M,G1,norm,M);
					}
				}
			}
			s--;
		} 
	}	
}

bool Molecule::ComputePhiLin(){
	if (debug) cout <<"ComputePhiLin for Mol " + name << endl;
	int M=Lat[0]->M;
	bool success=true;
	int blocks=mon_nr.size(); 
	int s=0; Cp(Gg_f,Seg[mon_nr[0]]->G1,M); 
	for (int i=0; i<blocks; i++) {
		propagate_forward(Gg_f,Seg[mon_nr[i]]->G1,s,n_mon[i],i);
	} 
	s=chainlength-1; Cp(Gg_b+(s%2)*M,Seg[mon_nr[blocks-1]]->G1,M);
	for (int i=blocks-1; i>-1; i--) propagate_backward(Gg_f,Gg_b,Seg[mon_nr[i]]->G1,s,n_mon[i],i);
	Lat[0]->remove_bounds(Gg_b);
	GN=Lat[0]->WeightedSum(Gg_b);
	return success;
}

