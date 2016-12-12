
class Molecule { 
public:
	Molecule(vector<Input*>,vector<Segment*>,string,int,int);

~Molecule();

	string name; 
	vector<Input*> In;
	vector<Segment*> Seg; 
	vector<string> MolMonList;
	int n_mol; 
	int M; 
	int mol_nr; 
	double theta;
	double phibulk; 
	string freedom;
	double n; 
	double GN; 
	int chainlength;
	bool save_memory; 
	string composition;
	vector<int> mon_nr;
	vector<int> n_mon; 
	float *phi;
	float *Gg_f;
	float *Gg_b; 
	int tag_segment; 
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(void);
	void PutParameter(string);
	bool Decomposition(string); 
	int GetChainlength(void); 
	double Theta(void);
	string GetValue(string); 
	int GetMonNr(string);
	bool MakeMonList(void);  
	bool IsPinned(void);
	bool IsTagged(void); 
	bool IsCharged(void); 
	void PrepareForCalculations(void);
	string GetFreedom(void); 


};
Molecule::Molecule(vector<Input*> In_,vector<Segment*> Seg_, string name_,int molnr,int N_mol) {
	In=In_; Seg=Seg_; name=name_; n_mol=N_mol; mol_nr=molnr; 
	KEYS.push_back("freedom"); 
	KEYS.push_back("composition"); //KEYS.push_back("chainlength");
	KEYS.push_back("theta");
	KEYS.push_back("phibulk");
	KEYS.push_back("n");
	KEYS.push_back("save_memory"); 
}
Molecule::~Molecule() {
}

bool Molecule::CheckInput() {
	bool success=true;
	if (!In[0]->CheckParameters("mol",name,KEYS,PARAMETERS,VALUES)) {
		success=false; 	
	} else { 
		if (GetValue("save_memory").size()>0) {
			save_memory=In[0]->Get_bool(GetValue("save_memory"),false); 
			if  (save_memory) { save_memory =false; cout << "For mol '" + name + "' the flag 'save_memory' is not yet activated. Input ignored " << endl;}
		}
		if (GetValue("composition").size()==0) {cout << "For mol '" + name + "' the definition of 'composition' is required" << endl; success = false;
		} else {
			if (!Decomposition(GetValue("composition"))) {cout << "For mol '" + name + "' the composition is rejected. " << endl; success=false;}

		}
		if (GetValue("freedom").size()==0 && !IsTagged()) {
			cout <<"For mol " + name + " the setting 'freedom' is expected. Problem terminated " << endl; success = false;
			} else { if (!IsTagged()) {
				string freedom_value;
				vector<string> free_list; 
				if (!IsPinned()) free_list.push_back("free"); free_list.push_back("restricted"); free_list.push_back("solvent"); 
				free_list.push_back("neutralizer");
				cout <<  GetValue("freedom") << endl; 
				if (!In[0]->Get_string(GetValue("freedom"),freedom_value,free_list,"In mol " + name + " the value for 'freedom' is not recognised ")) success=false;
				if (freedom_value == "solvent") {
					if (IsPinned()) {success=false; cout << "Mol '" + name + "' is 'pinned' and therefore this molecule can not be the solvent" << endl; } 
				}
				if (freedom_value == "neutralizer") {
					if (IsPinned()) {success=false; cout << "Mol '" + name + "' is 'pinned' and therefore this molecule can not be the neutralizer" << endl; } 
				}
					if (IsCharged()) {success=false; cout << "Mol '" + name + "' is not 'charged' and therefore this molecule can not be the neutralizer" << endl; } 
				if (freedom_value == "free") {
					if (GetValue("theta").size()>0 || GetValue("n").size() > 0) {
						cout << "In mol " + name + ", the setting of 'freedom = free', can not not be combined with 'theta' or 'n': use 'phibulk' instead." << endl; success=false;  
					} else {
						if (GetValue("phibulk").size() ==0) {
							cout <<"In mol " + name + ", the setting 'freedom = free' should be combined with a value for 'phibulk'. "<<endl; success=false;
						} else {
							phibulk=In[0]->Get_float(GetValue("phibulk"),-1); 
							if (phibulk < 0 || phibulk >1) {
								cout << "In mol " + name + ", the value of 'phibulk' is out of range 0 .. 1." << endl; success=false;
							}
						} 
					}
				}
				if (freedom_value == "restricted") {
					if (GetValue("phibulk").size()>0) {
						cout << "In mol " + name + ", the setting of 'freedom = restricted', can not not be combined with 'phibulk'  use 'theta' or 'n'  instead." << endl; success=false;  
					} else {
						if (GetValue("theta").size() ==0 && GetValue("n").size()==0) {
							cout <<"In mol " + name + ", the setting 'freedom = restriced' should be combined with a value for 'theta' or 'n'; do not use both settings! "<<endl; success=false;
						} else {
							float volume =  Seg[0]->MX*Seg[0]->MY*Seg[0]->MZ;
							if (GetValue("n").size()>0) n=In[0]->Get_float(GetValue("n"),10*volume);
							if (GetValue("theta").size()>0) theta = In[0]->Get_float(GetValue("theta"),10*volume);
							if (theta ==10*volume ) theta=n*chainlength; 
							if (chainlength >0) {if (n==10*volume) n=theta/chainlength;}
							if (theta < 0 || theta > volume) {
								cout << "In mol " + name + ", the value of 'n' or 'theta' is out of range 0 .. 'volume', cq 'volume'/N." << endl; success=false;
							
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

int Molecule:: GetMonNr(string s){
	int n_segments=In[0]->MonList.size();
	int found=-1;
	int i=0;
	while(i<n_segments) {
		if (Seg[i]->name ==s) found=i; 
		i++;
	}
	return found; 
}

bool Molecule:: Decomposition(string s){
	bool success = true;
	vector<int> open;
	vector<int> close;
	bool done=false;
	while (!done) { done = true; 
		open.clear(); close.clear();
		if (!In[0]->EvenBrackets(s,open,close)) {cout << "In composition of mol '" + name + "' the backets are not balanced."<<endl; success=false; 
		} else {
			int length=open.size();
			int i=0;
			int pos_low=-1;
			int pos_high=-1;
			int Pos_high=-1;
			while (i<length-1) {
				if (pos_low>0 && pos_high<0) {if (open[i+1]>close[i]) {pos_high=close[i]; Pos_high=open[i+1];}}
				if (pos_low<0) {if (open[i+1]<close[i]) pos_low=open[i]; }
				i++;
			} 
			if (pos_low >0 && pos_high > 0) { done=false; 
				string sA,sB,sC;
				sA=s.substr(0,pos_low);
				sB=s.substr(pos_low+1,pos_high-pos_low-1);
				sC=s.substr(Pos_high,length-Pos_high);
				int x = In[0]->Get_int(s.substr(pos_high+1,length-pos_high),0);
				s=sA;
				for (int k=0; k<x; k++) s.append(sB); 
				s=s.append(sC);
			} 
		}
	}
	int length = open.size();
	int k=0; chainlength=0;
	while (k<length) {
		string segname=s.substr(open[k]+1,close[k]-open[k]-1); 
		int mnr=GetMonNr(segname); 
		if (mnr <0)  {cout <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised"  << endl; success=false; 
		} else mon_nr.push_back(mnr);  
		int nn = In[0]->Get_int(s.substr(close[k]+1,s.size()-close[k]),0); 
		if (nn<1) {cout <<"In composiiton of mol '" + name + "' the number of repeats should have values larger than unity " << endl; success=false;
		} else {n_mon.push_back(nn); }
		chainlength +=nn;
		k++;
	}
	success = MakeMonList(); 
		
	return success; 
}

int Molecule::GetChainlength(void){
	return chainlength; 
} 

bool Molecule:: MakeMonList(void) {
	bool success=true;
	int length = mon_nr.size();
	int i=0;
	while (i<length) {
		string Seg_name=Seg[mon_nr[i]]->name;
		if (!In[0]->InSet(MolMonList,Seg_name)) {
			if (Seg[mon_nr[i]]->GetFreedom()=="frozen") {
				success = false;
				cout << "In 'composition of mol " + name + ", a segment was found with freedom 'frozen'. This is not permitted. " << endl; 

			}
			MolMonList.push_back(Seg_name);
		}
		i++;
	}
	return success;
}
	
string Molecule::GetFreedom(void){
	return freedom; 
} 

bool Molecule::IsPinned() {
	bool success=false;
	int seg_nr=0; 
	int NumSegs=In[0]->MonList.size(); 
	int length=MolMonList.size();
	int i=0;
	while (i<length) {
		string segmentname=MolMonList[i]; 
		for (int j=0; j<NumSegs; j++) if (Seg[j]->name==segmentname) seg_nr=j;
		if (Seg[seg_nr]->GetFreedom()=="pinned") success = true;
		i++;
	}
	return success;
}
bool Molecule::IsTagged() {
	bool success=false;
	tag_segment = -1; 
	int seg_nr=0;
	int NumSegs=In[0]->MonList.size(); 
	int length=MolMonList.size();
	int i=0;
	while (i<length) {
		string segmentname=MolMonList[i]; 
		for (int j=0; j<NumSegs; j++) {if (Seg[j]->name==segmentname) seg_nr=j;}
		if (Seg[seg_nr]->GetFreedom()=="tagged") {success = true; tag_segment=seg_nr;}
		i++;
	}
	return success;
}
bool Molecule::IsCharged() {
	float charge =0;
	int length = n_mon.size(); 
	int i=0;
	while (i<length) {
		charge +=n_mon[i]*Seg[mon_nr[i]]->valence; 
		i++;
	}
	return charge<-1e-5 || charge > 1e-5; 
}
void Molecule::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

string Molecule::GetValue(string parameter) {
	int length = PARAMETERS.size(); 
	int i=0;
	while (i<length) {
		if (PARAMETERS[i]==parameter) { return VALUES[i];}		
		i++;
	} 
	return ""; 
}

void Molecule:: PrepareForCalculations() {
	M=(MX+2)*(MY+2)*(MX+2);

//define on CPU
	H_phi= new float[M*MolMonList.size()]; 

#ifdef CUDA
//define on GPU
	phi=(float*)AllOnDev(M*MolMonList.size());
	Gg_f=(float*)AllOnDev(M*chainlength);
	Gg_b=(float*)AllOnDev(M*2);
#else
//set ref for the rho equal to H_rho etc. or define the ref if not done above.
	phi = H_phi;
	Gg_f= new float[M*chainlength];
	Gg_b= new float[M*2];
#endif

}
 
