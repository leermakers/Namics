#include "segment.h" 
#include <fstream>
Segment::Segment(vector<Input*> In_,vector<Lattice*> Lat_, string name_,int segnr,int N_seg) {
	In=In_; Lat=Lat_; name=name_; n_seg=N_seg; seg_nr=segnr; 
if (debug) cout <<"Segment constructor" + name << endl;
	KEYS.push_back("freedom"); 
	KEYS.push_back("valence");
	KEYS.push_back("epsilon");;
	KEYS.push_back("pinned_range");
	KEYS.push_back("frozen_range");
	KEYS.push_back("tagged_range");
	KEYS.push_back("pinned_filename"); 	
	KEYS.push_back("frozen_filename");
	KEYS.push_back("tagged_filename"); 
	KEYS.push_back("clamp_filename");
}
Segment::~Segment() {
if (debug) cout <<"Segment destructor " + name << endl;
	DeAllocateMemory();
}

void Segment::DeAllocateMemory(void){
if (debug) cout << "In Segment, Deallocating memory"<< endl;
	if(n_pos>0) free(H_P);
	if (freedom != "free") free(r);
	free(H_u); free(H_phi); free(H_MASK);
#ifdef CUDA
	if(n_pos>0) cudaFree(P);
	cudaFree(u); cudaFree(phi); cudaFree(G1); cudaFree(MASK); cudaFree(phi_side);
#else
	free(G1); free(phi_side);
#endif
}

void Segment::AllocateMemory() {
if (debug) cout <<"Allocate Memory in Segment " + name << endl; 
	int M=Lat[0]->M;

	phibulk =0; //initialisatie van monomer phibulk.
	//r=(int*) malloc(6*sizeof(int));
	//if (n_pos>0) {
	//	H_P = (int*) malloc(n_pos*sizeof(int));

	//}
	H_u = (Real*) malloc(M*sizeof(Real));
	H_phi = (Real*) malloc(M*sizeof(Real));
	H_Zero(H_u,M); 
	H_Zero(H_phi,M);
	if (freedom=="free") {
		H_MASK = (int*) malloc(M*sizeof(int));
		H_Zero(H_MASK,M); 
	}
#ifdef CUDA
	if (n_pos>0) {
		Px=(int*)AllIntOnDev(n_pos);
	}
	G1=(Real*)AllOnDev(M); Zero(G1,M);
	u=(Real*)AllOnDev(M); Zero(u,M); 
	MASK=(int*)AllIntOnDev(M); Zero(MASK,M);
	phi=(Real*)AllOnDev(M); Zero(phi,M);
	phi_side=(Real*)AllOnDev(M); Zero(phi_side,M);
#else
	if (n_pos>0) {
		P=H_P;
	}	
	MASK=H_MASK;  
	phi =H_phi;
	u = H_u;
	G1 = (Real*) malloc(M*sizeof(Real));
	phi_side = (Real*) malloc(M*sizeof(Real));
	Zero(G1,M);
	Zero(phi_side,M); 	
#endif
}
bool Segment::PrepareForCalculations(int* KSAM) {
if (debug) cout <<"PrepareForCalcualtions in Segment " +name << endl;

	int M=Lat[0]->M; 
#ifdef CUDA
	TransferIntDataToDevice(H_MASK, MASK, M);
	//TransferDataToDevice(H_u, u, M);
	//TransferIntDataToDevice(H_Px, Px, n_pos);
	//TransferIntDataToDevice(H_Py, Py, n_pos);
	//TransferIntDataToDevice(H_Pz, Pz, n_pos);
#endif

	bool success=true; 
	phibulk=0;
	if (freedom=="frozen") {
		Cp(phi,MASK,M); 
	} else Zero(phi,M); 
	if (freedom=="tagged") Zero(u,M); 
	Boltzmann(G1,u,M);
	if (freedom=="pinned") Times(G1,G1,MASK,M);
	if (freedom=="tagged") Cp(G1,MASK,M);
	Lat[0]->set_bounds(G1);
	if (!(freedom ==" frozen" || freedom =="tagged")) Times(G1,G1,KSAM,M); 
	return success;
}

bool Segment::CheckInput(int start) {
if (debug) cout <<"CheckInput in Segment " + name << endl;
	bool success;
	string s; 
	vector<string>options; 
	success = In[0]->CheckParameters("mon",name,start,KEYS,PARAMETERS,VALUES);
	if(success) {
		valence=In[0]->Get_Real(GetValue("valence"),0);
		if (valence<-5 || valence > 5) cout <<"For mon " + name + " valence value out of range -5 .. 5. Default value used instead" << endl; 
		epsilon=In[0]->Get_Real(GetValue("epsilon"),80);
		if (epsilon<1 || epsilon > 250) cout <<"For mon " + name + " relative epsilon value out of range 1 .. 255. Default value 80 used instead" << endl;

		options.push_back("free"); 
		options.push_back("pinned");
		options.push_back("frozen");
		options.push_back("tagged");
		options.push_back("clamp"); 

		
		freedom = In[0]->Get_string(GetValue("freedom"),"free"); 
		if (!In[0]->InSet(options,freedom)) cout << "Freedom for mon " + name + " not recognized. Default value 'freedom' is used"<< endl;

		if (freedom =="free") {
			if (GetValue("frozen_range").size()>0||GetValue("pinned_range").size()>0 || GetValue("tagged_range").size()>0 ||
			GetValue("frozen_filename").size()>0 || GetValue("pinned_filename").size()>0 || GetValue("tagged_filename").size()>0) {
				success=false; cout <<"In mon " + name + " you should not combine 'freedom : free' with 'frozen_range' or 'pinned_range' or 'tagged_range' or corresponding filenames." << endl; 
			}
		} else {
			r=(int*) malloc(6*sizeof(int));
		}

		if (freedom =="clamps" ) {
			
			if (GetValue("clamp_filename").size()>0) {
				if (!GetClamp(GetValue("clamp_filename"))) {
					success=false; cout <<"Failed to read 'clamp_filename'. Problem terminated" << endl; 
				}
			} else {success=false; cout << "When freedom of segment is 'clamps' we expect to find the quantity 'clamp_filename' which refers to a file that contains the data for the sub_box with clamping coordinates " << endl; } 
		}
		
	
		if (freedom == "pinned") {
			phibulk=0;
			if (GetValue("frozen_range").size()>0 || GetValue("tagged_range").size()>0 || GetValue("frozen_filename").size()>0 || GetValue("tag_filename").size()>0) {
			cout<< "For mon " + name + ", you should exclusively combine 'freedom : pinned' with 'pinned_range' or 'pinned_filename'" << endl;  success=false;}
			if (GetValue("pinned_range").size()>0 && GetValue("pinned_filename").size()>0) {
				cout<< "For mon " + name + ", you can not combine 'pinned_range' with 'pinned_filename' " <<endl; success=false;
			}
			if (GetValue("pinned_range").size()==0 && GetValue("pinned_filename").size()==0) {
				cout<< "For mon " + name + ", you should provide either 'pinned_range' or 'pinned_filename' " <<endl; success=false;
			}
			if (GetValue("pinned_range").size()>0) { s="pinned_range";
				n_pos=0;
				if (success) success=Lat[0]->ReadRange(r, H_P, n_pos, block, GetValue("pinned_range"),name,s);
				if (n_pos>0) {
					H_P=(int*) malloc(n_pos*sizeof(int));
					if (success) success=Lat[0]->ReadRange(r, H_P, n_pos, block, GetValue("pinned_range"),name,s);
				}
			}
			if (GetValue("pinned_filename").size()>0) { s="pinned"; 
				filename=GetValue("pinned_filename"); 
				n_pos=0;
				if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s);
				if (n_pos>0) {
					H_P=(int*) malloc(n_pos*sizeof(int));
					if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s);
				}
			}
		}

		if (freedom == "frozen") {
			phibulk=0;
			if (GetValue("pinned_range").size()>0 || GetValue("tagged_range").size()>0 || GetValue("pinned_filename").size()>0 || GetValue("tag_filename").size()>0) {
			cout<< "For mon " + name + ", you should exclusively combine 'freedom : frozen' with 'frozen_range' or 'frozen_filename'" << endl;  success=false;}
			if (GetValue("frozen_range").size()>0 && GetValue("frozen_filename").size()>0) {
				cout<< "For mon " + name + ", you can not combine 'frozen_range' with 'frozen_filename' " <<endl; success=false;
			}
			if (GetValue("frozen_range").size()==0 && GetValue("frozen_filename").size()==0) {
				cout<< "For mon " + name + ", you should provide either 'frozen_range' or 'frozen_filename' " <<endl; success=false;
			}
			if (GetValue("frozen_range").size()>0) { s="frozen_range";
				n_pos=0;
				success=Lat[0]->ReadRange(r, H_P, n_pos, block, GetValue("frozen_range"),name,s);
				if (n_pos>0) {
					H_P=(int*) malloc(n_pos*sizeof(int));
					success=Lat[0]->ReadRange(r, H_P, n_pos, block, GetValue("frozen_range"),name,s);
				}
			}
			if (GetValue("frozen_filename").size()>0) { s="frozen";
				filename=GetValue("frozen_filename");
				n_pos=0;
				if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s);
				if (n_pos>0) {
					H_P=(int*) malloc(n_pos*sizeof(int));
					if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s);
				}
			} 			 
		} 
	
		if (freedom == "tagged") { 
			phibulk=0;
			if (GetValue("pinned_range").size()>0 || GetValue("frozen_range").size()>0 || GetValue("pinned_filename").size()>0 || GetValue("frozen_filename").size()>0) {
			cout<< "For mon " + name + ", you should exclusively combine 'freedom : tagged' with 'tagged_range' or 'tagged_filename'" << endl;  success=false;}
			if (GetValue("tagged_range").size()>0 && GetValue("tagged_filename").size()>0) {
				cout<< "For mon " + name + ", you can not combine 'tagged_range' with 'tagged_filename' " <<endl; success=false;
			}
			if (GetValue("tagged_range").size()==0 && GetValue("tagged_filename").size()==0) {
				cout<< "For mon " + name + ", you should provide either 'tagged_range' or 'tagged_filename' " <<endl; success=false;
			}
			if (GetValue("tagged_range").size()>0) { s="tagged_range";
				n_pos=0;
				if (success) success=Lat[0]->ReadRange(r, H_P, n_pos, block, GetValue("tagged_range"),name,s);
				if (n_pos>0) {
					H_P=(int*) malloc(n_pos*sizeof(int));
					if (success) success=Lat[0]->ReadRange(r, H_P, n_pos, block, GetValue("tagged_range"),name,s);
				}
			} 
			if (GetValue("tagged_filename").size()>0) {s="tagged";
				filename=GetValue("tagged_filename");
				n_pos=0;
			 	if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s);
				if (n_pos>0) {
					H_P=(int*) malloc(n_pos*sizeof(int));
			 		if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s);
				}
			} 	
		}
	
		if (freedom!="free") { 
			H_MASK = (int*) malloc(Lat[0]->M*sizeof(int));
			Lat[0]->CreateMASK(H_MASK,r,H_P,n_pos,block); 
		}

	}
	
	return success; 
}

bool Segment::GetClamp(string filename) {
	bool success=true;
	string line_;
	string NN,X,Y,Z;
	int pos;
	int N=0;
	mx=my=mz=0;
	n_box=0;
	int i=0;
	ifstream in_file;
	in_file.open(filename.c_str());
	if (in_file.is_open()) {
		while (in_file) {
			n_box++; i++; 
			if (i==1) {
				in_file >> line_ >> NN >> line_ >> X >> Y >> Z ;				
				if (!In[0]->Get_int(NN,pos,0)) { success=false; cout<<" length of 'fragment' not an integer " << endl; }
				else {if (N>0) {if (N!=pos) {cout <<"lengths of fregment are not equal " << endl; success=false;}} else N = pos;}
				if (!In[0]->Get_int(X,pos,0)) {success=false; cout <<"X-value for sub_box size is not integer " << endl;  } 
				else {if (mx>0) {if (mx !=pos) {success=false; cout <<"We can deal with only one sub-box size" << endl;}} else mx=pos; }
				if (!In[0]->Get_int(Y,pos,0)) {success=false; cout <<"Y-value for sub_box size is not integer " << endl;  } 
				else {if (my>0) {if (my !=pos) {success=false; cout <<"We can deal with only one sub-box size" << endl;}} else my=pos; }
				if (!In[0]->Get_int(Z,pos,0)) {success=false; cout <<"Z-value for sub_box size is not integer " << endl;  } 
				else {if (mz>0) {if (mz !=pos) {success=false; cout <<"We can deal with only one sub-box size" << endl;}} else mz=pos; }
			}
			if (i==2) {in_file>>X>>Y>>Z;
				if (!In[0]->Get_int(X,pos,0)) {success =false; cout << " X-pos of particle sub_box nr " << n_box << " not integer"  << endl; }
				else bx.push_back(pos);
				if (!In[0]->Get_int(Y,pos,0)) {success =false; cout << " Y-pos of particle sub_box nr " << n_box << " not integer"  << endl; }
				else by.push_back(pos);
				if (!In[0]->Get_int(Z,pos,0)) {success =false; cout << " Z-pos of particle sub_box nr " << n_box << " not integer"  << endl; }
				else bz.push_back(pos);
			}
			if (i==3) {in_file>>X>>Y>>Z;
				if (!In[0]->Get_int(X,pos,0)) {success =false; cout << " X-pos of particle pos 1 of sub_box nr " << n_box << " not integer"  << endl; }
				else px1.push_back(pos);
				if (!In[0]->Get_int(Y,pos,0)) {success =false; cout << " Y-pos of particle pos 1 of sub_box nr " << n_box << " not integer"  << endl; }
				else py1.push_back(pos);
				if (!In[0]->Get_int(Z,pos,0)) {success =false; cout << " Z-pos of particle pos 1 of sub_box nr " << n_box << " not integer"  << endl; }
				else pz1.push_back(pos);
			}
			if (i==4) {i=0;in_file>>X>>Y>>Z;
				if (!In[0]->Get_int(X,pos,0)) {success =false; cout << " X-pos of particle pos 2 of sub_box nr " << n_box << " not integer"  << endl; }
				else px2.push_back(pos);
				if (!In[0]->Get_int(Y,pos,0)) {success =false; cout << " Y-pos of particle pos 2 of sub_box nr " << n_box << " not integer"  << endl; }
				else py2.push_back(pos);
				if (!In[0]->Get_int(Z,pos,0)) {success =false; cout << " Z-pos of particle pos 2 of sub_box nr " << n_box << " not integer"  << endl; }
				else pz2.push_back(pos);
			}
		}
		
	} else {
		success=false;
		cout <<"To read 'clamp_file', the file " << filename << " is not available." << endl; 
	}
	in_file.close();
	return success; 
}

int* Segment::GetMASK() {
if (debug) cout <<"Get Mask for segment" + name << endl;
	if (MASK==NULL) {cout <<"MASK not yet created. Task to point to MASK in segment is rejected. " << endl; return NULL;} 
	else return MASK; 
}

Real* Segment::GetPhi() {
if (debug) cout <<"GetPhi in segment " + name << endl;
	int M=Lat[0]->M;
	if (freedom=="frozen") Cp(phi,MASK,M); 
	return phi; 
}

string Segment::GetFreedom(void){
if (debug) cout <<"GetFreedom for segment " + name << endl;
	return freedom; 
}
bool Segment::IsClamp(void) {
if (debug) cout <<"Is free for " + name << endl;
	return freedom == "clamp"; 
}
bool Segment::IsFree(void) {
if (debug) cout <<"Is free for " + name << endl;
	return freedom == "free"; 
}
bool Segment::IsPinned(void) {
if (debug) cout <<"IsPinned for segment " + name << endl;
	return freedom == "pinned"; 
}
bool Segment::IsFrozen(void) {
if (debug) cout <<"IsFrozen for segment " + name << endl;
	phibulk =0;
	return freedom == "frozen"; 
}
bool Segment::IsTagged(void) {
if (debug) cout <<"IsTagged for segment " + name << endl;
	phibulk =0;
	return freedom == "tagged"; 
}

void Segment::PutChiKEY(string new_name) {
if (debug) cout <<"PutChiKey " + name << endl;
	if(name != new_name) KEYS.push_back("chi-" + new_name); 
}

string Segment::GetValue(string parameter) {
if (debug) cout <<"GetValue for segment " + name + " for parameter " + parameter << endl;
	int length = PARAMETERS.size(); 
	int i=0;
	while (i<length) {
		if (PARAMETERS[i]==parameter) { return VALUES[i];}		
		i++;
	} 
	return ""; 
}

void Segment::push(string s, Real X) {
if (debug) cout <<"Push in Segment (Real) " + name << endl;
	Reals.push_back(s);
	Reals_value.push_back(X); 
}
void Segment::push(string s, int X) {
if (debug) cout <<"Push in Segment (int) " + name << endl;
	ints.push_back(s);
	ints_value.push_back(X); 
}
void Segment::push(string s, bool X) {
if (debug) cout <<"Push in Segment (bool) " + name << endl;
	bools.push_back(s);
	bools_value.push_back(X); 
}
void Segment::push(string s, string X) {
if (debug) cout <<"Push in Segment (string) " + name << endl;
	strings.push_back(s);
	strings_value.push_back(X); 	
}
void Segment::PushOutput() {
if (debug) cout <<"PushOutput for segment " + name << endl;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();  
	push("freedom",freedom);
	if (freedom=="pinned") push("range",GetValue("pinned_range"));
	if (freedom=="frozen") push("range",GetValue("frozen_range"));
	if (freedom=="tagged") push("range",GetValue("tagged_range"));
	string profile="profile;0"; push("phi",profile);
	profile="profile;1"; push("u",profile);
	//profile="profile;2"; push("mask",profile);
#ifdef CUDA
	int M = Lat[0]->M;
	TransferDataToHost(H_phi, phi, M);
	TransferDataToHost(H_u, u, M);
#endif
}

Real* Segment::GetPointer(string s) {
if (debug) cout <<"Get Pointer for segment " + name << endl;
	vector<string> sub;
	In[0]->split(s,';',sub);
	if (sub[1]=="0") return H_phi;
	if (sub[1]=="1") return H_u;
	//if (sub[1]=="2") return H_MASK;
	return NULL; 
}

int Segment::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
if (debug) cout <<"GetValue long for segment " + name << endl;
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


/*
	for (int p=0; p<n_box; p++){
		H_mask[p*M + jx*(H_Px[p]-H_Bx[p])+jy*(H_Py[p]-H_By[p])+(H_Pz[p]-H_Bz[p])]=1;
		H_MASK[((H_Px[p]-1)%MX+1)*JX + ((H_Py[p]-1)%MY+1)*JY + (H_Pz[p]-1)%MZ+1]=1;
	}

	H_mask = new Real[M*n_box];
#ifdef CUDA
	mask= (Real*)AllOnDev(M*n_box);
#else
	mask=H_mask;
#endif

	H_Zero(H_mask,M*n_box);
	H_Bx  = new int[n_box];  H_By  = new int[n_box];   H_Bz  = new int[n_box];
	H_Px = new int[n_box];  H_Py = new int[n_box];   H_Pz = new int[n_box];
#ifdef CUDA
	Bx=(int*)AllIntOnDev(n_box);
	By=(int*)AllIntOnDev(n_box);
	Bz=(int*)AllIntOnDev(n_box);
	Px=(int*)AllIntOnDev(n_box);
	Py=(int*)AllIntOnDev(n_box);
	Pz=(int*)AllIntOnDev(n_box);

#else
	Bx=H_Bx; By=H_By; Bz=H_Bz;
	Px=H_Px; Py=H_Py; Pz=H_Pz;

#endif
*/
