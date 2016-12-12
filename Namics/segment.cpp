class Segment {
public:
	Segment(vector<Input*>,string,int,int,int,int,int);

~Segment();

	string name; 
	vector<Input*> In;  
	int n_seg; 
	int seg_nr;  
	double epsilon;
	double valence; 
	double phibulk; 
	string freedom;

	string filename; 
	bool block;
	int n_pos;  
	int MX,MY,MZ;
	int M; 
	int xl,yl,zl,xh,yh,zh;
	int* H_Px; 
	int* H_Py;
	int* H_Pz;
	float* H_MASK;
	float* H_KSAM; 
	float* H_rho;
	float* H_u; //is this necessary?
	float* H_G1;
	float* H_phi;
	int* Px; 
	int* Py;
	int* Pz;
	float* MASK;
	float* KSAM; 
	float* rho;
	float* G1;
	float* phi;
	float* phi_side;  


	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES; 
	bool CheckInput(); 
	void PutChiKEY(string); 
	string GetValue(string); 
	string GetFreedom();
	bool IsFree();
	bool IsPinned();
	bool IsFrozen();
	bool IsTagged(); 
	bool CreateMASK();
	float* GetMASK(); 
	void PrepareForCalculations();  

};
Segment::Segment(vector<Input*> In_, string name_,int segnr,int N_seg, int MX_, int MY_, int MZ_) {
	In=In_, name=name_; n_seg=N_seg; seg_nr=segnr; MX=MX_; MY=MY_; MZ=MZ_; 
	KEYS.push_back("freedom"); 
	KEYS.push_back("valence");
	KEYS.push_back("epsilon");;
	KEYS.push_back("pinned_range");
	KEYS.push_back("frozen_range");
	KEYS.push_back("tagged_range");
	KEYS.push_back("pinned_filename"); 	
	KEYS.push_back("frozen_filename");
	KEYS.push_back("tag_filename"); 
}
Segment::~Segment() {
}

bool Segment::CheckInput() {
	bool success;
	vector<std::string>set;
	vector<std::string>xyz;
	vector<std::string>coor; 
	vector<string>options; 
	success = In[0]->CheckParameters("mon",name,KEYS,PARAMETERS,VALUES);
	if(success) {
		valence=In[0]->Get_float(GetValue("valence"),0);
		if (valence<-5 || valence > 5) cout <<"For mon " + name + " valence value out of range -5 .. 5. Default value used instead" << endl; 
		epsilon=In[0]->Get_float(GetValue("epsilon"),80);
		if (epsilon<1 || epsilon > 250) cout <<"For mon " + name + " relative epsilon value out of range 1 .. 255. Default value 80 used instead" << endl;

		options.push_back("free"); options.push_back("pinned");options.push_back("frozen");options.push_back("tagged");
		freedom = In[0]->Get_string(GetValue("freedom"),"free"); 
		if (!In[0]->InSet(options,freedom)) cout << "Freedom for mon " + name + " not recognized. Default value 'freedom' is used"<< endl;
//--------------------------free
		if (freedom =="free") {
			if (GetValue("frozen_range").size()>0||GetValue("pinned_range").size()>0 || GetValue("tagged_range").size()>0 ||
			GetValue("frozen_filename").size()>0 || GetValue("pinned_filename").size()>0 || GetValue("tagged_filename").size()>0) {
				success=false; cout <<"In mon " + name + " you should not combine 'freedom : free' with 'frozen_range' or 'pinned_range' or 'tagged_range' or xxx_filename." << endl; 
			}
		}
//--------------------------pinned		
		if (freedom == "pinned") {
			if (GetValue("frozen_range").size()>0 || GetValue("tagged_range").size()>0 || GetValue("frozen_filename").size()>0 || GetValue("tag_filename").size()>0) {
			cout<< "For mon " + name + ", you should exclusively combine 'freedom : pinned' with 'pinned_range' or 'pinned_filename'" << endl;  success=false;}
			if (GetValue("pinned_range").size()>0 && GetValue("pinned_filename").size()>0) {
				cout<< "For mon " + name + ", you can not combine 'pinned_range' with 'pinned_filename' " <<endl; success=false;
			}
			if (GetValue("pinned_range").size()==0 && GetValue("pinned_filename").size()==0) {
				cout<< "For mon " + name + ", you should provide either 'pinned_range' or 'pinned_filename' " <<endl; success=false;
			}
			if (GetValue("pinned_range").size()>0) { 
				In[0]->split(GetValue("pinned_range"),';',set);
				if (set.size()==2) { coor.clear(); 
					block=true; In[0]->split(set[0],',',coor);

					if (coor.size()!=3) {cout << "In mon " + name + ", for 'pos 1', in 'pinned_range' the coordiantes do not come in three: (x,y,z)" << endl; success=false;}
					else {

					xl=In[0]->Get_int(coor[0],0);
					if (xl < 1 || xl > MX) {cout << "In mon " + name + ", for 'pos 1', the x-coordinate in 'pinned_range' is out of bounds: 1.." << MX << endl; success =false;}
					yl=In[0]->Get_int(coor[1],0);
					if (yl < 1 || yl > MY) {cout << "In mon " + name + ", for 'pos 1', the y-coordinate in 'pinned_range' is out of bounds: 1.." << MY << endl; success =false;}
					zl=In[0]->Get_int(coor[2],0);
					if (zl < 1 || zl > MZ) {cout << "In mon " + name + ", for 'pos 1', the z-coordinate in 'pinned_range' is out of bounds: 1.." << MZ << endl; success =false;} 
					}
					coor.clear(); In[0]->split(set[1],',',coor);

					if (coor.size()!=3) {cout << "In mon " + name + ", for 'pos 2', in 'pinned_range', the coordinates do not come in three: (x,y,z)" << endl; success=false;}
					else {


					xh=In[0]->Get_int(coor[0],0);
					if (xh < 1 || xh > MX) {cout << "In mon " + name + ", for 'pos 2', the x-coordinate in 'pinned_range' is out of bounds; 1.." << MX << endl; success =false;}
					yh=In[0]->Get_int(coor[1],0);
					if (yh < 1 || yh > MY) {cout << "In mon " + name + ", for 'pos 2', the y-coordinate in 'pinned_range' is out of bounds; 1.." << MY << endl; success =false;}
					zh=In[0]->Get_int(coor[2],0);
					if (zh < 1 || xh > MZ) {cout << "In mon " + name + ", for 'pos 2', the z-coordinate in 'pinned_range' is out of bounds; 1.." << MZ << endl; success =false;}
					if (xl > xh) {cout << "In mon " + name + ", for 'pos 1', the x-coordinate in 'pinned_range' should be less than that of 'pos 2'" << endl; success =false;}
					if (yl > yh) {cout << "In mon " + name + ", for 'pos 1', the y-coordinate in 'pinned_range' should be less than that of 'pos 2'" << endl; success =false;}
					if (zl > zh) {cout << "In mon " + name + ", for 'pos 1', the z-coordinate in 'pinned_range' should be less than that of 'pos 2'" << endl; success =false;}
					}
				} else {
					block=false;  
					In[0]->split(set[0],')',coor);
					int length=coor.size(); n_pos=length; 
					H_Px = new int[length];
					H_Py = new int[length]; 
					H_Pz = new int[length]; 
					int i=0;	
					while (i<length) { 
						string s=coor[i].substr(1,coor[i].size()-1); 
						In[0]->split(s,',',xyz);
						int length_xyz=xyz.size();
						if (length_xyz!=3) { 
							cout << "In mon " + name + " pinned_range  the expected 'triple coordinate' structure (x,y,z) was not found. " << endl;  success = false;
						} else {   
							H_Px[i]=In[0]->Get_int(xyz[0],0);
							if (H_Px[i] < 1 || H_Px[i] > MX) {cout << "In mon " + name + ", for 'pos' "<< i << ", the x-coordinate in pinned_range out of bounds: 1.." << MX << endl; success =false;}
							H_Py[i]=In[0]->Get_int(xyz[1],0);
							if (H_Py[i] < 1 || H_Py[i] > MY) {cout << "In mon " + name + ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 1.." << MY << endl; success =false;}								
							H_Pz[i]=In[0]->Get_int(xyz[2],0);
							if (H_Pz[i] < 1 || H_Pz[i] > MZ) {cout << "In mon " + name + ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 1.." << MZ << endl; success =false;}	
						}
						i++;
					}
				}
			} 
			if (GetValue("pinned_filename").size()>0) filename=GetValue("pinned_filename"); 			 
		}//pinned
//--------------------------frozen

if (freedom == "frozen") {
			if (GetValue("pinned_range").size()>0 || GetValue("tagged_range").size()>0 || GetValue("pinned_filename").size()>0 || GetValue("tag_filename").size()>0) {
			cout<< "For mon " + name + ", you should exclusively combine 'freedom : frozen' with 'frozen_range' or 'frozen_filename'" << endl;  success=false;}
			if (GetValue("frozen_range").size()>0 && GetValue("frozen_filename").size()>0) {
				cout<< "For mon " + name + ", you can not combine 'frozen_range' with 'frozen_filename' " <<endl; success=false;
			}
			if (GetValue("frozen_range").size()==0 && GetValue("frozen_filename").size()==0) {
				cout<< "For mon " + name + ", you should provide either 'frozen_range' or 'frozen_filename' " <<endl; success=false;
			}
			if (GetValue("frozen_range").size()>0) { 
				In[0]->split(GetValue("frozen_range"),';',set);
				if (set.size()==2) { coor.clear(); 
					n_pos=2; 
					H_Px = new int[2];
					H_Py = new int[2]; 
					H_Pz = new int[2]; 
					block=true; In[0]->split(set[0],',',coor);

					if (coor.size()!=3) {cout << "In mon " + name + ", for 'pos 1', in 'frozen_range' the coordiantes do not come in three: (x,y,z)" << endl; success=false;}
					else {

					H_Px[0]=In[0]->Get_int(coor[0],0);
					if (H_Px[0] < 1 || H_Px[0] > MX) {cout << "In mon " + name + ", for 'pos 1', the x-coordinate in 'frozen_range' is out of bounds: 1.." << MX << endl; success =false;}
					H_Py[0]=In[0]->Get_int(coor[1],0);
					if (H_Py[0] < 1 || H_Py[0] > MY) {cout << "In mon " + name + ", for 'pos 1', the y-coordinate in 'frozen_range' is out of bounds: 1.." << MY << endl; success =false;}
					H_Pz[0]=In[0]->Get_int(coor[2],0);
					if (H_Pz[0] < 1 || H_Pz[0] > MZ) {cout << "In mon " + name + ", for 'pos 1', the z-coordinate in 'frozen_range' is out of bounds: 1.." << MZ << endl; success =false;} 
					}
					coor.clear(); In[0]->split(set[1],',',coor);

					if (coor.size()!=3) {cout << "In mon " + name + ", for 'pos 2', in 'frozen_range', the coordinates do not come in three: (x,y,z)" << endl; success=false;}
					else {


					H_Px[1]=In[0]->Get_int(coor[0],0);
					if (H_Px[1] < 1 || H_Px[1] > MX) {cout << "In mon " + name + ", for 'pos 2', the x-coordinate in 'frozen_range' is out of bounds; 1.." << MX << endl; success =false;}
					H_Py[1]=In[0]->Get_int(coor[1],0);
					if (H_Py[1] < 1 || H_Py[1] > MY) {cout << "In mon " + name + ", for 'pos 2', the y-coordinate in 'frozen_range' is out of bounds; 1.." << MY << endl; success =false;}
					H_Pz[1]=In[0]->Get_int(coor[2],0);
					if (H_Pz[1] < 1 || H_Pz[1] > MZ) {cout << "In mon " + name + ", for 'pos 2', the z-coordinate in 'frozen_range' is out of bounds; 1.." << MZ << endl; success =false;}
					if (H_Px[0] >  H_Px[1]) {cout << "In mon " + name + ", for 'pos 1', the x-coordinate in 'frozen_range' should be less than that of 'pos 2'" << endl; success =false;}
					if (H_Py[0] >  H_Py[1]) {cout << "In mon " + name + ", for 'pos 1', the y-coordinate in 'frozen_range' should be less than that of 'pos 2'" << endl; success =false;}
					if (H_Pz[0] >  H_Pz[1]) {cout << "In mon " + name + ", for 'pos 1', the z-coordinate in 'frozen_range' should be less than that of 'pos 2'" << endl; success =false;}
					}
				} else {
					block=false;  
					In[0]->split(set[0],')',coor);
					int length=coor.size();
					H_Px = new int[length];
					H_Py = new int[length]; 
					H_Pz = new int[length]; 
					int i=0;	
					while (i<length) { 
						string s=coor[i].substr(1,coor[i].size()-1); 
						In[0]->split(s,',',xyz);
						int length_xyz=xyz.size();
						if (length_xyz!=3) { 
							cout << "In mon " + name + " frozen_range  the expected 'triple coordinate' structure (x,y,z) was not found. " << endl;  success = false;
						} else {   
							H_Px[i]=In[0]->Get_int(xyz[0],0);
							if (H_Px[i] < 1 || H_Px[i] > MX) {cout << "In mon " + name + ", for 'pos' "<< i << ", the x-coordinate in frozen_range out of bounds: 1.." << MX << endl; success =false;}
							H_Py[i]=In[0]->Get_int(xyz[1],0);
							if (H_Py[i] < 1 || H_Py[i] > MY) {cout << "In mon " + name + ", for 'pos' "<< i << ", the y-coordinate in frozen_range out of bounds: 1.." << MY << endl; success =false;}								
							H_Pz[i]=In[0]->Get_int(xyz[2],0);
							if (H_Pz[i] < 1 || H_Pz[i] > MZ) {cout << "In mon " + name + ", for 'pos' "<< i << ", the y-coordinate in frozen_range out of bounds: 1.." << MZ << endl; success =false;}	
						}
						i++;
					}
				}
			} 
			if (GetValue("frozen_filename").size()>0) filename=GetValue("frozen_filename"); 			 
		}//frozen

//------------------------- tagged	 

if (freedom == "tagged") { 
			if (GetValue("pinned_range").size()>0 || GetValue("frozen_range").size()>0 || GetValue("pinned_filename").size()>0 || GetValue("frozen_filename").size()>0) {
			cout<< "For mon " + name + ", you should exclusively combine 'freedom : tagged' with 'tagged_range' or 'tagged_filename'" << endl;  success=false;}
			if (GetValue("tagged_range").size()>0 && GetValue("tagged_filename").size()>0) {
				cout<< "For mon " + name + ", you can not combine 'tagged_range' with 'tagged_filename' " <<endl; success=false;
			}
			if (GetValue("tagged_range").size()==0 && GetValue("tagged_filename").size()==0) {
				cout<< "For mon " + name + ", you should provide either 'tagged_range' or 'tagged_filename' " <<endl; success=false;
			}
			if (GetValue("tagged_range").size()>0) { 
				In[0]->split(GetValue("tagged_range"),';',set);
				if (set.size()==2) { coor.clear(); 
					n_pos=2; 
					H_Px = new int[2];
					H_Py = new int[2]; 
					H_Pz = new int[2]; 
					block=true; In[0]->split(set[0],',',coor);

					if (coor.size()!=3) {cout << "In mon " + name + ", for 'pos 1', in 'tagged_range' the coordiantes do not come in three: (x,y,z)" << endl; success=false;}
					else {

					H_Px[0]=In[0]->Get_int(coor[0],0);
					if (H_Px[0] < 1 || H_Px[0] > MX) {cout << "In mon " + name + ", for 'pos 1', the x-coordinate in 'tagged_range' is out of bounds: 1.." << MX << endl; success =false;}
					H_Py[0]=In[0]->Get_int(coor[1],0);
					if (H_Py[0] < 1 || H_Py[0] > MY) {cout << "In mon " + name + ", for 'pos 1', the y-coordinate in 'tagged_range' is out of bounds: 1.." << MY << endl; success =false;}
					H_Pz[0]=In[0]->Get_int(coor[2],0);
					if (H_Pz[0] < 1 || H_Pz[0] > MZ) {cout << "In mon " + name + ", for 'pos 1', the z-coordinate in 'tagged_range' is out of bounds: 1.." << MZ << endl; success =false;} 
					}
					coor.clear(); In[0]->split(set[1],',',coor);

					if (coor.size()!=3) {cout << "In mon " + name + ", for 'pos 2', in 'tagged_range', the coordinates do not come in three: (x,y,z)" << endl; success=false;}
					else {


					H_Px[1]=In[0]->Get_int(coor[0],0);
					if (H_Px[1] < 1 || H_Px[1] > MX) {cout << "In mon " + name + ", for 'pos 2', the x-coordinate in 'tagged_range' is out of bounds; 1.." << MX << endl; success =false;}
					H_Py[1]=In[0]->Get_int(coor[1],0);
					if (H_Py[1] < 1 || H_Py[1] > MY) {cout << "In mon " + name + ", for 'pos 2', the y-coordinate in 'tagged_range' is out of bounds; 1.." << MY << endl; success =false;}
					H_Pz[1]=In[0]->Get_int(coor[2],0);
					if (H_Pz[1] < 1 || H_Pz[1] > MZ) {cout << "In mon " + name + ", for 'pos 2', the z-coordinate in 'tagged_range' is out of bounds; 1.." << MZ << endl; success =false;}
					if (H_Px[0] >  H_Px[1]) {cout << "In mon " + name + ", for 'pos 1', the x-coordinate in 'tagged_range' should be less than that of 'pos 2'" << endl; success =false;}
					if (H_Py[0] >  H_Py[1]) {cout << "In mon " + name + ", for 'pos 1', the y-coordinate in 'tagged_range' should be less than that of 'pos 2'" << endl; success =false;}
					if (H_Pz[0] >  H_Pz[1]) {cout << "In mon " + name + ", for 'pos 1', the z-coordinate in 'tagged_range' should be less than that of 'pos 2'" << endl; success =false;}
					}
				} else {
					block=false;  
					In[0]->split(set[0],')',coor);
					n_pos=coor.size();
					H_Px = new int[n_pos];
					H_Py = new int[n_pos]; 
					H_Pz = new int[n_pos]; 
					int i=0;	
					while (i<n_pos) { 
						string s=coor[i].substr(1,coor[i].size()-1); 
						In[0]->split(s,',',xyz);
						int length_xyz=xyz.size();
						if (length_xyz!=3) { 
							cout << "In mon " + name + " tagged_range  the expected 'triple coordinate' structure '(x,y,z)' was not found. " << endl;  success = false;
						} else {   
							H_Px[i]=In[0]->Get_int(xyz[0],0);
							if (H_Px[i] < 1 || H_Px[i] > MX) {cout << "In mon " + name + ", for 'pos' "<< i << ", the x-coordinate in tagged_range out of bounds: 1.." << MX << endl; success =false;}
							H_Py[i]=In[0]->Get_int(xyz[1],0);
							if (H_Py[i] < 1 || H_Py[i] > MY) {cout << "In mon " + name + ", for 'pos' "<< i << ", the y-coordinate in tagged_range out of bounds: 1.." << MY << endl; success =false;}								
							H_Pz[i]=In[0]->Get_int(xyz[2],0);
							if (H_Pz[i] < 1 || H_Pz[i] > MZ) {cout << "In mon " + name + ", for 'pos' "<< i << ", the y-coordinate in tagged_range out of bounds: 1.." << MZ << endl; success =false;}	
						}
						i++;
					}
				}
			} 
			if (GetValue("tagged_filename").size()>0) filename=GetValue("tagged_filename"); 			 
		}//tagged
//---------------------------------now read file.
		if (filename.size()>0) {
			string content;
			vector<string> lines;
			vector<string> sub;  
			string Infilename=In[0]->name;
			In[0]->split(Infilename,'.',sub);
			 
			 
			if (!In[0]->ReadFile(sub[0].append(".").append(filename),content)) {success=false;} else {
				In[0]->split(content,'#',lines);
				int length = lines.size();
				if (length == MX*MY*MZ) { //expect to read 'mask file';
					n_pos=0;
					int i=0;
					while (i<length){
						if (In[0]->Get_int(lines[i],0)==1) n_pos++; 
						i++; 
					}
					if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any unities." << endl;}
					else { int p_i; 
						H_Px = new int[n_pos];
						H_Py = new int[n_pos]; 
						H_Pz = new int[n_pos];
						i=0; p_i=-1;
						for (int x=0; x<MX; x++)  for (int y=0; y<MY; y++) for (int z=0; z<MZ; z++) {
							if (In[0]->Get_int(lines[i],0)==1) {p_i++; H_Px[p_i]=x; H_Py[p_i]=y; H_Pz[p_i]=z; }
							i++; 
						}
					}
				} else {
					//expect to read x,y,z
					int i=0;
					n_pos=length; 
					H_Px = new int[n_pos];
					H_Py = new int[n_pos]; 
					H_Pz = new int[n_pos];
					while (i<length) {  
						xyz.clear(); 
						In[0]->split(lines[i],',',xyz);
						int length_xyz=xyz.size();
						if (length_xyz!=3) { 
							cout << "In mon " + name + " 'xxx_filename'  the expected 'triple coordinate' structure 'x,y,z' was not found. " << endl;  success = false;
						} else {   
							H_Px[i]=In[0]->Get_int(xyz[0],0);
							if (H_Px[i] < 1 || H_Px[i] > MX) {cout << "In mon " + name + ", for 'pos' "<< i << ", the x-coordinate in 'xxx_filename' out of bounds: 1.." << MX << endl; success =false;}
							H_Py[i]=In[0]->Get_int(xyz[1],0);
							if (H_Py[i] < 1 || H_Py[i] > MY) {cout << "In mon " + name + ", for 'pos' "<< i << ", the y-coordinate in tagged_range out of bounds: 1.." << MY << endl; success =false;}								
							H_Pz[i]=In[0]->Get_int(xyz[2],0);
							if (H_Pz[i] < 1 || H_Pz[i] > MZ) {cout << "In mon " + name + ", for 'pos' "<< i << ", the y-coordinate in tagged_range out of bounds: 1.." << MZ << endl; success =false;}	
						}
						i++;
					}
					
				}	
			} 
		}
	
	if (!IsFree()) CreateMASK(); 

	}//success
	
//for (int i=0; i<n_pos; i++) cout << "pos[x] = " << H_Px[i] << " pos[y] = " << H_Py[i] << " pos[z] = " << H_Pz[i] << endl; 

	return success; 
}



bool Segment::CreateMASK() {
	bool success=true; 
	if (!IsFree()) {
		int M = (MX+2)*(MY+2)*(MZ+2);
		int JX = (MX+2)*(MY+2);
		int JY = (MY+2);
		if (H_MASK==NULL) {H_MASK = new float[M];} else {cout << "MASK already exists; task to create MASK rejected. " << endl; success=false; }
		H_Zero(H_MASK,M);
		if (block) {
			for (int x=1; x<MX+1; x++) for (int y=1; y<MY+1; y++) for (int z=1; z<MZ+1; z++) 
			if (x >= H_Px[0] && y >= H_Py[0] && z >= H_Pz[0] && x <= H_Px[1] && y <= H_Py[1] && z <= H_Pz[1]) {H_MASK[x*JX+y*JY+z]=1;} else {H_MASK[x*JX+y*JY+z]=0;}
		} else {
			for (int i=0; i<n_pos; i++) H_MASK[H_Px[i]*JX+H_Py[i]*JY+H_Pz[i]]=1;
		}
	} else {
		success=false; 
	}
	return success; 
}

float* Segment::GetMASK() {
	if (MASK==NULL) {cout <<"MASK not yet created. Task to point to MASK in segment is rejected. " << endl; return NULL;} 
	else return MASK; 
}

string Segment::GetFreedom(void){
	return freedom; 
}

bool Segment::IsFree(void) {
	return freedom == "free"; 
}
bool Segment::IsPinned(void) {
	return freedom == "pinned"; 
}
bool Segment::IsFrozen(void) {
	return freedom == "frozen"; 
}
bool Segment::IsTagged(void) {
	return freedom == "tagged"; 
}

void Segment::PutChiKEY(string new_name) {
	if(name != new_name) KEYS.push_back("chi-" + new_name); 
}

string Segment::GetValue(string parameter) {
	int length = PARAMETERS.size(); 
	int i=0;
	while (i<length) {
		if (PARAMETERS[i]==parameter) { return VALUES[i];}		
		i++;
	} 
	return ""; 
}

void Segment::PrepareForCalculations() {
	M=(MX+2)*(MY+2)*(MX+2);
	phibulk =0; //initialisatie van monomer phibulk.

//define on CPU
	if (H_Px==NULL && n_pos>0) H_Px=new int[n_pos];
	if (H_Py==NULL && n_pos>0) H_Py=new int[n_pos];
	if (H_Pz==NULL && n_pos>0) H_Pz=new int[n_pos];
	H_G1= new float[M]; 
	H_rho = new float[M];
	if (H_MASK == NULL) H_MASK=new float[M]; //or H_MASK = new int[]?
	H_KSAM = new float[M];
#ifdef CUDA
//define on GPU
	if (n_pos>0) {
		Px=(int*)AllIntOnDev(n_pos);
		Py=(int*)AllIntOnDev(n_pos);
		Pz=(int*)AllIntOnDev(n_pos);
	}
	rho =(float*)AllOnDev(M);
	g1=(float*)AllOnDev(M);
	MASK=(float*)AllOnDev(M);
	KSAM=(float*)AllOnDev(M);
	phi_side(float*)AllOnDev(M);
#else
//set ref for the rho equal to H_rho etc.
	if (n_pos>0) {
		Px=H_Px;
		Py=H_Py;
		Pz=H_Pz;
	}	
	MASK=H_MASK;
	KSAM=H_KSAM;  
	rho = H_rho;
	G1 = H_G1;
	phi_side = new float[M]; 
#endif

}
