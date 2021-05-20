#include <iostream> 
#include <string> 
#include "lattice.h" 
#include "LGrad3.h" 

LGrad3::LGrad3(vector<Input*> In_,string name_): Lattice(In_,name_) {}

LGrad3::~LGrad3() {
if (debug) cout <<"LGrad3 destructor " << endl;
}

void LGrad3:: ComputeLambdas() {
}

bool LGrad3::PutM() {
if (debug) cout << "PutM in LGrad3 " << endl;
	bool success=true;
	volume = MX*MY*MZ;
	JX=(MZ+2*fjc)*(MY+2*fjc); JY=MZ+2*fjc; JZ=1; M = (MX+2*fjc)*(MY+2*fjc)*(MZ+2*fjc);

	Accesible_volume=volume;
	return success;
}

void LGrad3::TimesL(Real* X){
if (debug) cout << "TimesL in LGrad3 " << endl;
}

void LGrad3::DivL(Real* X){
if (debug) cout << "DivL in LGrad3 " << endl;
}

Real LGrad3:: Moment(Real* X,Real Xb, int n) {
if (debug) cout << "Moment in LGrad3 " << endl;
	Real Result=0;
	cout <<"Moment analysis not implemented in LGrad3 " << endl;
	return Result/fjc;
}

Real LGrad3::WeightedSum(Real* X){
if (debug) cout << "weighted sum in LGrad3 " << endl;
	Real sum{0};
	remove_bounds(X);
	Sum(sum,X,M);
	return sum;
}

void LGrad3::vtk(string filename, Real* X, string id,bool writebounds) {
if (debug) cout << "vtk in LGrad3 " << endl;
	FILE *fp;
	int i;
	fp = fopen(filename.c_str(),"w+");
	fprintf(fp,"# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\n",MZ,MY,MX);

	if (writebounds) {
		fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",(MX+2*fjc)*(MY+2*fjc)*(MZ+2*fjc));
	} else {
		fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",MX*MY*MZ);
	}
	fprintf(fp,"SCALARS %s double\nLOOKUP_TABLE default\n",id.c_str());
	if (writebounds) for(i=0; i<M; i++) fprintf(fp,"%f\n",X[i]);
	else {
		for (int x=1; x<MX+1; x++)
		for (int y=1; y<MY+1; y++)
		for (int z=1; z<MZ+1; z++)
		fprintf(fp,"%e\n",X[P(x,y,z)]);
	}
	fclose(fp);
}

void LGrad3::PutProfiles(FILE* pf,vector<Real*> X,bool writebounds){
if (debug) cout <<"PutProfiles in LGrad3 " << endl;
	int x,y,z,i;
	int length=X.size();
	int a;
	if (writebounds) a=fjc; else a = 0;
	for (x=1-a; x<MX+1+a; x++)
	for (y=1-a; y<MY+1+a; y++)
	for (z=1-a; z<MZ+1+a; z++) {
		fprintf(pf,"%e\t%e\t%e\t",1.0*x*fjc-1/(2.0*fjc),1.0*y*fjc-1/(2.0*fjc),1.0*z*fjc-1/(2.0*fjc));
		for (i=0; i<length; i++) fprintf(pf,"%.20g\t",X[i][x*JX+y*JY+fjc-1+z]);
		fprintf(pf,"\n");
	}
}


void LGrad3::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
if (debug) cout <<" Side in LGrad3 " << endl;
	if (ignore_sites) {
		Cp(X_side,X,M); return;
	}
	Zero(X_side,M);set_bounds(X);


	Add(X_side+JX,X,M-JX);
	Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY);
	Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);
	Add(X_side,X+1, M-1);
	if (stencil_full) {
		if (lattice_type == "simple_cubic") {
			Norm(X_side,4.0,M);
		} else {
			Norm(X_side,2.0,M);
		}
		Add(X_side+JX+JY, X,         M-JX-JY);
		Add(X_side,         X+JX+JY, M-JX-JY);
		Add(X_side+JY,     X+JX,     M-JY-JX);
		Add(X_side+JX,      X+JY,    M-JY-JX);
		Add(X_side+JX+1,   X,        M-JX-1);
		Add(X_side,         X+JX+1,  M-JX-1);
		Add(X_side+JX,     X+1,      M-JX);
		Add(X_side+1,       X+JX,    M-JX);
		Add(X_side+JY+1,   X,        M-JY-1);
		Add(X_side,         X+JY+1,  M-JX-1);
		Add(X_side+JY,     X+1,      M-JY);
		Add(X_side+1,       X+JY,    M-JY);
		if (lattice_type == "simple_cubic") {
			Norm(X_side,4.0,M);
		} else {
			Norm(X_side,2.0,M);
		}
		Add(X_side+JX+JY+1,  X,	    M-JX-JY-1);
		Add(X_side,          X+JX+JY+1, M-JX-JY-1);
		Add(X_side+JX+JY,    X+1,       M-JX-JY-1);
		Add(X_side+1,        X+JX+JY,   M-JX-JY-1);
		Add(X_side+JX+1,     X+JY,      M-JX-JY-1);
		Add(X_side+JY,       X+JX+1,    M-JX-JY-1);
		Add(X_side+JY+1,     X+JX,      M-JX-JY-1);
		Add(X_side+JX,       X+JY+1,    M-JX-JY-1);
		if (lattice_type == "simple_cubic") {
			Norm(X_side,1.0/152.0,M);
		} else {
			Norm(X_side,1.0/56.0,M);
		}
	} else Norm(X_side,1.0/6.0,M);

}

void LGrad3::propagate(Real *G, Real *G1, int s_from, int s_to,int M) { //this procedure should function on simple cubic lattice.
if (debug) cout <<" propagate in LGrad3 " << endl;
	Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	int JX_=JX, JY_=JY; 
	int k=sub_box_on;

	Zero(gs,M); set_bounds(gs_1);

	if (k>0) {
		JX_=jx[k];
		JY_=jy[k];
	}

	if (stencil_full) {
		Add(gs+JX_,gs_1,M-JX_);
		Add(gs,gs_1+JX_,M-JX_);
		Add(gs+JY_,gs_1,M-JY_);
		Add(gs,gs_1+JY_,M-JY_);
		Add(gs+1,gs_1,M-1);
		Add(gs,gs_1+1, M-1);
		if (lattice_type == "simple_cubic") {
			Norm(gs,4.0,M);
		} else {
			Norm(gs,2.0,M);
		}
		Add(gs+JX_+JY_,gs_1,M-JX_-JY_);
		Add(gs,gs_1+JX_+JY_,M-JX_-JY_);
		Add(gs+JY_,gs_1+JX,M-JY_-JX_);
		Add(gs+JX,gs_1+JY_,M-JY_-JX_);
		Add(gs+JX_+1,gs_1,M-JX_-1);
		Add(gs,gs_1+JX_+1,M-JX_-1);
		Add(gs+JX_,gs_1+1,M-JX_);
		Add(gs+1,gs_1+JX_,M-JX_);
		Add(gs+JY_+1,gs_1,M-JY_-1);
		Add(gs,gs_1+JY_+1,M-JX_-1);
		Add(gs+JY_,gs_1+1,M-JY_);
		Add(gs+1,gs_1+JY_,M-JY_);
		if (lattice_type == "simple_cubic") {
			Norm(gs,4.0,M);
		} else {
			Norm(gs,2.0,M);
		}
		Add(gs+JX_+JY_+1,gs_1,M-JX_-JY_-1);
		Add(gs,gs_1+JX_+JY_+1,M-JX_-JY_-1);
		Add(gs+JX_+JY_,gs_1+1,M-JX_-JY_-1);
		Add(gs+1,gs_1+JX_+JY_,M-JX_-JY_-1);
		Add(gs+JX_+1,gs_1+JY_,M-JX_-JY_-1);
		Add(gs+JY_,gs_1+JX_+1,M-JX_-JY_-1);
		Add(gs+JY_+1,gs_1+JX_,M-JX_-JY_-1);
		Add(gs+JX_,gs_1+JY_+1,M-JX_-JY_-1);
		if (lattice_type == "simple_cubic") {
			Norm(gs,1.0/152.0,M);
		} else {
			Norm(gs,1.0/56.0,M);
		}
		Times(gs,gs,G1,M);
	} else {
#ifdef CUDA
		Propagate_gs_locality(gs, gs_1, G1, JX, JY, JZ, M);
#else
		Add(gs+JX_,gs_1,M-JX_);
		Add(gs,gs_1+JX_,M-JX_);
		Add(gs+JY_,gs_1,M-JY_);
		Add(gs,gs_1+JY_,M-JY_);
		Add(gs+1,gs_1,M-1);
		Add(gs,gs_1+1, M-1);
		Norm(gs,1.0/6.0,M);
		Times(gs,gs,G1,M);
#endif
	}
}


bool LGrad3::ReadRange(int* r, int* H_p, int &n_pos, bool &block, string range, int var_pos, string seg_name, string range_type) {
if (debug) cout <<"ReadRange in LGrad3 " << endl;
	bool success=true;
	vector<string>set;
	vector<string>coor;
	vector<string>xyz;
	string diggit;
	bool recognize_keyword;
	//int a; if (range_type=="frozen_range") a=1; else a=0;
	int a=0;
	In[0]->split(range,';',set);
	if (set.size()==2) {
		coor.clear();
		block=true; In[0]->split(set[0],',',coor);

		if (coor.size()!=3) {cout << "In mon " + 	seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes do not come in set of three: 'x,y,z'" << endl; success=false;}
		else {
			diggit=coor[0].substr(0,1);
			if (In[0]->IsDigit(diggit)) r[0]=In[0]->Get_int(coor[0],0); else {
				recognize_keyword=false;
				if (coor[0]=="var_pos") {recognize_keyword=true; r[0]=var_pos;}
				if (coor[0]=="firstlayer") {recognize_keyword=true; r[0] = 1;}
						//if (coor[0]=="lowerbound") {recognize_keyword=true; r[0] = 0;}
						//if (coor[0]=="upperbound") {recognize_keyword=true; r[0] = MX+1;}
				if (coor[0]=="lastlayer")  {recognize_keyword=true; r[0] = MX;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
					success=false;
				}
			}
			if (r[0] < 1-a || r[0] > MX+a) {cout << "In mon " + seg_name + ", for 'pos 1', the x-coordinate in '" + range_type + "' is out of bounds: "<< 1-a <<" .."<< MX+a << endl; success =false;}
			diggit=coor[1].substr(0,1);
			if (In[0]->IsDigit(diggit)) r[1]=In[0]->Get_int(coor[1],0); else {
				recognize_keyword=false;
				if (coor[1]=="var_pos") {recognize_keyword=true; r[1]=var_pos;}
				if (coor[1]=="firstlayer") {recognize_keyword=true; r[1] = 1;}
						//if (coor[1]=="lowerbound") {recognize_keyword=true; r[1] = 0;}
						//if (coor[1]=="upperbound") {recognize_keyword=true; r[1] = MY+1;}
				if (coor[1]=="lastlayer")  {recognize_keyword=true; r[1] = MY;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
					success=false;
				}
			}
			if (r[1] < 1-a || r[1] > MY+a) {cout << "In mon " + seg_name+ ", for 'pos 1', the y-coordinate in '" + range_type + "' is out of bounds: "<< 1-a <<" .." << MY+a << endl; success =false;}
			diggit=coor[2].substr(0,1);
			if (In[0]->IsDigit(diggit)) r[2]=In[0]->Get_int(coor[2],0); else {
				recognize_keyword=false;
				if (coor[2]=="var_pos") {recognize_keyword=true; r[2]=var_pos;}
				if (coor[2]=="firstlayer") {recognize_keyword=true; r[2] = 1;}
						//if (coor[2]=="lowerbound") {recognize_keyword=true; r[2] = 0;}
						//if (coor[2]=="upperbound") {recognize_keyword=true; r[2] = MZ+1;}
				if (coor[2]=="lastlayer")  {recognize_keyword=true; r[2] = MZ;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer'" << endl;
					success=false;
				}
			}

			if (r[2] < 1-a || r[2] > MZ+a) {cout << "In mon " + seg_name+ ", for 'pos 1', the z-coordinate in '" + range_type + "' is out of bounds: "<< 1-a <<" .." << MZ+a << endl; success =false;}
		}
		coor.clear(); In[0]->split(set[1],',',coor);

		if (coor.size()!=3) {cout << "In mon " + seg_name+ ", for 'pos 2', in '" + range_type + "', the coordinates do not come in set of three: 'x,y,z'" << endl; success=false;}
		else {
			diggit=coor[0].substr(0,1);
			if (In[0]->IsDigit(diggit)) r[3]=In[0]->Get_int(coor[0],0); else {
				recognize_keyword=false;
				if (coor[0]=="var_pos") {recognize_keyword=true; r[3]=var_pos;}
				if (coor[0]=="firstlayer") {recognize_keyword=true; r[3] = 1;}
						//if (coor[0]=="lowerbound") {recognize_keyword=true; r[3] = 0;}
						//if (coor[0]=="upperbound") {recognize_keyword=true; r[3] = MX+1;}
				if (coor[0]=="lastlayer")  {recognize_keyword=true; r[3] = MX;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
					success=false;
				}
			}
			if (r[3] < 1-a || r[3] > MX+a) {cout << "In mon " + seg_name+ ", for 'pos 2', the x-coordinate in '" + range_type + "' is out of bounds; "<< 1-a <<" .."<< MX+a << endl; success =false;}
			diggit=coor[1].substr(0,1);
			if (In[0]->IsDigit(diggit)) r[4]=In[0]->Get_int(coor[1],0); else {
				recognize_keyword=false;
				if (coor[1]=="var_pos") {recognize_keyword=true; r[4]=var_pos;}
				if (coor[1]=="firstlayer") {recognize_keyword=true; r[4] = 1;}
						//if (coor[1]=="lowerbound") {recognize_keyword=true; r[4] = 0;}
						//if (coor[1]=="upperbound") {recognize_keyword=true; r[4] = MY+1;}
				if (coor[1]=="lastlayer")  {recognize_keyword=true; r[4] = MY;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
					success=false;
				}
			}
			if (r[4] < 1-a || r[4] > MY+a) {cout << "In mon " + seg_name+ ", for 'pos 2', the y-coordinate in '" + range_type + "' is out of bounds; "<< 1-a <<" .." << MY+a << endl; success =false;}
			diggit=coor[2].substr(0,1);
			if (In[0]->IsDigit(diggit)) r[5]=In[0]->Get_int(coor[2],0); else {
				recognize_keyword=false;
				if (coor[2]=="var_pos") {recognize_keyword=true; r[5]=var_pos;}
				if (coor[2]=="firstlayer") {recognize_keyword=true; r[5] = 1;}
						//if (coor[2]=="lowerbound") {recognize_keyword=true; r[5] = 0;}
						//if (coor[2]=="upperbound") {recognize_keyword=true; r[5] = MZ+1;}
				if (coor[2]=="lastlayer")  {recognize_keyword=true; r[5] = MZ;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
					success=false;
				}
			}
			if (r[5] < 1-a || r[5] > MZ+a) {cout << "In mon " + seg_name+ ", for 'pos 2', the z-coordinate in '" + range_type + "' is out of bounds; "<< 1-a <<" .." << MZ+a << endl; success =false;}
			if (r[0] > r[3]) {cout << "In mon " + seg_name+ ", for 'pos 1', the x-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
			if (r[1] > r[4]) {cout << "In mon " + seg_name+ ", for 'pos 1', the y-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
			if (r[2] > r[5]) {cout << "In mon " + seg_name+ ", for 'pos 1', the z-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
	}
	} else {
		string s;
		In[0]->split(set[0],')',coor);
		s=coor[0].substr(0,1);
		if (s!="(") { //now expect one of the keywords
			block=true;
			cout << "In mon " + seg_name + " and  range_type " + range_type + ", the info was not recognised because when 'gradients>1' the lonely keywords 'firstlayer' 'lastlayers' do not work." << endl;
			success=false;
		} else {
			int px{0},py{0},pz{0};
			string s;

			if (coor.size()==0)
				block=false;
			else {
				for (size_t i = 0 ; i < coor.size() ; ++i) {
					s=coor[i].substr(1,coor[i].size()-1);
					In[0]->split(s,',',xyz);
					if (xyz.size()!=3) {
						cout << "In mon " + seg_name+ " pinned_range  the expected 'triple coordinate' -with brackets- structure '(x,y,z)' was not found. " << endl;  success = false;
					} else {
						px=In[0]->Get_int(xyz[0],0);
						if (px < 1 || px > MX) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the x-coordinate in pinned_range out of bounds: 1.." << MX << endl; success =false;}
						py=In[0]->Get_int(xyz[1],0);
						if (py < 1 || py > MY) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 1.." << MY << endl; success =false;}
						pz=In[0]->Get_int(xyz[2],0);
						if (pz < 1 || pz > MZ) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 1.." << MZ << endl; success =false;}
						H_p[i]=px*JX+py*JY+fjc-1+pz;
					}
				}
			}

		}
	}
	return success;
}

bool LGrad3::ReadRangeFile(string filename,int* H_p, int &n_pos, string seg_name, string range_type) {
if (debug) cout <<"ReadRangeFile in LGrad3 " << endl;
	bool success=true;
	string content;
	vector<string> lines;
	vector<string> sub;
	vector<string> xyz;
	string Infilename=In[0]->name;
	In[0]->split(Infilename,'.',sub);

	int length;
	int length_xyz;
	int px,py,pz,p_i,x,y,z;
	int i=0;
	if (!In[0]->ReadFile(sub[0].append(".").append(filename),content)) {
		success=false;
		return success;
	}

	In[0]->split(content,'#',lines);
	length = lines.size();
	if (length == MX*MY*MZ) { //expect to read 'mask file';
		i=0;
		if (n_pos==0) {
			while (i<length){
				if (In[0]->Get_int(lines[i],0)==1) n_pos++;
				i++;
			};
			if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any unities." << endl;}
		} else {
			i=0; p_i=0;
			for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++) for (z=1; z<MZ+1; z++) {
				if (In[0]->Get_int(lines[i],0)==1) {H_p[p_i]=x*JX+y*JY+fjc-1+z; p_i++;}
				i++;
			}
		}
	} else { //expect to read x,y,z
		px=0,py=0,pz=0; i=0;
		if (n_pos==0) n_pos=length;
		else {
			while (i<length) {
				xyz.clear();
				In[0]->split(lines[i],',',xyz);
				length_xyz=xyz.size();
				if (length_xyz!=3) {
					cout << "In mon " + seg_name + " " +range_type+"_filename  the expected 'triple coordinate' structure 'x,y,z' was not found. " << endl;  success = false;
				} else {
					px=In[0]->Get_int(xyz[0],0);
					if (px < 1 || px > MX) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the x-coordinate in "+range_type+"_filename out of bounds: 1.." << MX << endl; success =false;}
					py=In[0]->Get_int(xyz[1],0);
					if (py < 1 || py > MY) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the y-coordinate in "+range_type+"_filename out of bounds: 1.." << MY << endl; success =false;}
					pz=In[0]->Get_int(xyz[2],0);
					if (pz < 1 || pz > MZ) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the y-coordinate in "+range_type+"_filename out of bounds: 1.." << MZ << endl; success =false;}
				}
				H_p[i]=px*JX+py*JY+fjc-1+pz;
				i++;
			}
		}
	}
	return success;
}

bool LGrad3::FillMask(int* Mask, vector<int>px, vector<int>py, vector<int>pz, string filename) {
	bool success=true;
	bool readfile=false;
	int length=0;
	int length_px = px.size();

	vector<string> lines;
	int p;
	if (px.size()==0) {
		readfile=true;
		string content;
		success=In[0]->ReadFile(filename,content);
		if (success) {
			In[0]->split(content,'#',lines);
			length = lines.size();
		}
	}

	if (readfile) {
		if (MX*MY*MZ!=length) {
			success=false; cout <<"inputfile for filling delta_range has not expected length in x,y,z-directions" << endl;
		} else {
			for (int x=1; x<MX+1; x++)
			for (int y=1; y<MY+1; y++)
			for (int z=1; z<MZ+1; z++) Mask[x*JX + y*JY + z]=In[0]->Get_int(lines[x*JX + y*JY + z],-1);
		}
	} else  {
		for (int i=0; i<length_px; i++) {
			p=px[i]; if (p<1 || p>MX) {success=false; cout<<" x-value in delta_range out of bounds; " << endl; }
			p=py[i]; if (p<1 || p>MY) {success=false; cout<<" y-value in delta_range out of bounds; " << endl; }
			p=pz[i]; if (p<1 || p>MZ) {success=false; cout<<" z-value in delta_range out of bounds; " << endl; }
			if (success) Mask[px[i]*JX + py[i]*JY + fjc-1+ pz[i]]=1;
		}
	}
	for (int i=0; i<M; i++) if (!(Mask[i]==0 || Mask[i]==1)) {success =false; cout <<"Delta_range does not contain '0' or '1' values. Check delta_inputfile values"<<endl; }
	return success;
}

bool LGrad3::CreateMASK(int* H_MASK, int* r, int* H_P, int n_pos, bool block) {
if (debug) cout <<"CreateMask for LGrad3 " + name << endl;
	bool success=true;
	H_Zero(H_MASK,M);
	if (block) {
		for (int x=r[0]; x<r[3]+1; x++)
		for (int y=r[1]; y<r[4]+1; y++)
		for (int z=r[2]; z<r[5]+1; z++)
			H_MASK[x*JX+y*JY+fjc-1+z]=1;
	} else {
		for (int i = 0; i<n_pos; i++) H_MASK[H_P[i]]=1;
	}
	return success;
}


Real LGrad3::ComputeTheta(Real* phi) {
	Real result=0; remove_bounds(phi);
	Dot(result,phi,L,M);
	return result;
}

void LGrad3::UpdateEE(Real* EE, Real* psi, Real* E) {
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_bounds(psi);

	Zero(EE,M);
	AddGradSquare(EE+1,psi,psi+1,psi+2,M-2);
	AddGradSquare(EE+JX,psi,psi+JX,psi+2*JX,M-2*JX);
	AddGradSquare(EE+JY,psi,psi+JY,psi+2*JY,M-2*JY);
	Norm(EE,pf,M);

/* In lattice refinement. this is possibly not working...old version restored.
			YisAminB(EE+1,psi,psi+2,M-2);
			Times(EE,EE,EE,M);
			Norm(EE,pf/4.0,M);
			YisAminB(E+JX,psi,psi+2*JX,M-2*JX);
			Times(E,E,E,M);
			YplusisCtimesX(EE,E,pf/4.0,M);
			YisAminB(E+JY,psi,psi+2*JY,M-2*JY);
			Times(E,E,E,M);
			YplusisCtimesX(EE,E,pf/4.0,M);
*/
}


void LGrad3::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, int* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	int x,y;
#ifndef CUDA
	int z;
#endif

#ifndef CUDA
	Real epsZplus, epsZmin;
#endif
	Real epsXplus, epsXmin, epsYplus,epsYmin;
	set_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);

   if (!fixedPsi0) {

#ifdef CUDA
	C *=6;
	Cp(X,psi,M);
	UpPsi(g+JX+JY+1,psi+JX+JY+1,X+JX+JY+1,eps+JX+JY+1,JX,JY,C,Mask+JX+JY+1,M-2*(JX+JY+1));
	if (fjc==2) cout << "in GPU FJC-choices > 3 not implemented yet " << endl;
#else

	for (x=1; x<MX+1; x++) {
		for (y=1; y<MY+1; y++) {
			epsZplus=eps[x*JX+y*JY]+eps[x*JX+y*JY+1];
			for (z=1; z<MZ+1; z++) {
				epsZmin=epsZplus;
				epsZplus=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z+1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
					epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
				}
			}
		}
	}
	for (x=MX; x>0; x--) {
		for (y=MY; y>0; y--) {
			epsZmin=eps[x*JX+y*JY+MZ+1]+eps[x*JX+y*JY+MZ];
			for (z=MZ; z>0; z--) {
				epsZplus=epsZmin;
				epsZmin=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z-1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
					epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
					g[x*JX+y*JY+z]-=psi[x*JX+y*JY+z];
				}
			}
		}
	}



/*
	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++) {
		X[x*JX+y*JY+z]=(psi[(x-1)*JX+y*JY+z]+psi[(x+1)*JX+y*JY+z]
				 +psi[x*JX+(y-1)*JY+z]+psi[x*JX+(y+1)*JY+z]
			        +psi[x*JX+y*JY+z-1]  +psi[x*JX+y*JY+z+1])/6.0
			        +q[x*JX+y*JY+z]*C/eps[x*JX+y*JY+z]/3.0;
	}

	if (grad_epsilon) for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)  for (z=fjc; z<MZ+fjc; z++){
		X[x*JX+y*JY+z]+=0.25*(eps[(x+1)*JX+y*JY+z]-eps[(x-1)*JX+y*JY+z])*(psi[(x+1)*JX+y*JY+z]-psi[(x-1)*JX+y*JY+z]+
                              eps[x*JX+(y+1)*JY+z]-eps[x*JX+(y-1)*JY+z])*(psi[x*JX+(y+1)*JY+z]-psi[x*JX+(y-1)*JY+z]+
	                       eps[x*JX+y*JY+z-1]  -eps[x*JX+y*JY+z+1])  *(psi[x*JX+y*JY+z-1]  -psi[x*JX+y*JY+z+1])
				  /eps[x*JX+y*JY+z];
	}
	Cp(psi,X,M);
	YisAminB(g,g,psi,M);
*/

#endif
   } else { //fixedPsi0 is true

#ifdef CUDA
	C *=6;
	Cp(X,psi,M);
	UpPsi(g+JX+JY+1,psi+JX+JY+1,X+JX+JY+1,eps+JX+JY+1,JX,JY,C,Mask+JX+JY+1,M-2*(JX+JY+1));
#else
	for (x=1; x<MX+1; x++) {
		for (y=1; y<MY+1; y++) {
			epsZplus=eps[x*JX+y*JY]+eps[x*JX+y*JY+1];
			for (z=1; z<MZ+1; z++) {
				epsZmin=epsZplus;
				epsZplus=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z+1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0)
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					    epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
				           epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					    C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
			}
		}
	}
	for (x=MX; x>0; x--) {
		for (y=MY; y>0; y--) {
			epsZmin=eps[x*JX+y*JY+MZ+1]+eps[x*JX+y*JY+MZ];
			for (z=MZ; z>0; z--) {
				epsZplus=epsZmin;
				epsZmin=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z-1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==0) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
					epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]+
					C*q[x*JX+y*JY+z])/(epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus);
					g[x*JX+y*JY+z]-=psi[x*JX+y*JY+z];
				}
			}
		}
	}


/*  in lattice refinement this part is appartently not working. previous lines is for previous version...

	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++) {
		if (Mask[x*JX+y*JY+z] == 0)
		X[x*JX+y*JY+z]=(psi[(x-1)*JX+y*JY+z]+psi[(x+1)*JX+y*JY+z]
				 +psi[x*JX+(y-1)*JY+z]+psi[x*JX+(y+1)*JY+z]
			        +psi[x*JX+y*JY+z-1]  +psi[x*JX+y*JY+z+1])/6.0
			          +q[x*JX+y*JY+z]*C/eps[x*JX+y*JY+z]/3.0;
	}

	if (grad_epsilon) for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)  for (z=fjc; z<MZ+fjc; z++){
		if (Mask[x*JX+y*JY+z] == 0)
		X[x*JX+y*JY+z]+=0.25*(eps[(x+1)*JX+y*JY+z]-eps[(x-1)*JX+y*JY+z])*(psi[(x+1)*JX+y*JY+z]-psi[(x-1)*JX+y*JY+z]+
                              eps[x*JX+(y+1)*JY+z]-eps[x*JX+(y-1)*JY+z])*(psi[x*JX+(y+1)*JY+z]-psi[x*JX+(y-1)*JY+z]+
			         eps[x*JX+y*JY+z-1]  -eps[x*JX+y*JY+z+1])  *(psi[x*JX+y*JY+z-1]  -psi[x*JX+y*JY+z+1])
				   /eps[x*JX+y*JY+z];
	}
	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)  for (z=fjc; z<MZ+fjc; z++)
	if (Mask[x*JX+y*JY+z] == 0) {
		psi[x*JX+y*JY+z]=X[x*JX+y*JY+z];
		g[x*JX+y*JY+z]-=psi[x*JX+y*JY+z];
	}
*/

#endif
   } 
}


void LGrad3::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, int* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
	int x,y;
	#ifndef CUDA
	int z;
	#endif
	Real epsXplus,epsXmin,epsYplus,epsYmin,epsZplus,epsZmin;

	Real C = -e*e/(eps0*k_BT*bond_length);
#ifdef CUDA
	C *=6;
	Cp(X,psi,M);
	UpQ(g+JX+JY+1, q+JX+JY+1, psi+JX+JY+1, eps+JX+JY+1, JX, JY, C, Mask+JX+JY+1, M-2*(JX+JY+1));
#else

	for (x=1; x<MX; x++) {
		for (y=1; y<MY; y++) {
			epsZplus=eps[x*JX+y*JY]+eps[x*JX+y*JY+1];
			for (z=1; z<MZ; z++) {
				epsZmin=epsZplus;
				epsZplus=eps[x*JX+y*JY+z]+eps[x*JX+y*JY+z+1];
				epsYmin= eps[x*JX+y*JY+z]+eps[x*JX+(y-1)*JY+z];
				epsYplus=eps[x*JX+y*JY+z]+eps[x*JX+(y+1)*JY+z];
				epsXmin = eps[x*JX+y*JY+z]+eps[(x-1)*JX+y*JY+z];
				epsXplus= eps[x*JX+y*JY+z]+eps[(x+1)*JX+y*JY+z];
				if (Mask[x*JX+y*JY+z]==1) {
					psi[x*JX+y*JY+z]= (epsXmin*psi[(x-1)*JX+y*JY+z]+epsXplus*psi[(x+1)*JX+y*JY+z]+
					                   epsYmin*psi[x*JX+(y-1)*JY+z]+epsYplus*psi[x*JX+(y+1)*JY+z]+
							     epsZmin*psi[x*JX+y*JY+z-1]+epsZplus*psi[x*JX+y*JY+z+1]-
							     (epsXmin+epsXplus+epsYmin+epsYplus+epsZmin+epsZplus)*psi[x*JX+y*JY+z])/C;
					g[x*JX+y*JY+z]=-q[x*JX+y*JY+z];
				}
			}
		}
	}


/* in lattice refinement charge in 3d not working trying to restore....

	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++) {
		if (Mask[x*JX+y*JY+z] == 1)
		q[x*JX+y*JY+z]=-0.5*(psi[(x-1)*JX+y*JY+z]+psi[(x+1)*JX+y*JY+z]
				 +psi[x*JX+(y-1)*JY+z]+psi[x*JX+(y+1)*JY+z]
			        +psi[x*JX+y*JY+z-1]  +psi[x*JX+y*JY+z+1]
			       -6.0*q[x*JX+y*JY+z])*fjc*fjc*eps[x*JX+y*JY+z]/C;
	}

	if (grad_epsilon) for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)  for (z=fjc; z<MZ+fjc; z++){
		if (Mask[x*JX+y*JY+z] == 1)
		q[x*JX+y*JY+z]-=0.25*(eps[(x+1)*JX+y*JY+z]-eps[(x-1)*JX+y*JY+z])*(psi[(x+1)*JX+y*JY+z]-psi[(x-1)*JX+y*JY+z]+
                              eps[x*JX+(y+1)*JY+z]-eps[x*JX+(y-1)*JY+z])*(psi[x*JX+(y+1)*JY+z]-psi[x*JX+(y-1)*JY+z]+
				  eps[x*JX+y*JY+z-1]  -eps[x*JX+y*JY+z+1])  *(psi[x*JX+y*JY+z-1]  -psi[x*JX+y*JY+z+1])*fjc*fjc/C;
	}
	for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)  for (z=fjc; z<MZ+fjc; z++)
	if (Mask[x*JX+y*JY+z] == 1) {
		g[x*JX+y*JY+z]=-q[x*JX+y*JY+z];
	}
*/

#endif

}

void LGrad3::remove_bounds(Real *X){
if (debug) cout <<"remove_bounds in LGrad3 " << endl;
	int x,y,z;
	int k;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			RemoveBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) RemoveBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); else {
			for (x=fjc; x<MX+fjc+1; x++) for (y=fjc; y<MY+fjc+1; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+(fjc-1)-k] = 0;
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc-k]  = 0;
			}
			for (y=fjc; y<MY+fjc+1; y++) for (z=fjc; z<MZ+fjc+1; z++)  {
				for (k=0; k<fjc; k++) X[(fjc-k-1)*JX+y*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[(MX+fjc-k)*JX+y*JY+z*JZ] = 0;
			}
			for (z=fjc; z<MZ+fjc+1; z++) for (x=fjc; x<MX+fjc+1; x++){
				for (k=0; k<fjc; k++) X[x*JX+(fjc-k-1)*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc-k)*JY+z*JZ] = 0;
			}
		}
	}
}

 
void LGrad3::set_bounds(Real* X){
if (debug) cout <<"set_bounds in LGrad3 " << endl;
	int x,y,z;
	int k=0;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			SetBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); else {
			for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+(fjc-1)-k] = X[x*JX+y*JY+B_Z1[k]];
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
			}
			for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
				for (k=0; k<fjc; k++) X[(fjc-k-1)*JX+y*JY+z*JZ] = X[B_X1[k]*JX+y*JY+z*JZ];
				for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = X[B_XM[k]*JX+y*JY+z*JZ];
			}
			for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){
				for (k=0; k<fjc; k++) X[x*JX+(fjc-k-1)*JY+z*JZ] = X[x*JX+B_Y1[k]*JY+z*JZ];
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = X[x*JX+B_YM[k]*JY+z*JZ];
			}
		}
	}
}

void LGrad3::remove_bounds(int *X){
if (debug) cout <<"remove_bounds in LGrad3 " << endl;
	int x,y,z;
	int k;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			RemoveBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) RemoveBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); else {
			for (x=fjc; x<MX+fjc+1; x++) for (y=fjc; y<MY+fjc+1; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+(fjc-1)-k] = 0;
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc-k]  = 0;
			}
			for (y=fjc; y<MY+fjc+1; y++) for (z=fjc; z<MZ+fjc+1; z++)  {
				for (k=0; k<fjc; k++) X[(fjc-k-1)*JX+y*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[(MX+fjc-k)*JX+y*JY+z*JZ] = 0;
			}
			for (z=fjc; z<MZ+fjc+1; z++) for (x=fjc; x<MX+fjc+1; x++){
				for (k=0; k<fjc; k++) X[x*JX+(fjc-k-1)*JY+z*JZ] = 0;
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc-k)*JY+z*JZ] = 0;
			}
		}
	}
}

 
void LGrad3::set_bounds(int* X){
if (debug) cout <<"set_bounds in LGrad3 " << endl;
	int x,y,z;
	int k=0;
	if (sub_box_on!=0) {
		int k=sub_box_on;
		for (int i=0; i<n_box[k]; i++)
			SetBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
	} else {
		if (fjc==1) SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); else {
			for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++){
				for (k=0; k<fjc; k++) X[x*JX+y*JY+(fjc-1)-k] = X[x*JX+y*JY+B_Z1[k]];
				for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc+k]  = X[x*JX+y*JY+B_ZM[k]];
			}
			for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
				for (k=0; k<fjc; k++) X[(fjc-k-1)*JX+y*JY+z*JZ] = X[B_X1[k]*JX+y*JY+z*JZ];
				for (k=0; k<fjc; k++) X[(MX+fjc+k)*JX+y*JY+z*JZ] = X[B_XM[k]*JX+y*JY+z*JZ];
			}
			for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){
				for (k=0; k<fjc; k++) X[x*JX+(fjc-k-1)*JY+z*JZ] = X[x*JX+B_Y1[k]*JY+z*JZ];
				for (k=0; k<fjc; k++) X[x*JX+(MY+fjc+k)*JY+z*JZ] = X[x*JX+B_YM[k]*JY+z*JZ];
			}
		}
	}
}



