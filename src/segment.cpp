#include "segment.h"
#include <random>
#include <fstream>

Segment::Segment(vector<Input*> In_,vector<Lattice*> Lat_, string name_,int segnr,int N_seg) {
	In=In_; Lat=Lat_; name=name_; n_seg=N_seg; seg_nr=segnr; prepared = 0;
if (debug) cout <<"Segment constructor" + name << endl;
	KEYS.push_back("freedom");
	KEYS.push_back("valence");
	KEYS.push_back("epsilon");
	KEYS.push_back("e.psi0/kT");
	KEYS.push_back("pinned_range");
	KEYS.push_back("frozen_range");
	KEYS.push_back("tagged_range");
	KEYS.push_back("pinned_filename");
	KEYS.push_back("frozen_filename");
	KEYS.push_back("tagged_filename");
	KEYS.push_back("clamp_filename");
	KEYS.push_back("sub_box_size");
	KEYS.push_back("clamp_info");
	KEYS.push_back("fluctuation_potentials");
	KEYS.push_back("fluctuation_amplitude");
	KEYS.push_back("fluctuation_wavelength");
	KEYS.push_back("seed");
	KEYS.push_back("var_pos");
	KEYS.push_back("phi");
	KEYS.push_back("n");
	KEYS.push_back("size");
	KEYS.push_back("pos");

	Amplitude=0; labda=0; seed=1;
	var_pos=0;
	ns=1;
	all_segment=false;
	constraints=false;
}
Segment::~Segment() {
if (debug) cout <<"Segment destructor " + name << endl;
	DeAllocateMemory();
}

void Segment::DeAllocateMemory(void){
if (debug) cout << "In Segment, Deallocating memory " + name << endl;
if (!all_segment) return;
if (n_pos>0) cout <<"problem for n_pos " <<endl;
	if(n_pos>0) free(H_P);
	free(r);
	free(H_u);
	free(H_u_ext); //cout <<"H_u_ext dismissed" << endl;
	free(H_phi);
	free(H_MASK);
	free(H_alpha);
	free(H_phi_state);
#ifdef CUDA
	if(n_pos>0) cudaFree(P);
	cudaFree(u);
	cudaFree(u_ext);
	cudaFree(phi);
	cudaFree(phi_state);
	cudaFree(G1);
	cudaFree(MASK);
	cudaFree(phi_side);
#else
	free(G1);
	free(phi_side);
#endif
	all_segment=false;
}

void Segment::AllocateMemory() {
if (debug) cout <<"Allocate Memory in Segment " + name << endl;
	DeAllocateMemory();
	int M=Lat[0]->M;
	ns=state_name.size(); if (ns==0) ns=1;
	r=(int*) malloc(6*sizeof(int)); std::fill(r,r+6,0);
	H_u = (Real*) malloc(M*ns*sizeof(Real));
	H_u_ext = (Real*) malloc(M*sizeof(Real));
	H_phi = (Real*) malloc(M*sizeof(Real));
	H_MASK = (int*) malloc(M*sizeof(int));
	H_alpha=(Real*) malloc(M*ns*sizeof(Real));
	H_phi_state = (Real*) malloc(M*ns*sizeof(Real));
	H_Zero(H_u,M*ns);
	H_Zero(H_phi,M);
	H_Zero(H_MASK,M);
#ifdef CUDA

	//if (n_pos>0) Px=(int*)AllIntOnDev(n_pos);
	G1=(Real*)AllOnDev(M); Zero(G1,M);
	u=(Real*)AllOnDev(M*ns); Zero(u,M*ns);
	phi_state=(Real*)AllOnDev(M*ns); Zero(phi_state,M*ns);
	MASK=(int*)AllIntOnDev(M); Zero(MASK,M);
	phi=(Real*)AllOnDev(M); Zero(phi,M);
	u_ext=(Real*)AllOnDev(M); Zero(u_ext,M);
	phi_side=(Real*)AllOnDev(ns*M); Zero(phi_side,ns*M);
	alpha=(Real*)AllOnDev(ns*M); Zero(alpha,ns*M);
#else
	if (n_pos>0) P=H_P;
	MASK=H_MASK;
	phi =H_phi;
	u = H_u;
	u_ext=H_u_ext; Zero(u_ext,M);
	phi_state = H_phi_state;
	alpha=H_alpha;
	G1 = (Real*)malloc(M*sizeof(Real));
	phi_side = (Real*)malloc(M*ns*sizeof(Real));
	Zero(G1,M);
	Zero(phi_side,ns*M);
#endif
	bool success=true;
	bool HMaskDone=false;
	success=ParseFreedoms(HMaskDone);


	if (freedom!="free"&& !HMaskDone) {
		if (freedom=="clamp") {
			int JX=Lat[0]->JX;
			int JY=Lat[0]->JY;
			int MX=Lat[0]->MX;
			int MY=Lat[0]->MY;
			int MZ=Lat[0]->MZ;
			for (int i=0; i<n_box; i++) {
				if (bx[i]<1) {bx[i] +=MX; px1[i] +=MX; px2[i] +=MX;}
				if (by[i]<1) {by[i] +=MY; py1[i] +=MY; py2[i] +=MY;}
				if (bz[i]<1) {bz[i] +=MZ; pz1[i] +=MZ; pz2[i] +=MZ;}
				if (bx[i]<1 || bx[i]>MX) {success=false; cout <<"For cleng particle nr " << i << "the coordinate 'x' of the subbox origin is out of bounds. " << endl; }
				if (by[i]<1 || by[i]>MY) {success=false; cout <<"For cleng particle nr " << i << "the coordinate 'y' of the subbox origin is out of bounds. " << endl; }
				if (bz[i]<1 || bz[i]>MZ) {success=false; cout <<"For cleng particle nr " << i << "the coordinate 'z' of the subbox origin is out of bounds. " << endl; }
				H_MASK[((px1[i]-1)%MX+1)*JX + ((py1[i]-1)%MY+1)*JY + (pz1[i]-1)%MZ+1]=1;
				H_MASK[((px2[i]-1)%MX+1)*JX + ((py2[i]-1)%MY+1)*JY + (pz2[i]-1)%MZ+1]=1;
			}

		} else {
			r[0]*=Lat[0]->fjc; r[1]*=Lat[0]->fjc; r[2]*=Lat[0]->fjc;
			r[3]=(r[3]+1)*Lat[0]->fjc-1;r[4]=(r[4]+1)*Lat[0]->fjc-1;r[5]=(r[5]+1)*Lat[0]->fjc-1;
			Lat[0]->CreateMASK(H_MASK,r,H_P,n_pos,block);
		}

	}
	if (!success) cout <<"errors occurred.... progress uncertain...." << endl;

	all_segment=true;
}

bool Segment::ParseFreedoms(bool& HMaskDone) {
if (debug) cout <<"ParseFreedoms " << endl;
	bool success=true;
	if (freedom =="clamp" ) {
		n_box=0; mx=0;
		int m_x;
		int MX=Lat[0]->MX;
		if (GetValue("sub_box_size").size()>0) {
			m_x=In[0]->Get_int(GetValue("sub_box_size"),-1);
			if (m_x <1 || m_x > MX) {success=false; cout <<"Value of sub_box_size is out of bounds: 1 ... " << MX << endl; }
			if (mx>0) {
				if (m_x!=mx) {
					cout <<"values for sub_box_size of input " << m_x << " not consistent with value found in clamp_filename " << mx << " input file data is taken " << endl;
					mx=m_x; my=m_x; mz=m_x;
				}
			} else { mx=m_x; my=m_x; mz=m_x; }
		}
		Lat[0]->PutSub_box(mx,my,mz,n_box);
		clamp_nr = Lat[0]->m.size()-1;

		if (GetValue("clamp_filename").size()>0) {
			if (!GetClamp(GetValue("clamp_filename"))) {
				success=false; cout <<"Failed to read 'clamp_filename'. Problem terminated" << endl;
				cout <<"Example of structure of clamp_filename is the following:" << endl;
				cout <<"N : 20 : subbox0 : 15 15 15 pcb:[0. 0. 0.] " << endl;					cout <<"-1 -6 -6 " << endl;
				cout <<"1 1 1 "  << endl;
				cout <<"11 1 1 " << endl;
				cout <<" explanation: 1st line: length of chain fragment (should coinside with the composition) " << endl;
				cout <<"              followed by subboxnr, Mx, My , Mz (sizes of subbox)" << endl;
				cout <<" 	      followed by pcb info. Note that the spaces beween numbers is essential info " << endl ;
				cout <<"              2nd line: lower coordinate of subbox " << endl;
				cout <<"              3rd line: coordinates of clamp point 1 " << endl;
				cout <<"              4th line: coordinates of clamp point 2 " << endl;
				cout <<" repeat these 4 lines for every sub-box " << endl;
				cout <<" note that currently all subboxes should be equal in size " << endl;
				cout <<" note as well that the chain lengths should also be the same.(redundent information because inputfile overrules this setting" << endl;
			}
		} else {
			string s;
			if (GetValue("clamp_info").size() >0) {
				s=GetValue("clamp_info");
				vector<string> sub;
				vector<string> set;
				vector<string> coor;
				In[0]->split(s,';',sub);
				n_box=sub.size();
				px1.clear();
				py1.clear();
				pz1.clear();
				px2.clear();
				py2.clear();
				pz2.clear();
				bx.clear();
				by.clear();
				bz.clear();
				for (int i=0; i<n_box; i++) {
					set.clear();
					In[0]->split(sub[i],'(',set);
					int length = set.size();
					if (length!=3) {
						success=false; cout <<" In 'clamp_info' for segment '"+name+"', for box number " << i << " the expected format (px1,py1,pz1)(px2,py2,pz2) was not found" << endl;
					} else {
						coor.clear();
						In[0]->split(set[1],',',coor);
						if (coor.size()!=3) {
							success=false; cout <<" In 'clamp_info' for segment '"+name+"' for box number "<< i <<" the coordinates for the p1 position (px1,py1,pz1) not correct format. " << endl;
						} else {
							px1.push_back(In[0]->Get_int(coor[0],-10000));
							py1.push_back(In[0]->Get_int(coor[1],-10000));
							pz1.push_back(In[0]->Get_int(coor[2],-10000));
						}
						coor.clear();
						In[0]->split(set[2],',',coor);
						if (coor.size()!=3) {
							success=false; cout <<" In 'clamp_info' for segment '"+name+"' for box number "<< i <<" the coordinates for the box position (px2,py2,pz2) not correct format. " << endl;
						} else {
							px2.push_back(In[0]->Get_int(coor[0],-10000));
							py2.push_back(In[0]->Get_int(coor[1],-10000));
							pz2.push_back(In[0]->Get_int(coor[2],-10000));
						}
						bx.push_back((px2[i]+px1[i]-mx)/2);
						by.push_back((py2[i]+py1[i]-mx)/2);
						bz.push_back((pz2[i]+pz1[i]-mx)/2);//box is equal in size in x y and z.
					}
				}
			} else {
				success=false;
				cout<<"Segment " + name + " with 'freedom: clamp' expects input from 'clamp_filename' or 'clamp_info' " << endl;
			}
		}
	}

	if (freedom == "pinned") {
		int fjc = Lat[0]->fjc;
		int n_layers_x=(Lat[0]->MX)/fjc;
		int n_layers_y=(Lat[0]->MY)/fjc;
		int n_layers_z=(Lat[0]->MZ)/fjc;


		if (GetValue("pinned_range").size()>0) {
			s_freedom="pinned_range";
			string p_range=GetValue("pinned_range");
			vector<string>sub;
			In[0]->split(p_range,';',sub);
			int Lsub=sub.size();
			vector<string>xyz;
			for (int k=0; k<Lsub; k++) {
				xyz.clear();
				In[0]->split(sub[k],',',xyz);
				int Lxyz=xyz.size();
				for (int kk=0; kk<Lxyz; kk++){
					if (xyz[kk]=="firstlayer") {
						if (kk==0) p_range ="firstlayer_x";
						if (kk==1) p_range ="firstlayer_y";
						if (kk==2) p_range ="firstlayer_z";
					}
					if (xyz[kk]=="lastlayer") {
						if (kk==0) p_range ="lastlayer_x";
						if (kk==1) p_range ="lastlayer_y";
						if (kk==2) p_range ="lastlayer_z";
					}
				}
			}
			if (p_range=="firstlayer_x") {
				if (Lat[0]->gradients==1) {p_range = "firstlayer;firstlayer"; }
				if (Lat[0]->gradients==2) {p_range = "firstlayer,1;firstlayer,";p_range.append(to_string(n_layers_y)); }
				if (Lat[0]->gradients==3) {p_range = "firstlayer,1,1;firstlayer,"; p_range.append(to_string(n_layers_y)).append(",").append(to_string(n_layers_z)); }
			}
			if (p_range=="lastlayer_x") {
				if (Lat[0]->gradients==1) {p_range = "lastlayer;lastlayer";}
				if (Lat[0]->gradients==2) {p_range = "lastlayer,1;lastlayer,"; p_range.append(to_string(n_layers_y)); }
				if (Lat[0]->gradients==3) {p_range = "lastlayer,1,1;lastlayer,"; p_range.append(to_string(n_layers_y)).append(",").append(to_string(n_layers_z)); }
			}
			if (p_range=="firstlayer_y") {
				if (Lat[0]->gradients==2) {p_range = "1,firstlayer;";p_range.append(to_string(n_layers_x)).append(",firstlayer"); }
				if (Lat[0]->gradients==3) {p_range = "1,firstlayer,1;";p_range.append(to_string(n_layers_x)).append(",firstlayer,").append(to_string(n_layers_z)); }
			}
			if (p_range=="lastlayer_y") {
				if (Lat[0]->gradients==2) {p_range = "1,lastlayer;";p_range.append(to_string(n_layers_x)).append(",lastlayer"); }
				if (Lat[0]->gradients==3) {p_range = "1,lastlayer,1;";p_range.append(to_string(n_layers_x)).append(",lastlayer,").append(to_string(n_layers_z)); }
			}
			if (p_range=="firstlayer_z") {
				if (Lat[0]->gradients==3) {p_range = "1,1,firstlayer;";p_range.append(to_string(n_layers_x)).append(",").append(to_string(n_layers_y)).append(",firstlayer"); }
			}
			if (p_range=="lastlayer_z") {
				if (Lat[0]->gradients==3) {p_range = "1,1,lastlayer;";p_range.append(to_string(n_layers_x)).append(",").append(to_string(n_layers_y)).append(",lastlayer"); }
			}
			phibulk=0;
			if (GetValue("frozen_range").size()>0 || GetValue("tagged_range").size()>0 || GetValue("frozen_filename").size()>0 || GetValue("tag_filename").size()>0) {
			cout<< "For mon :" + name + ", you should exclusively combine freedom : pinned with pinned_range or pinned_filename" << endl;  success=false;}
			if (GetValue("pinned_range").size()>0 && GetValue("pinned_filename").size()>0) {
				cout<< "For mon " + name + ", you can not combine pinned_range with 'pinned_filename' " <<endl; success=false;
			}
			if (GetValue("pinned_range").size()==0 && GetValue("pinned_filename").size()==0) {
				cout<< "For mon " + name + ", you should provide either pinned_range or pinned_filename " <<endl; success=false;
			}
			sub.clear();
			In[0]->split(p_range,';',sub);
			p_range.clear();
			Lsub=sub.size();

			if (Lsub!=2) {
				cout <<"For mon " + name + ", the parsing of 'pinned_range' failed. Use x1,y1,z1;x2,y2,z2, x1,y1;x2,y2, or x1;x2 for 3, 2, or 1  gradient computations, respectively. x, y and z can also be  keys: 'firstlayer', 'lastlayer'" << endl;
				success=false;
				return success;
			}

			for (int k=0; k<Lsub; k++) {
				xyz.clear();
				In[0]->split(sub[k],',',xyz);
				int Lxyz=xyz.size();
				if (Lxyz<1 || Lxyz>3){
					cout <<"For mon " + name + ", the parsing of 'pinned_range' failed. Number of coordinates should be 1, 2 or 3: e.g., x1,y1,z1;x2,y2,z2, x1,y1;x2,y2, x1;x2 for 1, 2 or 3 gradients, respectively.  " << endl;
					success=false;
					return success;
				}
				for (int kk=0; kk<Lxyz; kk++) {
					if (xyz[kk]=="firstlayer") {
						p_range.append("1");
					} else if (xyz[kk]=="lastlayer") {
						if (kk==0) p_range.append(to_string(n_layers_x));
						if (kk==1) p_range.append(to_string(n_layers_y));
						if (kk==2) p_range.append(to_string(n_layers_z));
					} else if (xyz[kk]=="var_pos") {
						if (((kk==0) && (var_pos<1 || var_pos> n_layers_x)) || ((kk==1) && (var_pos<1 || var_pos> n_layers_y))  ||((kk==2) && (var_pos<1 || var_pos> n_layers_z))) {
							cout <<"In pinned_range, 'var_pos' entry is out of bounds..." << endl; success=false; return success;
						}
						p_range.append(to_string(var_pos));

					} else {
						int cor=In[0]->Get_int(xyz[kk],-1);
						if (((kk==0) && (cor <1 || cor > n_layers_x)) || ((kk==1) && (cor <1 || cor > n_layers_y))  ||((kk==2) && (cor <1 || cor > n_layers_z))) {
							cout <<" For mon " + name+ ", the 'pinned_range' is not parsed properly! Coordinates either out of bounds or keywords 'var_pos', 'firstlayer', 'lastlayer' were not found" << endl;
							success=false;
							return success;
						} else p_range.append(xyz[kk]);
					}
					if (kk<Lxyz-1) p_range.append(",");
				}
				if (k<Lsub-1) p_range.append(";");
			}

//cout <<"p_range = " << p_range << endl;

			n_pos=0;
			if (success) success=Lat[0]->ReadRange(r, H_P, n_pos, block, p_range,var_pos,name,s_freedom);
			if (n_pos>0) {
				H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
				if (success) success=Lat[0]->ReadRange(r, H_P, n_pos, block, p_range,var_pos,name,s_freedom);
			}
		}
		if (GetValue("pinned_filename").size()>0) { s_freedom="pinned";
			filename=GetValue("pinned_filename");
			n_pos=0;
			if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s_freedom);
			if (n_pos>0) {
				H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
				if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s_freedom);
			}
		}
	}

	if (freedom == "frozen") {
		int n_layers_x=(Lat[0]->MX)/Lat[0]->fjc;
		int n_layers_y=(Lat[0]->MY)/Lat[0]->fjc;
		int n_layers_z=(Lat[0]->MZ)/Lat[0]->fjc;

		frozen_at_bound=-1;
		phibulk=0;
		if (GetValue("pinned_range").size()>0 || GetValue("tagged_range").size()>0 || GetValue("pinned_filename").size()>0 || GetValue("tag_filename").size()>0) {
		        cout<< "For mon " + name + ", you should exclusively combine 'freedom : frozen' with 'frozen_range' or 'frozen_filename'" << endl;  success=false;
		}
		if (GetValue("frozen_range").size()>0 && GetValue("frozen_filename").size()>0) {
			cout<< "For mon " + name + ", you can not combine 'frozen_range' with 'frozen_filename' " <<endl; success=false;
		}
		if (GetValue("frozen_range").size()==0 && GetValue("frozen_filename").size()==0 && GetValue("n").size()==0) {
			cout<< "For mon " + name + ", you should provide either 'frozen_range' or 'frozen_filename' or specify the number of particles by  'n' " <<endl; success=false;
		}
		if (GetValue("frozen_range").size()>0) {
			s_freedom="frozen_range";
			string f_range=GetValue("frozen_range");
			vector<string>sub;
			In[0]->split(f_range,';',sub);
			int Lsub=sub.size();
			vector<string>xyz;
			for (int k=0; k<Lsub; k++) {
				xyz.clear();
				In[0]->split(sub[k],',',xyz);
				int Lxyz=xyz.size();
				for (int kk=0; kk<Lxyz; kk++){
					if (xyz[kk]=="lowerbound") {
						if (kk==0) f_range ="lowerbound_x";
						if (kk==1) f_range ="lowerbound_y";
						if (kk==2) f_range ="lowerbound_z";
					}
					if (xyz[kk]=="upperbound") {
						if (kk==0) f_range ="upperbound_x";
						if (kk==1) f_range ="upperbound_y";
						if (kk==2) f_range ="upperbound_z";
					}
				}
			}

			if (f_range=="lowerbound_x") {
				if (Lat[0]->gradients==1) {f_range = "lowerbound;lowerbound"; frozen_at_bound=0;}
				if (Lat[0]->gradients==2) {f_range = "lowerbound,1;lowerbound,"; f_range.append(to_string(n_layers_y)); frozen_at_bound=0;}
				if (Lat[0]->gradients==3) {f_range = "lowerbound,1,1;lowerbound"; f_range.append(to_string(n_layers_y)).append(",").append(to_string(n_layers_z)); frozen_at_bound=0;}
			}
			if (f_range=="upperbound_x") {
				if (Lat[0]->gradients==1) {f_range = "upperbound;upperbound";frozen_at_bound=3;}
				if (Lat[0]->gradients==2) {f_range = "upperbound,1;upperbound,"; f_range.append(to_string(n_layers_y)); frozen_at_bound=3;}
				if (Lat[0]->gradients==3) {f_range = "upperbound,1,1;upperbound"; f_range.append(to_string(n_layers_y)).append(",").append(to_string(n_layers_z)); frozen_at_bound=3;}
			}
			if (f_range=="lowerbound_y") {
				if (Lat[0]->gradients==2) {f_range = "1,lowerbound;";f_range.append(to_string(n_layers_x)).append(",lowerbound"); frozen_at_bound=1;}
				if (Lat[0]->gradients==3) {f_range = "1,lowerbound,1;";f_range.append(to_string(n_layers_x)).append(",lowerbound,").append(to_string(n_layers_z)); frozen_at_bound=1;}
			}
			if (f_range=="upperbound_y") {
				if (Lat[0]->gradients==2) {f_range = "1,upperbound;";f_range.append(to_string(n_layers_x)).append(",upperbound");  frozen_at_bound=4;}
				if (Lat[0]->gradients==3) {f_range = "1,upperbound,1;";f_range.append(to_string(n_layers_x)).append(",upperbound,").append(to_string(n_layers_z)); frozen_at_bound=4;}
			}

			if (f_range=="lowerbound_z") {
				if (Lat[0]->gradients==3) {f_range = "1,1,lowerbound;";f_range.append(to_string(n_layers_x)).append(",").append(to_string(n_layers_y)).append(",lowerbound"); frozen_at_bound=2; }
			}
			if (f_range=="upperbound_z") {
				if (Lat[0]->gradients==3) {f_range = "1,1,upperbound;";f_range.append(to_string(n_layers_x)).append(",").append(to_string(n_layers_y)).append(",upperbound"); frozen_at_bound=5;}
			}

			sub.clear();
			In[0]->split(f_range,';',sub);
			f_range.clear();

			Lsub=sub.size();
			if (Lsub!=2) {
				cout <<"For mon " + name + ", the parsing of 'frozen_range' failed. Use x1,y1,z1;x2,y2,z2 in 3 gradients, x1,y1;x2,y2 in two gradients x1;x2 for one gradient computations. x, y and z can also be  keys: 'firstlayer', 'lastlayer', 'lowerbound', or 'upperbound'." << endl;
				cout <<"Alternatively - you can try a single keyword such as 'lowerbound', 'lowerbound_x', 'lowerbound_y', 'lowerbound_z', 'upperbound', 'upperbound_x', 'upperbound_y', 'upperbound_z'" << endl;
				success=false;
				return success;
			}
			for (int k=0; k<Lsub; k++) {
				xyz.clear();
				In[0]->split(sub[k],',',xyz);
				int Lxyz=xyz.size();
				if (Lxyz<1 || Lxyz>3){
					cout <<"For mon " + name + ", the parsing of 'frozen_range' failed. Number of coordinates should be 1, 2 or 3: e.g., x1,y1,z1;x2,y2,z2, x1,y1;x2,y2, x1;x2 for 1, 2 or 3 gradients, respectively.  " << endl;
					success=false;
					return success;
				}

				for (int kk=0; kk<Lxyz; kk++) {
					if (xyz[kk]=="firstlayer")
						f_range.append("1");
					if (xyz[kk]=="lastlayer") {
						if (kk==0) f_range.append(to_string(n_layers_x));
						if (kk==1) f_range.append(to_string(n_layers_y));
						if (kk==2) f_range.append(to_string(n_layers_z));
					}

					if (xyz[kk]=="lowerbound") {
						if ( (kk==0 && Lat[0]->BC[0] !="surface") ||(kk==1 && Lat[0]->BC[1] !="surface") ||(kk==2 && Lat[0]->BC[2] !="surface") ) {
							cout<<"In lattice you need boundary condition 'surface' in combination with frozen_range containing 'lowerbound' " << endl; success=false;
							return success;
						}
						f_range.append("0");
					}
					if (xyz[kk]=="upperbound") {
						if ( (kk==0 && Lat[0]->BC[3] !="surface") ||(kk==1 && Lat[0]->BC[4] !="surface") ||(kk==2 && Lat[0]->BC[5] !="surface") ) {
							cout<<"In lattice you need boundary condition 'surface' in combination with frozen_range containing 'upperbound' " << endl; success=false;
							return success;
						}
						if (kk==0) f_range.append(to_string(n_layers_x+1));
						if (kk==1) f_range.append(to_string(n_layers_y+1));
						if (kk==2) f_range.append(to_string(n_layers_z+1));
					}
					if (xyz[kk]=="var_pos") {
						if (((kk==0) && (var_pos<0 || var_pos> n_layers_x+1)) || ((kk==1) && (var_pos<0 || var_pos> n_layers_y+1))  ||((kk==2) && (var_pos<0 || var_pos> n_layers_z+1))) {
							cout <<"In frozen_range, 'var_pos' entry is out of bounds..." << endl; success=false; return success;
						}

						f_range.append(to_string(var_pos)); //check if var_pos is within lattice-range not implemented.
					}

					if (xyz[kk]!="firstlayer" && xyz[kk]!="lastlayer" && xyz[kk]!="lowerbound" && xyz[kk]!="upperbound" && xyz[kk]!="var_pos") {
						int cor=In[0]->Get_int(xyz[kk],-1);
						if ((kk==0 && (cor <0 || cor > n_layers_x+1)) || (kk==1 && (cor <0 || cor > n_layers_y+1))  ||(kk==2 && (cor <0 || cor > n_layers_z+1))) {
							cout <<" For mon " + name+ ", the 'frozen_range' is not parsed properly! Coordinates either out of bounds or keywords  'var_pos', 'firstlayer', 'lastlayer', 'lowerbound', 'upperbound' were not found" << endl;
							success=false;
							return success;
						} else {
							if (kk==0 && cor ==0 ) {
								if (Lat[0]->BC[0] !="surface") {cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << endl; success=false; return success;}
								frozen_at_bound=0;
							}
							if (kk==0 && cor == n_layers_x+1 ) {
								if (Lat[0]->BC[3] !="surface") {cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << endl; success=false; return success;}
								frozen_at_bound=3;
							}
							if (kk==1 && cor ==0 ) {
								if (Lat[0]->BC[1] !="surface") {cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << endl; success=false; return success;}
								frozen_at_bound=1;
							}
							if (kk==1 && cor == n_layers_y+1 ) {
								if (Lat[0]->BC[4] !="surface") {cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << endl; success=false; return success;}
								frozen_at_bound=4;
							}
							if (kk==2 && cor ==0 ) {
								if (Lat[0]->BC[2] !="surface") {cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << endl; success=false; return success;}
								frozen_at_bound=2;
							}
							if (kk==2 && cor == n_layers_z+1 ) {
								if (Lat[0]->BC[5] !="surface") {cout <<" Frozen segment " + name + " put at boundary, but lattice bc is not set to 'surface'" << endl; success=false; return success;}
								frozen_at_bound=5;
							}
							f_range.append(xyz[kk]);
						}
					}
					if (kk<Lxyz-1) f_range.append(",");
				}
				if (k<Lsub-1) f_range.append(";");
			}

//cout <<"f_range = " << f_range << endl;

			n_pos=0;
			success=Lat[0]->ReadRange(r, H_P, n_pos, block, f_range,var_pos,name,s_freedom);
			if (n_pos>0) {
				H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
				success=Lat[0]->ReadRange(r, H_P, n_pos, block, f_range,var_pos,name,s_freedom);
			}
		}
		if (GetValue("frozen_filename").size()>0) { s_freedom="frozen";
			filename=GetValue("frozen_filename");
			n_pos=0;
			if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s_freedom);
			if (n_pos>0) {
				H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
				if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s_freedom);
			}
		}
		if (n_pos<1 && Lat[0]->gradients==3&& px.size()==0) {
			n=0; R=0; //px.clear(); py.clear(); pz.clear();
			bool found=false;
			if (GetValue("n").size()==0 || GetValue("pos").size()==0 || GetValue("size").size()==0) {
				success=false; cout <<"Expecting values for 'n', 'pos' and 'size' for the definition of the set of spherical particles in the system. " << endl;
			} else {
				n=In[0]->Get_int(GetValue("n"),-1); if (n<0) {success=false; cout <<" expecting positive integer for 'n'" << endl;}
				R=In[0]->Get_int(GetValue("size"),-1); if (R<0) {success =false ; cout <<" expecting positive integer for 'size' "<<endl; }
				//if (GetValue("pos")=="random") {
				//	found=true; srand(1);
				//	for (int i=0; i<n; i++) {
				//		px.push_back(1+rand()%Lat[0]->MX);
				//		py.push_back(1+rand()%Lat[0]->MY);
				//		pz.push_back(1+rand()%Lat[0]->MZ);
				//	}
				//}
				if (GetValue("pos")=="?") {
					cout << "Info: 'pos' variable should contain either the keywords 'regular' or 'random' or a sequence of n times '(x,y,z)' coordinates" << endl;
				       	success=false;
				}
				if (GetValue("pos")=="regular"||GetValue("pos")=="random") {
					found=true;
					int MX=Lat[0]->MX;
					int MY=Lat[0]->MY;
					int MZ=Lat[0]->MZ;
					int V=MX*MY*MZ;
					int a = pow(V/n,1.0/3.0);
					if (a<2*R) {
						cout << "Warning: overlap can not be avoided " << endl;
					}
					int n_placed=0;
					if (a>2*R+1) a--;

					for (int i=1; i<MX-a+2; i+=a)
					for (int j=1; j<MY-a+2; j+=a)
					for (int k=1; k<MZ-a+2; k+=a) {
						if (n_placed < n) {
							px.push_back(i); py.push_back(j); pz.push_back(k);
							n_placed++;
							//cout <<"x,y,z= " << i << "," << j << "," << k << endl;
						}

					}
					if (n_placed < n) { cout <<"Problem: could not place all n particles in volume. Placed only  "<< n_placed << " partictles " << endl; }
					//cout <<"In frozen_range regular positions for particles not implemented " <<endl;
				} else {
					cout <<"pos should be 'regular' or 'random' " << endl; success = false;

				}

				if (GetValue("pos")=="random" && success) {
					srand(1);
					int num_of_moves =100*n*n;
					int I;
					int X;
					int Y;
					int Z;
					int MX=Lat[0]->MX;
					int MY=Lat[0]->MY;
					int MZ=Lat[0]->MZ;
					for (int i=0; i<num_of_moves; i++) { //randomise...
						I=rand()%n;
						X=px[I]; Y=py[I]; Z=pz[I];
						px[I]+=rand()%3-1; if (px[I]<1) px[I]+=MX; if (px[I]>MX) px[I]-=MX;
						py[I]+=rand()%3-1; if (py[I]<1) py[I]+=MY; if (py[I]>MY) py[I]-=MY;
						pz[I]+=rand()%3-1; if (pz[I]<1) pz[I]+=MZ; if (pz[I]>MZ) pz[I]-=MZ;
						if (Overlap(I,R)) {
							px[I]=X; py[I]=Y; pz[I]=Z;
						}
					}
				}
				if (!found) {
					vector<int>open;
					vector<int>close;
					vector<string>sub;
					string t=GetValue("pos");
					if (!In[0]->EvenBrackets(t,open,close)) {cout << "Brackets in 'pos' not balanced "<<endl; success=false;};
					int opensize=open.size();
			        	if (opensize !=n ) {
							success=false; cout <<"number of positions not equal to 'n' " << endl;
						} else {
							for (int i=0; i<n; i++) {
								sub.clear();
							string tt=t.substr(open[i]+1,close[i]-open[i]-1); //cout << tt << endl;
							In[0]->split(tt,',',sub);
							if (sub.size() !=3) {
								success=false;
								cout <<"Format error: 'pos' did not contain the keywords 'regular' nor 'random'," << endl;
								cout <<"nor did 'pos'  contain expected (x,y,z) sets. Problem occurred for for particle nr " << i << "We found: " +tt << endl;
							} else {
								px.push_back(In[0]->Get_int(sub[0],-1));
								py.push_back(In[0]->Get_int(sub[1],-1));
								pz.push_back(In[0]->Get_int(sub[2],-1));
								if (px[i] <0 || px[i]>Lat[0]->MX) {success=false; cout << "pos x for particle " << i << " out of bounds or not an integer: " << px[i] << endl; }
								if (py[i] <0 || py[i]>Lat[0]->MY) {success=false; cout << "pos y for particle " << i << " out of bounds or not an integer: " << py[i] << endl; }
								if (pz[i] <0 || pz[i]>Lat[0]->MZ) {success=false; cout << "pos z for particle " << i << " out of bounds or not an integer: " << pz[i] << endl; }
								//cout <<"px,py,pz =" <<px[i]<<","<<py[i]","<<pz[i] << endl;
							}
						}
					}
				}
			}
		}
		if (px.size()>0) {
			HMaskDone=true;
			if (success) {
				if (!Lat[0]->PutMask(H_MASK,px,py,pz,R)) cout <<"overlap occurred"<<endl;
			}
		}
	}

	n_pos=-1;
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
		if (GetValue("tagged_range").size()>0) { s_freedom="tagged_range";
			string t_range=GetValue("tagged_range");
			vector<string>sub;
			In[0]->split(t_range,';',sub);
			t_range.clear();
			vector<string>xyz;
			int Lsub=sub.size();
			if (Lsub!=2) {
				cout <<"For mon " + name + ", the parsing of 'tagged_range' failed. Use x1,y1,z1;x2,y2,z2, x1,y1;x2,y2, or x1;x2 for 3, 2, or 1  gradient computations, respectively. x, y and z can also be  keys: 'firstlayer', 'lastlayer'" << endl;
				success=false;
				return success;
			}

			int n_layers_x=(Lat[0]->MX+1)/Lat[0]->fjc;
			int n_layers_y=(Lat[0]->MY+1)/Lat[0]->fjc;
			int n_layers_z=(Lat[0]->MZ+1)/Lat[0]->fjc;
			for (int k=0; k<Lsub; k++) {
				xyz.clear();
				In[0]->split(sub[k],',',xyz);
				int Lxyz=xyz.size();
				if (Lxyz<1 || Lxyz>3){
					cout <<"For mon " + name + ", the parsing of 'tagged_range' failed. Number of coordinates should be 1, 2 or 3: e.g., x1,y1,z1;x2,y2,z2, x1,y1;x2,y2, x1;x2 for 1, 2 or 3 gradients, respectively.  " << endl;
					success=false;
					return success;
				}
				for (int kk=0; kk<Lxyz; kk++) {
					if (xyz[kk]=="firstlayer") {
						t_range.append("1");
					} else if (xyz[kk]=="lastlayer") {
						if (kk==0) t_range.append(to_string(n_layers_x));
						if (kk==1) t_range.append(to_string(n_layers_y));
						if (kk==2) t_range.append(to_string(n_layers_z));
					} else {
						int cor=In[0]->Get_int(xyz[kk],-1);
						if ((kk==0 && (cor <1 || cor > n_layers_x)) || (kk==1 && (cor <1 || cor > n_layers_y))  ||(kk==2 && (cor <1 || cor > n_layers_z))) {
							cout <<" For mon " + name+ ", the 'tagged_range' is not parsed properly! Coordinates either out of bounds or keywords 'firstlayer', 'lastlayer' were not found" << endl;
							success=false;
							return success;
						} else t_range.append(xyz[kk]);
					}
					if (kk<Lxyz-1) t_range.append(",");
				}
				if (k<Lsub-1) t_range.append(";");
			}

//cout <<"t_range = " << t_range << endl;

			n_pos=0;
			if (success) success=Lat[0]->ReadRange(r, H_P, n_pos, block, t_range,var_pos,name,s_freedom);
			if (n_pos>0) {
				H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
				if (success) success=Lat[0]->ReadRange(r, H_P, n_pos, block, t_range,var_pos,name,s_freedom);
			}
		}
		if (GetValue("tagged_filename").size()>0) {
			s_freedom="tagged";
			filename=GetValue("tagged_filename");
			n_pos=0;
			if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s_freedom);
			if (n_pos>0) {
				H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
				if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s_freedom);
			}
		}
	}
	return success;
}

bool Segment::Overlap(int I, int R) {
	int X=px[I];
	int Y=py[I];
	int Z=pz[I];
	int n=px.size();
	int RR=4*R*R;
	int dx,dy,dz;
	int MX=Lat[0]->MX;
	int MY=Lat[0]->MY;
	int MZ=Lat[0]->MZ;

	for (int i=0; i<n; i++) {
		if (i != I) {
			dx=px[i]-X;
			if (dx>MX/2)  dx-=MX;
			if (dx<-MX/2) dx+=MX;
			dy=py[i]-Y;
			if (dy>MY/2)  dy-=MY;
			if (dy<-MY/2) dy+=MY;
			dz=pz[i]-Z;
			if (dz>MZ/2)  dz-=MZ;
			if (dz<-MZ/2) dz+=MZ;
			if (dx*dx+dy*dy+dz*dz<RR+4) return true;
		}
	}
	return false;
}

bool Segment::PrepareForCalculations(int* KSAM, bool first_time) {
if (debug) cout <<"PrepareForCalcualtions in Segment " +name << endl;

	int M=Lat[0]->M;
#ifdef CUDA
	if (In[0]->MesodynList.empty() or prepared == false) {
	TransferDataToDevice(H_MASK, MASK, M);
		if (In[0]->MesodynList.empty())
			TransferDataToDevice(H_u, u, M); //Wrong: This clears u for every CUDA application and messes up mesodyn
		prepared = true;
}

//}
	//TransferDataToDevice(H_Px, Px, n_pos);
	//TransferDataToDevice(H_Py, Py, n_pos);
	//TransferDataToDevice(H_Pz, Pz, n_pos);
#endif

	bool success=true;
	phibulk=0;
	if (freedom=="frozen") {
		Cp(phi,MASK,M);
	} else Zero(phi,M);

	if (freedom=="tagged" || freedom=="clamp" ) Zero(u,M); //no internal states for these segments.

	if (ns==1) {
		Lat[0]->set_bounds(u);
		Boltzmann(G1,u,M);
	} else {
		Zero(G1,M);
		for (int i=0; i<ns; i++) {
			Lat[0]->set_bounds(u+M*i);
			Boltzmann(alpha+M*i,u+M*i,M);
			Norm(alpha+M*i,state_alphabulk[i],M);
			Add(G1,alpha+M*i,M);
		}
		for (int i=0; i<ns; i++) Div(alpha+i*M,G1,M);
	}

	if (freedom=="pinned") Times(G1,G1,MASK,M);
	if (freedom=="tagged") {
		Cp(G1,MASK,M);
	}
	if (!(freedom ==" frozen" || freedom =="tagged")) Times(G1,G1,KSAM,M);
	if (GetValue("seed").size()>0) {
		seed=In[0]->Get_int(GetValue("seed"),1);
	}
	if (GetValue("fluctuation_potentials").size()>0&& first_time)
	{
		Zero(u_ext,M); srand(seed);
		int gradients=Lat[0]->gradients;
		vector<string> sub;
		string s;
		int my;
		int labda_y;
		int MX=Lat[0]->MX;
		int JX=Lat[0]->JX;
		int MY=Lat[0]->MY;
		int JY=Lat[0]->JY;
		Real shift_x,shift_y,shift_z;
		switch (gradients)
		{
			case 1:
				s = GetValue("fluctuation_potentials");
				In[0]->split(s, ',', sub);
				if (sub.size() !=1) {
					success=false; cout <<"expecting in 'mon : " + name + " : fluctuation_potentials : '  coordinate info in 1d, such as: x"<<endl;
				}
				labda=MX;
				cout <<"putting u_ext" << endl;
				for (int x=1; x<MX; x++) u_ext[x]+=Amplitude*(sin(2.0*PIE*x/labda));
				break;
			case 2:
				s = GetValue("fluctuation_potentials");
				In[0]->split(s, ',', sub);
				if (sub.size() !=2) {success=false; cout <<"expecting in 'mon : " + name + " : fluctuation_potentials : '  coordinate info in 2d, such as: x,5"<<endl; }
				my=In[0]->Get_int(sub[1],0);
				if (my<0 || my>Lat[0]->MY) {success =false; cout << "in fluctuation potentials the y-coordinate is out of bounds."<< endl; }
				labda_y=Lat[0]->MY;
				JX=Lat[0]->JX;
				for (int x=1; x<MX; x++) u_ext[x*JX+my]+=Amplitude*(sin(2.0*PIE*x/labda_y));
				break;
			case 3:
			s = GetValue("fluctuation_potentials");
			In[0]->split(s, ',', sub);
			if (sub.size()<3)
			{
				success=false;
				cout <<"expecting in 'mon : " + name + " : fluctuation_potentials : '  coordinate info in 3d, such as: x,y,5 or x,y,z"<<endl;
			}	else
			{
				if (sub[0] != "x" || sub[1] != "y" )
				{
					success=false;
					cout <<"expecting in 'mon : " + name + " : fluctuation_potentials : '  first two coordinates to be : x,y  "<<endl;
				}

				if (sub[2] == "z")
				{
					int MZ=Lat[0]->MZ;
					if (!(MZ==2 || MZ==4 || MZ==8 || MZ==16 ||MZ==32 || MZ==64 ||MZ==128 || MZ==256 || MZ==512 || MZ==1024))
					{
						success=false;
						cout << "Expecting n_layers_z to have a value 2^a with a = 1..10" << endl;
					}
					if (success)
					{
						for (int lambda=2; lambda <=MX; lambda*=2)
						{
							shift_x = rand() % lambda; //cout <<"setting shift_x: " << shift_x<< endl;
							shift_y = rand() % lambda; //cout <<"setting shift_y: " << shift_y <<endl;
							shift_z = rand() % lambda; //cout <<"setting shift_z: " << shift_z <<endl;
							for (int x=0; x<MX; x++) for (int y=0; y<MY; y++) for (int z=0; z<MZ; z++) u_ext[x*JX+y*JY+z]+=Amplitude*(sin(2.0*PIE*(x+shift_x)/lambda)+sin(2.0*PIE*(y+shift_y)/lambda)+sin(2.0*PIE*(z+shift_z)/lambda));
						}
					}
				} else
				{
					int mz=In[0]->Get_int(sub[2],0);
					if (mz<1 || mz>Lat[0]->MZ)
					{
						success=false;
						cout <<"expecting in 'mon : " + name + " : fluctuation_potentials : '  z-coordinate to be in z-range "<<endl;
						if (success && labda ==0)
						{
							cout <<"fluctutions set " << Amplitude << endl;
							Real shift_x,shift_y;
							for (int lambda_x=2; lambda_x <=MX; lambda_x*=2)
							for (int lambda_y=2; lambda_y <=MY; lambda_y*=2)
							{
								shift_x = rand() % lambda_x;
								shift_y = rand() % lambda_y;
								for (int x=0; x<MX; x++) for (int y=0; y<MY; y++)
								u_ext[x*JX+y*JY+mz]+=Amplitude*(sin(2.0*PIE*(x+shift_x)/lambda_x)+sin(2.0*PIE*(y+shift_y)/lambda_y));
							}
						} else //5
						{
							if (labda>0)
							{
								cout <<"fluctuation wavelength set to " << labda << " and amplitude to " << Amplitude << endl;
								for (int x=0; x<MX; x++) for (int y=0; y<MY; y++) u_ext[x*JX+y*JY+mz]+=Amplitude*(sin(2.0*PIE*(x)/labda)+sin(2.0*PIE*(y)/labda));
							}
						}
					}
				}
			}
			break;
			default:
			break;
		}
	}
	return success;
}


bool Segment::CheckInput(int start_) {
if (debug) cout <<"CheckInput in Segment " + name << endl;
	bool success;
	start=start_;
	block=false;
	unique=true;
	chi_var_seg=-1;
	chi_var_state=-1;
	seg_nr_of_copy=-1;
	state_nr_of_copy=-1;
	ns=1;
	//string s;
	vector<string>options;
	guess_u=0;
	n_pos=0;

	fixedPsi0=false;
	success = In[0]->CheckParameters("mon",name,start,KEYS,PARAMETERS,VALUES);
	if(success) {
		if (GetValue("var_pos").size()>0) var_pos=In[0]->Get_int(GetValue("var_pos"),0);


		options.push_back("free");
		options.push_back("pinned");
		options.push_back("frozen");
		options.push_back("tagged");
		options.push_back("clamp");
		freedom = In[0]->Get_string(GetValue("freedom"),"free");
		if (!In[0]->InSet(options,freedom)) {
			cout << "Freedom: '"<< freedom  <<"' for mon " + name + " not recognized. "<< endl;
			cout << "Freedom choices: free, pinned, frozen, tagged, clamp " << endl; success=false;
		}

		if (freedom =="free") {
			if (GetValue("frozen_range").size()>0||GetValue("pinned_range").size()>0 || GetValue("tagged_range").size()>0 ||
			GetValue("frozen_filename").size()>0 || GetValue("pinned_filename").size()>0 || GetValue("tagged_filename").size()>0) {
					if (start==1) {success=false; cout <<"In mon " + name + " you should not combine 'freedom : free' with 'frozen_range' or 'pinned_range' or 'tagged_range' or corresponding filenames." << endl;
				}
			}
		}


		valence =0;
		if (GetValue("valence").size()>0) {
			valence=In[0]->Get_Real(GetValue("valence"),0);
			if (valence<-10 || valence > 10) cout <<"For mon " + name + " valence value out of range -10 .. 10. Default value used instead" << endl;
		}
		epsilon=80;
		if (GetValue("epsilon").size()>0) {
			epsilon=In[0]->Get_Real(GetValue("epsilon"),80);
			if (epsilon<1 || epsilon > 250) cout <<"For mon " + name + " relative epsilon value out of range 1 .. 250. Default value 80 used instead" << endl;
		}
		if (valence !=0) {
			if (Lat[0]->bond_length <1e-12 || Lat[0]->bond_length > 1e-8) {
				success=false;
				if (Lat[0]->bond_length==0) cout << "When there are charged segments, you should set the bond_length in lattice to a reasonable value, e.g. between 1e-102... 1e-8 m " << endl;
				else cout <<"Bond length is out of range: 1e-12..1e-8 m " << endl;
			}
		}
		if (GetValue("e.psi0/kT").size()>0) {
			PSI0=0;
			fixedPsi0=true;
			PSI0=In[0]->Get_Real(GetValue("e.psi0/kT"),0);
			if (PSI0!=0 && valence !=0) {
				success=false;
				cout <<"You can set only 'valence' or 'e.psi0/kT', but not both " << endl;
			}
			if (PSI0!=0 && freedom!="frozen") {
				success=false;
				cout <<"You can not set potential on segment that has not freedom 'frozen' " << endl;
			}
			if (PSI0 <-25 || PSI0 > 25) {
				success=false;
				cout <<"Value for dimensionless surface potentials 'e.psi0/kT' is out of range -25 .. 25. Recall the value of 1 at room temperature is equivalent to approximately 25 mV " << endl;
			}
		}
	}

	int length = state_name.size();
	if (length >0 && (freedom == "frozen"||freedom=="tagged"||freedom=="clamp")) {
		success=false;
		cout <<" When freedom = {frozen,tagged,clamp} a 'mon' can not have multiple internal states; status violated for mon " << name << endl;
	}

	length=chi_name.size();

	Real Chi;
	for (int i=0; i<length; i++) {
		Chi=-999;
		if (GetValue("chi-"+chi_name[i]).size()>0) {
			Chi=In[0]->Get_Real(GetValue("chi-"+chi_name[i]),Chi);
			if (Chi==-999) {success=false; cout <<" chi value: chi("<<name<<","<<chi_name[i]<<") = "<<GetValue("chi-"+chi_name[i]) << "not valid." << endl; }
			if (name==chi_name[i] && Chi!=0) {if (Chi!=-999) cout <<" chi value for chi("<<name<<","<<chi_name[i]<<") = "<<GetValue("chi-"+chi_name[i]) << "value ignored: set to zero!" << endl; Chi=0;}

		}
		chi[i]=Chi;
	}

	if (GetValue("fluctuation_potentials").size()>0) {
		if (GetValue("fluctuation_wavelength").size()>0) {
				labda=In[0]->Get_int(GetValue("fluctuation_wavelength"),0);
				if (labda<1 || labda>Lat[0]->MX || labda > Lat[0]->MY || labda > Lat[0]->MZ) {
					success = false;cout <<"fluctuation_wavelength must be a positive number smaller or equal to the 'box' size" << endl;
				}
				if (!(labda ==2 || labda ==4 || labda ==8 || labda ==16 || labda ==32 || labda ==64 || labda ==128 || labda ==256 || labda ==512 ||labda ==1024)) {
					cout <<"fluctuation wavelength should be an integer 2^x, with x = 1..10" << endl;
				}
		}
		if (Lat[0]->gradients==2 || Lat[0]->gradients==1) {
			labda = Lat[0]->MY;
			labda=In[0]->Get_int(GetValue("fluctuation_wavelength"),labda);
			if (labda !=Lat[0]->MY) {
				labda=Lat[0]->MY; cout <<"fluctuation_wavelength is set to n_layers_y." << endl;
			}
			if (Lat[0]->geometry !="planar") {
				success=false; cout <<"fluctuation_potentials in 2 or 1 gradient(s) calculations only for 'planar' case." << endl;
			}
		}
		if (Lat[0]->gradients==3) {
				int MX=Lat[0]->MX;
				int MY=Lat[0]->MY;
				if (!(MX==2 || MX==4 || MX==8 || MX==16 ||MX==32 || MX==64 ||MX==128 || MX==256)) {success=false; cout << "Expecting n_layers_x to have a value 2^a with a = 1..8" << endl; }
				if (!(MY==2 || MY==4 || MY==8 || MY==16 ||MY==32 || MY==64 ||MY==128 || MY==256)) {success=false; cout << "Expecting n_layers_y to have a value 2^a with a = 1..8" << endl; }
			}
	}
	if (GetValue("fluctuation_amplitude").size()>0) {
		Amplitude = In[0]->Get_Real(GetValue("fluctuation_amplitude"),1);
		if (GetValue("fluctuation_potentials").size()==0) {
			success = false; cout <<"fluctuation_amplitude should be combined with fluctuation_potentials and optionally with fluctuation_wavelength" << endl;
		}
		if (Amplitude < 0 || Amplitude > 10) {
			success=false;  cout <<"fluctuation_amplidude sould have a value between 0 (no fluctuations) and 10. " << endl;
		}
	}




	//vector<int> constraint_z;
	//vector<Real> constraint_phi;
	//vector<Real> constraint_beta;
	if (GetValue("phi").size()>0){
		constraints=true;
		string s = GetValue("phi");
		if (Lat[0]->gradients > 1 ) {success=false; cout <<"Constraints in phi only allowed in one-gradient calculations. (for the time being)" << endl; }
		if (s=="?") {
				cout <<"Expect for the argument of 'phi' a sequence of z_values and corresponding volume fractions: " << endl;
				cout <<"(z_value_1,phi_value_1)(z_value_2,phi_value_2) ....(z_value_last,phi_value_last)" << endl;
				cout <<"'z_value' should be in range 1 .. " + Lat[0]->MX << endl;
				cout <<"'phi_value' should be in range 0 .. 1 " << endl;
			  success=false;
		}
		vector<int> open;
		vector<int> close;
		vector<string>sub;
		open.clear(); close.clear();
		if (!In[0]->EvenBrackets(s,open,close)) {
			cout << "s : " << s << endl;
			cout << "In constraints for segment " + name + " the backets are not balanced. For help use: 'mon : " +name + " : phi : ?' " << endl; success=false;
		}
		int length=open.size();
		int k=0;
		while (k<length and success) {
			string sA=s.substr(open[k]+1,close[k]-open[k]-1);
			sub.clear();
			In[0]->split(sA,',',sub);
				int zz=In[0]->Get_int(sub[0],0);
				Real RHO=In[0]->Get_Real(sub[1],-1);
				if (zz<1 || zz>Lat[0]->MX) {
					cout <<"In constraints for segment " + name + "failed to understand '" + sA + "' no valid integer found for first argument. For help use: 'mon : " +name + " : phi : ?' " << endl; success=false;
				} else constraint_z.push_back(zz);
				if (RHO<0 ||RHO>1) {
					cout <<"In constraints for segment " + name + "failed to understand '" + sA + "' no valid Real found for second argument. For help use: 'mon : " +name + " : phi : ?' " << endl; success=false;
				} else constraint_phi.push_back(RHO);
				constraint_beta.push_back(0);
			k++;
		}
	}

	if (constraints) {
		int length = constraint_z.size();
		cout <<"constraints for monomer type " + name << endl;
		for (int i=0; i<length; i++){
			cout << "constraint z =" << constraint_z[i] << " constraint phi = " << constraint_phi[i] << endl;
		}
	}

	//valence=In[0]->Get_Real(GetValue("valence"),0);
	bool HMD=false;
	r=(int*) malloc(6*sizeof(int)); std::fill(r,r+6,0);
	if (success) success=ParseFreedoms(HMD);
	if(n_pos>0) free(H_P);
	free(r);
	return success;
}

Real Segment::Get_g(int ii) {
	if (debug) cout <<"Get_g" + name << endl;
	return constraint_phi[ii]/phi[constraint_z[ii]]-1.0;
}

void Segment::Put_beta(int ii, Real BETA) {
	constraint_beta[ii]=BETA; //just to enable it to be outputted.
	u[constraint_z[ii]] +=BETA;
}

Real Segment::Volume_particles() {
	int volume=0;
	if (freedom=="frozen") Sum(volume,H_MASK,Lat[0]->M);
	//cout <<"volume_particles of type "+name + "= " << volume << endl;
	return 1.0*volume;
}

bool Segment::PutAdsorptionGuess(Real chi,int* Mask) {
if (debug) cout <<"PutAdsorptionGuess" + name << endl;
	bool success=true;
	Real lambda;
	if (Lat[0]->lattice_type==hexagonal) lambda=0.25; else lambda=1.0/6.0;
	int gradients=Lat[0]->gradients;
	int M=Lat[0]->M;
	int MX=Lat[0]->MX;
	int MY=Lat[0]->MY;
	int MZ=Lat[0]->MZ;
	int JX=Lat[0]->JX;
	int JY=Lat[0]->JY;
	switch(gradients) {
		case 1:
			for (int x=1; x<MX+1; x++)
				if (Mask[x-1]==1 ||Mask[x+1]==1)
					u[x]=-lambda*chi;
			break;
		case 2:
			for (int x=1; x<MX+1; x++) for (int y=1; y<MY+1; y++)
				if (Mask[(x-1)*JX+y]==1 ||Mask[(x+1)*JX+y]==1 || Mask[x*JX+y-1]==1 || Mask[x*JX+y+1]==1)
					u[x*JX+y]=-lambda*chi;
			break;
		case 3:
			for (int x=1; x<MX+1; x++) for (int y=1;y<MY+1; y++) for (int z=1; z<MZ+1; z++)
				if (Mask[(x-1)*JX+y*JY+z]==1 ||Mask[(x+1)*JX+y*JY+z]==1 || Mask[x*JX+(y-1)*JY+z]==1 || Mask[x*JX+(y+1)*JY+z]==1 || Mask[x*JX+y*JY+z-1]==1 || Mask[x*JX+y*JY+z+1]==1)
					u[x*JX+y*JY+z]=-lambda*chi;
			break;
		default:
			break;
	}
	Boltzmann(G1,u,M);
	return success;
}

bool Segment::PutTorusPotential(int sign) {
if (debug) cout <<"PutTorusPotential " + name << endl;
	bool success=true;
	Real distance=0;
	int count=0;
	Real L=0;
	int M=Lat[0]->M;
	int MX=Lat[0]->MX;
	int MY=Lat[0]->MY;
	int JX=Lat[0]->JX;
	int R_offset=Lat[0]->offset_first_layer;
	Real R_center = R_offset + MX/2.0;
	Real R=R_center/sqrt(2.0);
	if (R<MX/2.0) {
		for (int x=1; x<MX+1; x++) for (int y=1; y<MY+1; y++) {
	 		distance = sqrt((MX/2.0 -x)*(MX/2.0-x) + y*y);
			if ((distance-R)*(distance-R)<8) {
				u[x*JX+y]=-log(1.8)*sign; count++;
				//cout << "at x " << x << "and y " << y << "potential is set" << endl;
			}
		}
		int ylast=0;
		int xlow=0,xhigh=0;
		for (int y=1; y<MY+1; y++) {
			bool neg_found=false;
			bool pos_found=false;
			for (int x=1; x<MX+1; x++) {
				distance = sqrt((MX/2.0 -x)*(MX/2.0-x) + y*y);
				if (!neg_found && (distance-R<0)) {
					xlow = x; ylast=y;
					neg_found=true; L+=Lat[0]->L[x*JX+y];
					//cout <<"x " << x << "y " << y << endl;
				}
				if (neg_found && !pos_found && (distance-R)>0) {
					xhigh=x; ylast=y;
					pos_found=true; L+=Lat[0]->L[x*JX+y];
					//cout <<"x " << x << "y " << y << endl;
				}
			}
		}
		for (int x=xlow+1; x<xhigh; x++) L+=Lat[0]->L[x*JX+ylast];
		if (sign>0) cout << "Measured area is " << L << endl;
		cout << "For segment " << name << ", 'torus potentials' set at " << count << "coordinates" << endl;
		Boltzmann(G1,u,M);
	} else {
		success=false; cout <<" Probably the 'offset_first_layer' is too large so that the torus does not fit into the system.... Inital guess for torus is failing...."<<endl;
	}
	return success;
}

bool Segment::PutMembranePotential(int sign) {
if (debug) cout <<"PutMembranePotential " + name << endl;
	bool success=true;
	int fjc=Lat[0]->fjc;
	for (int x=1; x<4*fjc; x++) u[x]=-log(1.8)*sign;
	Boltzmann(G1,u,Lat[0]->M);
	return success;
}

void Segment::SetPhiSide(){
if (debug) cout <<"SetPhiSide in Segment " + name << endl;
	int M=Lat[0]->M;
	if (ns==1) {
		if (freedom !="frozen") Lat[0]->set_bounds(phi);
		Lat[0]->Side(phi_side,phi,M);
	} else {
		for (int i=0; i<ns; i++) {
			Times(phi_state+i*M,alpha+i*M,phi,M);
			Lat[0]->Side(phi_side+i*M,phi_state+i*M,M);
			state_phibulk[i]=phibulk*state_alphabulk[i];
		}
	}
	if (ns>1) {
		//cout <<" seg " << name << endl;
		//for (int i=0; i<ns; i++) {
		//	cout <<"alphabulk " << state_alphabulk[i] << "valence " << state_valence[i] << endl;
		//	for (int z=0; z<M; z++) cout <<  "phi("<< z<< ")= " << phi_state[z+i*M] << endl;
		//}
	}

}

bool Segment::PutVarInfo(string Var_type_, string Var_target_, Real Var_target_value_){
if (debug) cout << "Segment::PutVarInfo " << endl;
	bool success=true;

	int length_mon,length_state;
	int i;

	Var_target=-1;
	chi_var_seg=-1;
	chi_var_state=-1;
	Var_type="";
	if (Var_type_=="scan"){
		Var_type="scan";
		if (Var_target_=="valence") {Var_target=0; Var_start_value=valence;}
		if (Var_target_=="ePsi0/kT") {Var_target=1; Var_start_value=PSI0;}
		if (Var_target_ =="fluctuation_amplitude") {Var_target=3; Var_start_value=Amplitude;}
		if (Var_target_=="var_pos") {Var_target=4; Var_start_value=var_pos;}
		if (Var_target ==-1) {
			vector<string>sub;
			In[0]->split(Var_target_,'-',sub);
			if (sub.size()==2) {
				if (sub[0]=="chi") {
					length_mon=In[0]->MonList.size();
					for (i=0; i<length_mon; i++) {
						if (sub[1]==In[0]->MonList[i]) {
							Var_target=2; Var_start_value=chi[i];
							chi_var_seg=i;
						}
					}
					if (Var_target!=2) {
						length_state=In[0]->StateList.size();
						for (i=0; i<length_state; i++) {
							if (sub[1]==In[0]->StateList[i]) {
								Var_target = 2; Var_start_value=chi[i+length_mon];
								chi_var_state=i;
							}
						}
					}
					if (Var_target !=2) cout <<"In var: trying to read " + Var_target_ + " failed, because neither Seg " + sub[1] + " nor State " + sub[1] + " were found" << endl;
				}
			}
		}
	}
	if (Var_target<0) {success=false; cout <<"In var: for segment you can 'scan' {valence, ePsi0/kT, amplitude, or a chi-value 'chi-X' with 'X' valid mon/state : name} "<<endl; }
	return success;
}

int Segment::PutVarScan(Real step, Real end_value, int steps, string scale_) {
if (debug) cout << "Segment::PutVarScan " << endl;
	num_of_steps=-1;
	scale=scale_;
	Var_end_value=end_value;
	if (scale=="exponential") {
		Var_steps=steps; Var_step = 0;
		if (steps==0) {
			cout <<"In var scan: the value of 'steps' is zero, this is not allowed" << endl;
			return -1;
		}
		if (Var_end_value*Var_start_value <0) {
			cout <<"In var scan: the product end_value*start_value < 0. This is not allowed. " << endl;
			return -1;
		}
		if (Var_end_value > Var_start_value)
			num_of_steps= steps* log10 (Var_end_value/Var_start_value);
		else
			num_of_steps= steps* log10 (Var_start_value/Var_end_value);

	} else {
		Var_steps=0; Var_step=step;
		if (step==0) {
			cout <<"In var scan : of segment variable, the value of step can not be zero" << endl;
			return -1;
		}
		num_of_steps=(Var_end_value-Var_start_value)/step;

		if (num_of_steps<0) {
			cout<<"In var scan : (end_value-start_value)/step is negative. This is not allowed. Try changing the sign of the 'step'." << endl;
			return -1;
		}
	}


	return num_of_steps;
}

bool Segment::UpdateVarInfo(int step_nr) {
if (debug) cout << "Segment::UpdateVarInfo() " << endl;
	bool success=true;
	int length;
	switch(Var_target) {
		case 0:
			if (scale=="exponential") {
				if (valence <0)	{
					valence=-pow(10,(1-1.0*step_nr/num_of_steps)*log10(-Var_start_value)+ (1.0*step_nr/num_of_steps)*log10(-Var_end_value));
				} else {
					valence= pow(10,(1-1.0*step_nr/num_of_steps)*log10( Var_start_value)+ (1.0*step_nr/num_of_steps)*log10( Var_end_value));
				}
			} else {
				valence=Var_start_value+step_nr*Var_step;
			}
			break;
		case 1:
			if (scale=="exponential") {
				if (PSI0 <0)	{
					PSI0=-pow(10,(1-1.0*step_nr/num_of_steps)*log10(-Var_start_value)+ (1.0*step_nr/num_of_steps)*log10(-Var_end_value));
				} else {
					PSI0= pow(10,(1-1.0*step_nr/num_of_steps)*log10( Var_start_value)+ (1.0*step_nr/num_of_steps)*log10( Var_end_value));
				}
			} else {
				PSI0=Var_start_value+step_nr*Var_step;
			}
			break;
		case 2:
			if (scale=="exponential") {
				cout <<"In var of chi-parameter, only linear scale is implemented" << endl; success=false;
			} else {
				length = In[0]->MonList.size();
				if (chi_var_seg>-1) {
					chi[chi_var_seg]= Var_start_value+step_nr*Var_step;
				}
				if (chi_var_state>-1) {
					chi[length+chi_var_state] =Var_start_value+step_nr*Var_step;
				}
				//chi_value = Var_start_value+step_nr*Var_step;
			}
			break;
		case 3:
			if (scale=="exponential") {
				Amplitude= pow(10,(1-1.0*step_nr/num_of_steps)*log10(Var_start_value)+(1.0*step_nr/num_of_steps)*log10(Var_end_value));
			} else {
				Amplitude=Var_start_value+step_nr*Var_step;
			}
			break;
		case 4:
			if (scale=="exponential") {
				var_pos=   pow(10,(1-1.0*step_nr/num_of_steps)*log10(Var_start_value)+(1.0*step_nr/num_of_steps)*log10(Var_end_value));
			} else {
				var_pos=Var_start_value+step_nr*Var_step;
			}
			break;
		default:
			break;
	}
	return success;
}

bool Segment::ResetInitValue() {
if (debug) cout << "Segment::ResetInitValue() " << endl;
	bool success=true;
	int length;
	switch(Var_target) {
		case 0:
			valence=Var_start_value;
			break;
		case 1:
			PSI0=Var_start_value;
			break;
		case 2:
			length = In[0]->MonList.size();
			if (chi_var_seg>-1) {
				chi[chi_var_seg]= Var_start_value;
			}
			if (chi_var_state>-1) {
				chi[length+chi_var_state] =Var_start_value;
			}
			break;
		case 3:
			Amplitude=Var_start_value;
			break;
		case 4:
			var_pos=Var_start_value;
			break;
		default:
			cout <<"program error in Seg:ResetInitValue "<<endl;
			break;
	}
	return success;
}

void Segment::PutValue(Real X) {
if (debug) cout << "Segment::PutValue() " << endl;
	int length;
	switch(Var_target) {
		case 0:
			valence=X;
			break;
		case 1:
			PSI0=X;
			break;
		case 2:
			length = In[0]->MonList.size();
			if (chi_var_seg>-1) {
				chi[chi_var_seg]= X;
			}
			if (chi_var_state>-1) {
				chi[length+chi_var_state] =X;
			}
			break;
		case 3:
			Amplitude=X;
			break;
		case 4:
			var_pos=X;
			break;
		default:
			cout <<"program error in Segment:PutValue "<<endl;
			break;
	}
}

Real Segment::GetValue() {
if (debug) cout << "Segment::GetValue() " << endl;
	Real X=0;
	int length;
	switch(Var_target) {
		case 0:
			X=valence;
			break;
		case 1:
			X=PSI0;
			break;
		case 2:
			length = In[0]->MonList.size();
			if (chi_var_seg>-1) {
				X=chi[chi_var_seg];
			}
			if (chi_var_state>-1) {
				X=chi[length+chi_var_state];
			}

			break;
		case 3:
			X=Amplitude;
			break;
		case 4:
			X=var_pos;
			break;
		default:
			cout <<"program error in Segment:GetValue "<<endl;
			break;
	}
	return X;
}


bool Segment::GetClamp(string filename) {
	bool success=true;
	string line_;
	string NN,X,Y,Z;
	int pos;
	int N=0;
	mx=my=mz=0;
	n_box=-1;
	int i=0;
	ifstream in_file;
	in_file.open(filename.c_str());
	if (in_file.is_open()) {
		while (in_file) {
			i++;
			if (i==1) {
				n_box++;
				NN.clear();
				in_file >> line_ >> NN >> line_ >> X >> Y >> Z >> line_ >> line_ >> line_ >> line_;
				if (NN.size()>0) {
				if (!In[0]->Get_int(NN,pos,"")) { success=false; cout<<" length of 'fragment' not an integer " << endl; }
				else {if (N>0) {if (N!=pos) {cout <<"lengths of fregment are not equal " << endl; success=false;}} else N = pos;}
				if (!In[0]->Get_int(X,pos,"")) {success=false; cout <<"X-value for sub_box size is not integer " << endl;  }
				else {if (mx>0) {if (mx !=pos) {success=false; cout <<"We can deal with only one sub-box size" << endl;}} else mx=pos; }
				if (!In[0]->Get_int(Y,pos,"")) {success=false; cout <<"Y-value for sub_box size is not integer " << endl;  }
				else {if (my>0) {if (my !=pos) {success=false; cout <<"We can deal with only one sub-box size" << endl;}} else my=pos; }
				if (!In[0]->Get_int(Z,pos,"")) {success=false; cout <<"Z-value for sub_box size is not integer " << endl;  }
				else {if (mz>0) {if (mz !=pos) {success=false; cout <<"We can deal with only one sub-box size" << endl;}} else mz=pos; }
				}
			}
			if (i==2) {in_file>>X>>Y>>Z;
				if (!In[0]->Get_int(X,pos,"")) {success =false; cout << " X-pos of particle sub_box nr " << n_box << " not integer"  << endl; }
				else bx.push_back(pos);
				if (!In[0]->Get_int(Y,pos,"")) {success =false; cout << " Y-pos of particle sub_box nr " << n_box << " not integer"  << endl; }
				else by.push_back(pos);
				if (!In[0]->Get_int(Z,pos,"")) {success =false; cout << " Z-pos of particle sub_box nr " << n_box << " not integer"  << endl; }
				else bz.push_back(pos);
			}
			if (i==3) {in_file>>X>>Y>>Z;
				if (!In[0]->Get_int(X,pos,"")) {success =false; cout << " X-pos of particle pos 1 of sub_box nr " << n_box << " not integer"  << endl; }
				else px1.push_back(pos);
				if (!In[0]->Get_int(Y,pos,"")) {success =false; cout << " Y-pos of particle pos 1 of sub_box nr " << n_box << " not integer"  << endl; }
				else py1.push_back(pos);
				if (!In[0]->Get_int(Z,pos,"")) {success =false; cout << " Z-pos of particle pos 1 of sub_box nr " << n_box << " not integer"  << endl; }
				else pz1.push_back(pos);
			}
			if (i==4) {i=0;in_file>>X>>Y>>Z;
				if (!In[0]->Get_int(X,pos,"")) {success =false; cout << " X-pos of particle pos 2 of sub_box nr " << n_box << " not integer"  << endl; }
				else px2.push_back(pos);
				if (!In[0]->Get_int(Y,pos,"")) {success =false; cout << " Y-pos of particle pos 2 of sub_box nr " << n_box << " not integer"  << endl; }
				else py2.push_back(pos);
				if (!In[0]->Get_int(Z,pos,"")) {success =false; cout << " Z-pos of particle pos 2 of sub_box nr " << n_box << " not integer"  << endl; }
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
	KEYS.push_back("chi-" + new_name);
	chi_name.push_back(new_name);
	chi.push_back(-999);
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
	int M = Lat[0]->M;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	ints.clear();
	ints_value.clear();
	Reals.clear();
	Reals_value.clear();
	push("freedom",freedom);
	push("valence",valence);
	Real theta = Lat[0]->WeightedSum(phi);
	push("theta",theta);
	theta_exc=0;
	theta_exc=theta-Lat[0]->Accesible_volume*phibulk;
	push("theta_exc",theta_exc);
	push("phibulk",phibulk);
	push("fluctuation_amplitude",Amplitude);
	push("fluctuation_wavelength",labda);
	push("var_pos",var_pos);
	if (theta_exc!=0) {
		M1=0;
		M1=Lat[0]->Moment(phi,phibulk,1)/theta_exc;
	 	M2=0;
		M2=Lat[0]->Moment(phi,phibulk,2)/theta_exc;
		push("1st_M_phi_z",M1);
		push("2nd_M_phi_z",M2);
		Fl = (M2-M1*M1);
		if (Fl >0) Fl = sqrt(Fl); else Fl=0;
		push("fluctuations",Fl);
	}
	if (ns>1) {
		state_theta.clear();
		for (int i=0; i<ns; i++){
			push("alphabulk-"+state_name[i],state_alphabulk[i]);
			push("valence-"+state_name[i],state_valence[i]);
			push("phibulk-"+state_name[i],state_phibulk[i]);
			theta=Lat[0]->WeightedSum(phi_state+i*M);
			state_theta.push_back(theta);
			push("theta-"+state_name[i],theta);
			push("theta_exc-"+state_name[i],theta-Lat[0]->Accesible_volume*state_phibulk[i]);
		}
	}
	int length=chi_name.size();
	for (int i=0; i<length; i++) push("chi-"+chi_name[i],chi[i]);
	if (fixedPsi0) push("Psi0",PSI0);
	if (freedom=="pinned") push("range",GetValue("pinned_range"));
	if (freedom=="frozen") push("range",GetValue("frozen_range"));
	if (freedom=="tagged") push("range",GetValue("tagged_range"));
	string profile="profile;0"; push("phi",profile);

	profile="profile;1"; push("G1",profile);
	if (Lat[0]->gradients==3) {
		profile="profile;2"; push("phi[z]",profile);

	}
	int k=1;
	string s;
	string str;
	if (ns >1) {
		for (int i=0; i<ns; i++) { k++;
			stringstream ss; ss<<k; str=ss.str();
			s="profile;"+str; push("phi-"+state_name[i],s);
		}
		for (int i=0; i<ns; i++) { k++;
			stringstream ss; ss<<k; str=ss.str();
			s="profile;"+str; push("alpha-"+state_name[i],s);
		}
		for (int i=0; i<ns; i++) { k++;
			stringstream ss; ss<<k; str=ss.str();
			s="profile;"+str; push("u-"+state_name[i],s);
		}
	}


#ifdef CUDA
	TransferDataToHost(H_phi, phi, M);
	//TransferDataToHost(H_u, u, M);
#endif
}

Real* Segment::GetPointer(string s, int &SIZE) {
if (debug) cout <<"Get Pointer for segment " + name << endl;
	vector<string> sub;
	int M=Lat[0]->M;
	SIZE=Lat[0]->M;
	In[0]->split(s,';',sub);
	if (sub[0]=="profile") {

	if (sub[1]=="0") {
		if (freedom=="frozen") {
			Cp(phi,MASK,M);
		} else Lat[0]->set_bounds(phi);
		return phi;
	}
	if (sub[1]=="1") return G1;
        if (sub[1]=="2") {
		int MX=Lat[0]->MX;
		int MY=Lat[0]->MY;
		int MZ=Lat[0]->MZ;
		int JX=Lat[0]->JX;
		int JY=Lat[0]->JY;
		Real Sum;
		Zero(phi_side,M); //phi_side is reused because this array is no longer needed (hopefully....).
		for (int z=0; z<MZ; z++) {
			Sum=0;
			for (int x=1; x<MX+1; x++) for (int y=1; y<MY+1; y++)
				Sum +=phi[x*JX+y*JY+z];
			Sum /=MX*MY;
			phi_side[JX+JY+z]=Sum;
		}
		return phi_side;
	}
	if (ns>1) {
		for (int i=0; i<ns; i++) {
			stringstream ss; ss<<i+2; string str=ss.str();
			if (sub[1]==str) return phi_state+i*M;
		}
		for (int i=0; i<ns; i++) {
			stringstream ss; ss<<i+ns+2; string str=ss.str();
			if (sub[1]==str) return alpha+i*M;
		}
		for (int i=0; i<ns; i++) {
			stringstream ss; ss<<i+2*ns+2; string str=ss.str();
			if (sub[1]==str) return u+i*M;
		}
	}


	} else {//sub[0]=="vector" ..do not forget to set SIZE before returning the pointer.
	}
	return NULL;
}
int* Segment::GetPointerInt(string s, int &SIZE) {
if (debug) cout <<"GetPointerInt for segment " + name << endl;
	vector<string> sub;
	SIZE=Lat[0]->M;
	In[0]->split(s,';',sub);
	if (sub[0]=="array") {// set SIZE and return int pointer.
	}
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
void Segment::UpdateValence(Real*g, Real* psi, Real* q, Real* eps,bool grad_epsilon) {
	int M=Lat[0]->M;
	if (fixedPsi0) {
		OverwriteC(psi,MASK,PSI0,M);
		Lat[0]->set_M_bounds(psi);
		Lat[0]->UpdateQ(g,psi,q,eps,MASK,grad_epsilon);
	}

}
int Segment::AddState(int id_,Real alphabulk,Real valence,bool fixed) {
if (debug) cout <<"AddState " << id_ <<" to seg " << name << endl;
	int length = state_name.size();
	int state_number=-1;
	bool found=false;
	int ID=id_;
	for (int k=0; k<length; k++) {
		if (state_id[k]==ID) {
			state_name[k]=In[0]->StateList[ID];
			found=true;
			state_alphabulk[k]=alphabulk;
			state_valence[k]=valence;  
			if (fixed) state_change[k]=false; else state_change[k]=true;
			state_number=k;
		}
	}
	if (!found) {
		state_id.push_back(ID);
		state_name.push_back(In[0]->StateList[ID]);
		state_alphabulk.push_back(alphabulk);
		state_phibulk.push_back(0);
		state_valence.push_back(valence);
		if (fixed) state_change.push_back(false); else state_change.push_back(true);
		state_number=state_change.size()-1;
		length=In[0]->StateList.size();
		for (int k=0; k<length; k++) {if (name==In[0]->StateList[k]) state_nr.push_back(k);}
	}

	if (valence !=0) {
		if (Lat[0]->bond_length <1e-10 || Lat[0]->bond_length > 1e-8) {
			if (Lat[0]->bond_length==0) cout << "When there are charged states, you should set the bond_length in lattice to a reasonable value, e.g. between 1e-10 ... 1e-8 m " << endl;
			else cout <<"Bond length is out of range: 1e-10..1e-8 m " << endl;
		}
	}

	return state_number;
}

/*
void Segment::PutAlpha(Real* x,int &xi){
if (debug) cout <<"PutAlpha " << name << endl;

	int n_s;
	if (ns==1) n_s=0; else n_s=ns;
	if (n_s==0) return ;
	Real sum_alpha=0;
	Real fixed_value=0;
	int niv=1;
	int n_free=n_s-1;

	for (int i=0; i<n_s; i++) {
		if (!state_change[i]) {fixed_value+=state_alphabulk[i]; n_free--;} else state_alphabulk[i]=0; 	}

	for (int i=0; i<n_s; i++) {

		if (state_change[i] && niv>0 && n_free>0) {
			state_alphabulk[i]=0.5*(1.0+tanh(2.0*x[xi]))*(1.0-fixed_value);
//cout <<"iv : "<< xi <<" for "<< name << endl;  
			xi++; niv--;
		}
		sum_alpha+=state_alphabulk[i];
	}
	for (int i=0; i<n_s; i++) {
		if (state_alphabulk[i]==0) state_alphabulk[i]=1.0-sum_alpha;
	}
}
*/

bool Segment::PutAlpha(Real alpha) { //expected to replace other method with same name.
	bool success=true;
	Real fixed_value=0;
	Real sum_alpha=0;
	int n_s;
	if (ns==1) n_s=0; else n_s=ns;
	if (n_s==0) return false;
 

	for (int i=0; i<n_s; i++) {
		if (!state_change[i]) fixed_value+=state_alphabulk[i];} 
	for (int i=0; i<n_s; i++) {
		if (ItState !=i && state_change[i]) state_alphabulk[i]=0; 
	}


	state_alphabulk[ItState]*=alpha; 

	for (int i=0; i<n_s; i++) {
		sum_alpha+=state_alphabulk[i];
	}
	for (int i=0; i<n_s; i++) {
		if (state_alphabulk[i]==0) state_alphabulk[i]=1.0-sum_alpha;
		if (state_alphabulk[i]<0 || state_alphabulk[i]>1) {
			success=false; cout <<"In Segment::PutAlpha, alphabulk out of bounds " << endl; 
		}
	}

	//cout << endl; 
	//cout <<"Seg " << name << endl;
	//for (int i=0; i<n_s; i++) cout << "state[" << i << "].alpha_bulk = " << state_alphabulk[i] << endl;  

	return success;	
}

bool Segment::CanBeReached(int x0, int y0, int z0, int Ds) {
	int gradients=Lat[0]->gradients;
	int MX = Lat[0]->MX;
	int MY = Lat[0]->MY;
	int MZ = Lat[0]->MZ;
	int JX = Lat[0]->JX;
	int JY = Lat[0]->JY;
	int JZ = Lat[0]->JZ;
	int x=0,y=0,z=0;
	switch (gradients) {
		case 3:
			for (z=1; z<MZ+1; z++)
		case 2:
			for (y=1; y<MY+1; y++)
		case 1:
			for (x=1; x<MX+1; x++){
				if (MASK[x*JX+y*JY+z*JZ]==1) {
					if (abs(x-x0)+abs(y-y0)+abs(z-z0) < Ds) return true;
				}

			}
	}
	return false;
}
