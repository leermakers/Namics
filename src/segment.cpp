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
	Amplitude=0;
}
Segment::~Segment() {
if (debug) cout <<"Segment destructor " + name << endl;
	DeAllocateMemory();
}

void Segment::DeAllocateMemory(void){
if (debug) cout << "In Segment, Deallocating memory " + name << endl;

//if (n_pos>0) cout <<"problem for n_pos " <<endl;
	if(n_pos>0)
	free(H_P);
	if (freedom != "free"){
		 free(r);
	}
	free(H_u);
	free(H_u_ext);
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
}

void Segment::AllocateMemory() {
if (debug) cout <<"Allocate Memory in Segment " + name << endl;
	int M=Lat[0]->M;
	ns=state_name.size(); if (ns==0) ns=1;
	H_u = (Real*) malloc(M*ns*sizeof(Real));
	H_phi_state = (Real*) malloc(M*ns*sizeof(Real));
	H_phi = (Real*) malloc(M*sizeof(Real));
	H_u_ext = (Real*) malloc(M*sizeof(Real));
	H_alpha=(Real*) malloc(M*ns*sizeof(Real));
	H_Zero(H_u,M*ns);
	H_Zero(H_phi,M);
	if (freedom=="free") {
		H_MASK = (int*) malloc(M*sizeof(int));
		H_Zero(H_MASK,M);
	}
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
	u_ext=H_u_ext;	KEYS.push_back("fluctuation_coordinates");
	phi_state = H_phi_state;
	alpha=H_alpha;
	G1 = (Real*)malloc(M*sizeof(Real));
	phi_side = (Real*)malloc(M*ns*sizeof(Real));
	Zero(G1,M);
	Zero(phi_side,ns*M);
#endif
}

bool Segment::PrepareForCalculations(int* KSAM) {
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
	if (freedom=="tagged") {Cp(G1,MASK,M);
	}
	if (!(freedom ==" frozen" || freedom =="tagged")) Times(G1,G1,KSAM,M);

	if (GetValue("fluctuation_potentials").size()>0) {
		string s = GetValue("fluctuation_potentials");
		vector<string> sub;
		In[0]->split(s, ',', sub);
		if (sub.size()<3) {success=false; cout <<"expecting in 'mon : " + name + " : fluctuation_potentials : '  coordinate info, such as: x,y,5 or x,y,z "<<endl; }
		else {
			if (sub[0] != "x" || sub[1] != "y" ) {success=false; cout <<"expecting in 'mon : " + name + " : fluctuation_potentials : '  first two coordinates to be : x,y  "<<endl; }
			int MX=Lat[0]->MX;
			int JX=Lat[0]->JX;
			int MY=Lat[0]->MY;
			int JY=Lat[0]->JY;
			int M=Lat[0]->M;
			Zero(u_ext,M);
			if (sub[2] == "z") {
				int MZ=Lat[0]->MZ;
				if (!(MZ==2 || MZ==4 || MZ==8 || MZ==16 ||MZ==32 || MZ==64 ||MZ==128 || MZ==256)) {success=false; cout << "Expecting n_layers_z to have a value 2^a with a = 1..8" << endl; }
				if (success) {
					Real shift_x,shift_y,shift_z;
					for (int lambda_x=2; lambda_x <=MX; lambda_x*=2)
					for (int lambda_y=2; lambda_y <=MY; lambda_y*=2)
					for (int lambda_z=2; lambda_z <=MZ; lambda_z*=2){
						shift_x = rand() % lambda_x;
						shift_y = rand() % lambda_y;
						shift_z = rand() % lambda_z;
						for (int x=0; x<MX; x++) for (int y=0; y<MY; y++) for (int z=0; z<MZ; z++) {
							u_ext[x*JX+y*JY+MZ]+=Amplitude*(sin(2.0*PIE*(x+shift_x)/lambda_x)+sin(2.0*PIE*(y+shift_y)/lambda_y)+sin(2.0*PIE*(z+shift_z)/lambda_z));
						}
					}
				}
			} else {
				int mz=In[0]->Get_int(sub[2],0);
				if (mz<1 || mz>Lat[0]->MZ) {success=false; cout <<"expecting in 'mon : " + name + " : fluctuation_potentials : '  z-coordinate to be in z-range "<<endl; }
				if (success) {
					Real shift_x,shift_y;
					for (int lambda_x=2; lambda_x <=MX; lambda_x*=2)
					for (int lambda_y=2; lambda_y <=MY; lambda_y*=2){
						shift_x = rand() % lambda_x;
						shift_y = rand() % lambda_y;
						for (int x=0; x<MX; x++) for (int y=0; y<MY; y++) {
							u_ext[x*JX+y*JY+mz]+=Amplitude*(sin(2.0*PIE*(x+shift_x)/lambda_x)+sin(2.0*PIE*(y+shift_y)/lambda_y));
						}
					}
				}
			}
		}

	}

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

bool Segment::CheckInput(int start) {
if (debug) cout <<"CheckInput in Segment " + name << endl;
	bool success;
	block=false;
	unique=true;
	chi_var_seg=-1;
	chi_var_state=-1;
	seg_nr_of_copy=-1;
	state_nr_of_copy=-1;
	ns=1;
	string s;
	vector<string>options;
	guess_u=0;
	n_pos=0;

	fixedPsi0=false;
	success = In[0]->CheckParameters("mon",name,start,KEYS,PARAMETERS,VALUES);
	if(success) {
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
		} else {	KEYS.push_back("fluctuation_coordinates");
			r=(int*) malloc(6*sizeof(int)); std::fill(r,r+6,0);
		}

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
					cout <<"N : 20 : subbox0 : 15 15 15 pcb:[0. 0. 0.] " << endl;
					cout <<"-1 -6 -6 " << endl;
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
					for (int i=0; i<n_box; i++) {
						set.clear();
						In[0]->split(sub[i],'(',set);
						int length = set.size();
						if (length!=3) {
							success=false; cout <<" In 'clamp_info' for segment '"+name+"', for box number " << i << " the expected format (px1,py1,pz1)(px2,py2,pz2) was not found" << endl;
						} else {
/*
							coor.clear();
							In[0]->split(set[1],',',coor);
							if (coor.size()!=3) {
								success=false; cout <<" In 'clamp_info' for segment '"+name+"' for box number "<< i <<" the coordinates for the box position (bx,by,bz) not correct format. " << endl;
							} else {
								bx.push_back(In[0]->Get_int(coor[0],-10000));
								by.push_back(In[0]->Get_int(coor[1],-10000));
								bz.push_back(In[0]->Get_int(coor[2],-10000));
							}
*/

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
			phibulk=0;
			if (GetValue("frozen_range").size()>0 || GetValue("tagged_range").size()>0 || GetValue("frozen_filename").size()>0 || GetValue("tag_filename").size()>0) {
			cout<< "For mon :" + name + ", you should exclusively combine freedom : pinned with pinned_range or pinned_filename" << endl;  success=false;}
			if (GetValue("pinned_range").size()>0 && GetValue("pinned_filename").size()>0) {
				cout<< "For mon " + name + ", you can not combine pinned_range with 'pinned_filename' " <<endl; success=false;
			}
			if (GetValue("pinned_range").size()==0 && GetValue("pinned_filename").size()==0) {
				cout<< "For mon " + name + ", you should provide either pinned_range or pinned_filename " <<endl; success=false;
			}
			if (GetValue("pinned_range").size()>0) { s="pinned_range";
				n_pos=0;
				if (success) success=Lat[0]->ReadRange(r, H_P, n_pos, block, GetValue("pinned_range"),name,s);
				if (n_pos>0) {
					H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
					if (success) success=Lat[0]->ReadRange(r, H_P, n_pos, block, GetValue("pinned_range"),name,s);
				}
			}
			if (GetValue("pinned_filename").size()>0) { s="pinned";
				filename=GetValue("pinned_filename");
				n_pos=0;
				if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s);
				if (n_pos>0) {
					H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
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
					H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
					success=Lat[0]->ReadRange(r, H_P, n_pos, block, GetValue("frozen_range"),name,s);
				}
			}
			if (GetValue("frozen_filename").size()>0) { s="frozen";
				filename=GetValue("frozen_filename");
				n_pos=0;
				if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s);
				if (n_pos>0) {
					H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
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
					H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
					if (success) success=Lat[0]->ReadRange(r, H_P, n_pos, block, GetValue("tagged_range"),name,s);
				}
			}
			if (GetValue("tagged_filename").size()>0) {s="tagged";
				filename=GetValue("tagged_filename");
				n_pos=0;
			 	if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s);
				if (n_pos>0) {
					H_P=(int*) malloc(n_pos*sizeof(int)); std::fill(H_P, H_P+n_pos, 0);
			 		if (success) success=Lat[0]->ReadRangeFile(filename,H_P,n_pos,name,s);
				}
			}
		}
		if (freedom!="free") {
			H_MASK = (int*) malloc(Lat[0]->M*sizeof(int)); std::fill(H_MASK, H_MASK+Lat[0]->M,0);
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

			} else Lat[0]->CreateMASK(H_MASK,r,H_P,n_pos,block);
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
			if (Lat[0]->bond_length <1e-10 || Lat[0]->bond_length > 1e-8) {
				success=false;
				if (Lat[0]->bond_length==0) cout << "When there are charged segments, you should set the bond_length in lattice to a reasonable value, e.g. between 1e-10 ... 1e-8 m " << endl;
				else cout <<"Bond length is out of range: 1e-10..1e-8 m " << endl;
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
		if (Lat[0]->gradients<3) {
				success=false; cout<<"fluction_potentials only allowed when 'lat: * : gradients : 3'. "<<endl;
			} else {
				int MX=Lat[0]->MX;
				int MY=Lat[0]->MY;
				if (!(MX==2 || MX==4 || MX==8 || MX==16 ||MX==32 || MX==64 ||MX==128 || MX==256)) {success=false; cout << "Expecting n_layers_x to have a value 2^a with a = 1..8" << endl; }
				if (!(MY==2 || MY==4 || MY==8 || MY==16 ||MY==32 || MY==64 ||MY==128 || MY==256)) {success=false; cout << "Expecting n_layers_y to have a value 2^a with a = 1..8" << endl; }
			}
	}
	if (GetValue("fluctuation_amplitude").size()>0) {
		Amplitude = In[0]->Get_Real(GetValue("fluctuation_amplitude"),0);
		if (Amplitude < 0 || Amplitude > 1) {
			success=false;  cout <<"fluctuation_amplidude sould have a value between 0 (no fluctuations) and 1 (full fluctuations). " << endl; 
		}
	}
	
			//valence=In[0]->Get_Real(GetValue("valence"),0);
	return success;
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
	if (Var_target<0) {success=false; cout <<"In var: for segment you can 'scan' {valence, ePsi0/kT, or a chi-value 'chi-X' with 'X' valid mon/state : name} "<<endl; }
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
		num_of_steps=(Var_end_value-Var_start_value)/step+1;

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
					valence=-pow(10,(1-step_nr/num_of_steps)*log10(-Var_start_value)+ (step_nr/num_of_steps)*log10(-Var_end_value));
				} else {
					valence= pow(10,(1-step_nr/num_of_steps)*log10( Var_start_value)+ (step_nr/num_of_steps)*log10( Var_end_value));
				}
			} else {
				valence=Var_start_value+step_nr*Var_step;
			}
			break;
		case 1:
			if (scale=="exponential") {
				if (PSI0 <0)	{
					PSI0=-pow(10,(1-step_nr/num_of_steps)*log10(-Var_start_value)+ (step_nr/num_of_steps)*log10(-Var_end_value));
				} else {
					PSI0= pow(10,(1-step_nr/num_of_steps)*log10( Var_start_value)+ (step_nr/num_of_steps)*log10( Var_end_value));
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
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	push("freedom",freedom);
	push("valence",valence);
	Real theta = Lat[0]->WeightedSum(phi);
	push("theta",theta);
	push("theta_exc",theta-Lat[0]->Accesible_volume*phibulk);
	push("phibulk",phibulk);
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

	if (sub[1]=="0") return phi;
	if (sub[1]=="1") return G1;
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
		Lat[0]->set_bounds(psi);
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
int Segment::PutAlpha(Real* x,int xi){
	int xxi=xi;
	int n_s;
	if (ns==1) n_s=0; else n_s=ns;
	if (n_s==0) return xxi;
	Real sum_alpha=0;
	Real fixed_value=0;
	int niv=1;
	int n_free=n_s-1;

	for (int i=0; i<n_s; i++) {
		if (!state_change[i]) {fixed_value+=state_alphabulk[i]; n_free--;} else state_alphabulk[i]=0; 	}

	for (int i=0; i<n_s; i++) {
		if (state_change[i] && niv>0 && n_free>0) {
			state_alphabulk[i]=0.5*(1.0+tanh(10*(x[xxi])))*(1.0-fixed_value); xxi++; niv--;
		}
		sum_alpha+=state_alphabulk[i];
	}
	for (int i=0; i<n_s; i++) {
		if (state_alphabulk[i]==0) state_alphabulk[i]=1.0-sum_alpha;
	}
	return xxi;
}
