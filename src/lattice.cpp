#include "lattice.h"
Lattice::Lattice(vector<Input*> In_,string name_) :
	BC(6) // resize the boundary condition vector to 6 for Mesodyn
{ //this file contains switch (gradients). In this way we keep all the lattice issues in one file!
if (debug) cout <<"Lattice constructor" << endl;
	In=In_; name=name_;
	KEYS.push_back("gradients"); KEYS.push_back("n_layers"); KEYS.push_back("offset_first_layer");
	KEYS.push_back("geometry");
	KEYS.push_back("n_layers_x");   KEYS.push_back("n_layers_y"); KEYS.push_back("n_layers_z");
	KEYS.push_back("lowerbound"); KEYS.push_back("upperbound");
	KEYS.push_back("lowerbound_x"); KEYS.push_back("upperbound_x");
	KEYS.push_back("lowerbound_y"); KEYS.push_back("upperbound_y");
	KEYS.push_back("lowerbound_z"); KEYS.push_back("upperbound_z");
 	KEYS.push_back("bondlength");
	KEYS.push_back("ignore_site_fraction");
	KEYS.push_back("fcc_site_fraction");
  	KEYS.push_back("lattice_type");
	KEYS.push_back("stencil_full");
	//KEYS.push_back("lambda");
	//KEYS.push_back("Z");
	KEYS.push_back("FJC_choices");
	KEYS.push_back("b/l");
	//KEYS.push_back("Markov");
	KEYS.push_back("k_stiff");
	sub_box_on = 0;
	all_lattice = false;
	ignore_sites=false;
	fcc_sites=false;
	stencil_full=false;
	gradients=1;
	fjc=1;
	MX=MY=MZ=0;
	offset_first_layer=0;
}

void Lattice::DeAllocateMemory(void) {
if (debug) cout <<"DeAllocateMemory in lat " << endl;
	if (all_lattice) {
		all_lattice=false;
		if (fcc_sites) {
			free(fcc_lambda_1);
			free(fcc_lambda1);
			free(fcc_lambda0);
		}
		if (fjc==1) {
			if (gradients<3) {
				free(lambda_1);
				free(lambda1);
				free(lambda0);
				free(L);
			}
		} else {
			free(B_X1);free(B_Y1);free(B_Z1);free(B_XM);free(B_YM);free(B_ZM);
			free(L); free(LAMBDA);
			#ifdef CUDA
			cudaFree(X);
			#else
			free(X);
			#endif
		}
		if (Markov==2) {
			free(l1);
			free(l11);
			free(l_1);
			free(l_11);
			free(H);
		}
	}
}


void Lattice::AllocateMemory(void) {
if (debug) cout <<"AllocateMemory in lat " << endl;

	DeAllocateMemory();
	all_lattice=true;
	PutM();
	if (fjc>1) {
		B_X1=(int*)malloc(fjc*sizeof(int));
		B_Y1=(int*)malloc(fjc*sizeof(int));
		B_Z1=(int*)malloc(fjc*sizeof(int));
		B_XM=(int*)malloc(fjc*sizeof(int));
		B_YM=(int*)malloc(fjc*sizeof(int));
		B_ZM=(int*)malloc(fjc*sizeof(int));
	}

	switch (gradients) {
		case 3:
			if (BC[2]=="mirror") {
				if (fjc==1) BZ1=1; else {
					for (int k=0; k<fjc; k++) B_Z1[k]=2*fjc-1-k;
				}
			}
			if (BC[2]=="periodic") {
				if (fjc==1) BZ1=MZ; else {
					for (int k=0; k<fjc; k++) B_Z1[k]=MZ+fjc-k-1;
				}
			}
			if (BC[2]=="surface") {
				if (fjc==1) BZ1=0; else {
					for (int k=0; k<fjc; k++) B_Z1[k]=k;
				}
			}

			if (BC[5]=="mirror") {
				if (fjc==1) BZM=MZ; else {
					for (int k=0; k<fjc; k++) B_ZM[k]=MZ+fjc-k-1;
				}
			}
			if (BC[5]=="periodic") {
				if (fjc==1) BZM=1; else {
					for (int k=0; k<fjc; k++) B_ZM[k]=fjc+k;
				}
			}
			if (BC[5]=="surface") {
				if (fjc==1) BZM=MZ+1; else {
					for (int k=0; k<fjc; k++) B_ZM[k]=MZ+fjc+k;
				}
			}

			//Fall through
		case 2:
			if (BC[1]=="mirror") {
				if (fjc==1) BY1=1; else {
					for (int k=0; k<fjc; k++) B_Y1[k]=2*fjc-1-k;
				}
			}
			if (BC[1]=="periodic") {
				if (fjc==1) BY1=MY; else {
					for (int k=0; k<fjc; k++) B_Y1[k]=MY+fjc-k-1;
				}
			}
			if (BC[1]=="surface") {
				if (fjc==1) BY1=0; else {
					for (int k=0; k<fjc; k++) B_Y1[k]=k;
				}
			}

			if (BC[4]=="mirror") {
				if (fjc==1) BYM=MY; else {
					for (int k=0; k<fjc; k++) B_YM[k]=MY+fjc-k-1;
				}
			}
			if (BC[4]=="periodic") {
				if (fjc==1) BYM=1; else {
					for (int k=0; k<fjc; k++) B_YM[k]=fjc+k;
				}
			}
			if (BC[4]=="surface") {
				if (fjc==1) BYM=MY+1; else {
					for (int k=0; k<fjc; k++) B_YM[k]=MY+fjc+k;
				}
			}

			//Fall through
		case 1:
			if (BC[0]=="mirror") {
				if (fjc==1) BX1=1; else {
					for (int k=0; k<fjc; k++) B_X1[k]=2*fjc-1-k;
				}
			}
			if (BC[0]=="periodic") {
				if (fjc==1) BX1=MX;else {
					for (int k=0; k<fjc; k++) B_X1[k]=MX+fjc-k-1;
				}
			}
			if (BC[0]=="surface") {
				if (fjc==1) BX1=0; else {
					for (int k=0; k<fjc; k++) B_X1[k]=k;
				}
			}

			if (BC[3]=="mirror") {
				if (fjc==1) BXM=MX;else {
					for (int k=0; k<fjc; k++) B_XM[k]=MX+fjc-k-1;
				}
			}
			if (BC[3]=="periodic") {
				if (fjc==1) BXM=1; else {
					for (int k=0; k<fjc; k++) B_XM[k]=fjc+k;
				}
			}
			if (BC[3]=="surface") {
				if (fjc==1) BXM=MX+1; else {
					for (int k=0; k<fjc; k++) B_XM[k]=MX+fjc+k;
				}
			}

	}
	if (fcc_sites) {
		fcc_lambda_1=(Real*)malloc(M*sizeof(Real)); Zero(fcc_lambda_1,M);
		fcc_lambda1=(Real*)malloc(M*sizeof(Real)); Zero(fcc_lambda1,M);
		fcc_lambda0=(Real*)malloc(M*sizeof(Real)); Zero(fcc_lambda0,M);
	}

	if (fjc==1) {
		if (gradients<3) {
		L=(Real*)malloc(M*sizeof(Real)); Zero(L,M);
		lambda_1=(Real*)malloc(M*sizeof(Real)); Zero(lambda_1,M);
		lambda1=(Real*)malloc(M*sizeof(Real)); Zero(lambda1,M);
		lambda0=(Real*)malloc(M*sizeof(Real)); Zero(lambda0,M);
		}
	} else {
		L=(Real*)malloc(M*sizeof(Real)); Zero(L,M);
		LAMBDA =(Real*)malloc(FJC*M*sizeof(Real)); Zero(LAMBDA,FJC*M);
	}
	if (Markov==2) {
		l1=(Real*)malloc(M*sizeof(Real)); Zero(l1,M);
		l_1=(Real*)malloc(M*sizeof(Real));  Zero(l_1,M);
		l11=(Real*)malloc(M*sizeof(Real)); Zero(l11,M);
		l_11=(Real*)malloc(M*sizeof(Real)); Zero(l_11,M);
		H=(Real*)malloc(M*sizeof(Real));
	}


#ifdef CUDA
	X=(Real*)AllOnDev(M);
#else
	X=(Real*)malloc(M*sizeof(Real));
#endif
	ComputeLambdas();
}



Lattice::~Lattice() {
if (debug) cout <<"lattice destructor " << endl;
	DeAllocateMemory();

}


int Lattice::P(int x, int y, int z) {
	return (x+fjc-1)*JX + (y+fjc-1)*JY +(z+fjc-1);
}
int Lattice::P(int x, int y) {
	return (x+fjc-1)*JX +(y+fjc-1);
}

int Lattice::P(int x) {
	return x+fjc-1;
}


bool Lattice::PutSub_box(int mx_, int my_, int mz_,int n_box_) {
	bool success = true;
	if (mx_<1 || my_<1 || mz_<1 || mx_>MX || my_>MY || mz_>MZ) {cout <<"subbox size out of bound: mx= " << mx_ << " my = " << my_ << " mz = " << mz_ << ", while MX = " << MX << " MY = " << MY << " MZ = " << MZ  << endl; success=false; }
	mx.push_back(mx_); my.push_back(my_); mz.push_back(mz_);
	m.push_back((mx_+2)*(my_+2)*(mz_+2));
	jx.push_back((mx_+2)*(my_+2)); jy.push_back(my_+2);
	n_box.push_back(n_box_);
	return success;
}

bool Lattice::CheckInput(int start) {
if (debug) cout <<"CheckInput in lattice " << endl;
	bool success=true;
	mx.push_back(0); my.push_back(0); mz.push_back(0); jx.push_back(0); jy.push_back(0); m.push_back(0); n_box.push_back(0);
	string Value;

	success = In[0]->CheckParameters("lat",name,start,KEYS,PARAMETERS,VALUES);
	if (!success) return success;

	vector<string> options;

		FJC=3;	fjc=1;
		if (success && GetValue("FJC_choices").length()>0) {
			if (!In[0]->Get_int(GetValue("FJC_choices"),FJC,"FJC_choices can adopt only few integer values: 3 + i*2, with i = 0, 1, 2, 3, ..."))
				success=false;
			else if ((FJC-3) %2 != 0) {
				cout << "FJC_choices can adopt only few integer values: 3 + i*2, with i = 0, 1, 2, 3, ...." <<endl;
				success=false;
			}

			fjc=(FJC-3)/2+1;
		}

		if (success && GetValue("b/l").length()>0) {
			int fjc_new=1;
			if (!In[0]->Get_int(GetValue("b/l"),fjc_new,"b/l can adopt only few integer values: 1, 2, 3, ..."))
				success=false;
			else {
				if (fjc_new <1 ) {
					cout << "b/l should be a positive integer: 1, 2, 3, ...." <<endl;
					success=false;
				}
				if (GetValue("FJC_choices").length()>0 && fjc_new !=fjc) {
					cout <<"You have set both 'FJC_choices' and 'b/l', but their values are not consistent with each other."<<endl;
					if (fjc_new<fjc && fjc_new >0) {
							cout <<"The value of 'b/l' is used, and that of FJC_choices is rejected." << endl;
					} else {
							if (fjc_new<1) success=true;
							cout <<"The value of 'FJC_choices' is used, and that of b/l is rejected." << endl;
							fjc_new=fjc;
					}
				}
			}
			fjc=fjc_new;
			FJC=3+(fjc-1)*2;
		}

		bond_length=0;
		if (GetValue("bondlength").size()>0) {
			bond_length =  In[0]->Get_Real(GetValue("bondlength"),5e-10);
			if (bond_length < 1e-11 || bond_length > 1e-8) {cout <<" bondlength out of range 1e-11..1e-8 " << endl; success=false;}
		}
		bond_length/=fjc;

		string lat_type;
		//lat_type="simple_cubic"; lattice_type=simple_cubic; lambda=1.0/6.0; Z=6;
		lattice_type=simple_cubic;
		options.push_back("simple_cubic"); options.push_back("hexagonal");
		Value=GetValue("lattice_type");
		if (Value.length()>0) {
			if (!In[0]->Get_string(Value,lat_type,options,"Input for 'lattice_type' not recognized. 'simple_cubic' or 'hexagonal'.")) success = false; else {
				if (lat_type == "simple_cubic") {lattice_type=simple_cubic; lambda=1.0/6.0; Z=6;}
				if (lat_type == "hexagonal") {lattice_type=hexagonal; lambda=1.0/4.0; Z=4;}
			}
		} else {
			success=false; cout <<"Namics can not run without input for 'lattice_type'" << endl;
		}

		offset_first_layer =0;
		gradients=1;
		gradients=In[0]->Get_int(GetValue("gradients"),1);
		if (gradients<0||gradients>3) {cout << "value of gradients out of bounds 1..3; default value '1' is used instead " << endl; gradients=1;}
		switch(gradients) {
			case 1:
				MX = In[0]->Get_int(GetValue("n_layers"),-123);

				if (MX==-123) {success=false; cout <<"In 'lat' the parameter 'n_layers' is required. Problem terminated" << endl;}
				else {
					if (MX<0 || MX >1e6) {
						success = false;
						cout <<"n_layers out of bounds, currently: 0..1e6; Problem terminated" << endl;
					}
				}
				MX=fjc*(MX);
				options.clear();
				options.push_back("spherical");
				options.push_back("cylindrical");
				options.push_back("flat");options.push_back("planar");

				if (GetValue("geometry").size()>0) {
					if (!In[0]->Get_string(GetValue("geometry"),geometry,options,"In lattice input for 'geometry' not recognized."))
						success=false;
				} else geometry = "planar";
				if (geometry=="flat") geometry="planar";
				if (geometry!="planar") {
					offset_first_layer=In[0]->Get_Real(GetValue("offset_first_layer"),0);
					if (offset_first_layer<0) {
						cout <<"value of 'offset_first_layer' can not be negative. Value ignored. " << endl;
						offset_first_layer=0;
					}
				}
				offset_first_layer *=fjc;

				options.clear();
				options.push_back("mirror");
				options.push_back("surface");
				if (GetValue("lowerbound_x").size()>0) {success=false; cout << "lowerbound_x is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("lowerbound_y").size()>0) {success=false; cout << "lowerbound_y is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("lowerbound_z").size()>0) {success=false; cout << "lowerbound_z is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("upperbound_x").size()>0) {success=false; cout << "upperbound_x is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("upperbound_y").size()>0) {success=false; cout << "upperbound_y is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("upperbound_z").size()>0) {success=false; cout << "upperbound_z is not allowed in 1-gradient calculations" << endl;}


				if (GetValue("lowerbound").size()==0) BC[0]="mirror";
				else if (!In[0]->Get_string(GetValue("lowerbound"),BC[0],options,"For 'lowerbound' boundary condition not recognized. ")) success=false;

				if (GetValue("upperbound").size()==0) BC[3]="mirror";
				else if (!In[0]->Get_string(GetValue("upperbound"),BC[3],options,"For 'upperbound' boundary condition not recognized."))
					success = false;

				break;
			case 2:
				if (GetValue("upperbound").size()>0) {success=false; cout << "upperbound is only allowed in 1-gradient calculations" << endl;}
				if (GetValue("lowerbound").size()>0) {success=false; cout << "lowerbound is only allowed in 1-gradient calculations" << endl;}

				MX = In[0]->Get_int(GetValue("n_layers_x"),-123);
				if (MX==-123) {
					success=false;
					cout <<"In 'lat' the parameter 'n_layers_x' is required. Problem terminated" << endl;
				}
				else {
					if (MX<0 || MX >1e6) {
						success = false;
						cout <<"n_layers_x out of bounds, currently: 0.. 1e6; Problem terminated" << endl;
					}
				}
				MX=fjc*(MX);
				MY = In[0]->Get_int(GetValue("n_layers_y"),-123);
				if (MY==-123) {
					success=false;
					cout <<"In 'lat' the parameter 'n_layers_y' is required. Problem terminated" << endl;
				}
				else {
					if (MY<0 || MY >1e6) {
						success = false;
						cout <<"n_layers_y out of bounds, currently: 0.. 1e6; Problem terminated" << endl;
					}
				}
				MY=fjc*(MY);
				options.clear();
				options.push_back("cylindrical");
				options.push_back("flat");options.push_back("planar");
				if (GetValue("geometry").size()>0) {
					if (!In[0]->Get_string(GetValue("geometry"),geometry,options,"In lattice input for 'geometry' not recognized."))
						success=false;
				} else geometry = "planar";
				if (geometry=="flat") geometry="planar";
				if (geometry=="planar") {volume = MX*MY;}

				if (geometry!="planar") {
					offset_first_layer=In[0]->Get_Real(GetValue("offset_first_layer"),0);
					if (offset_first_layer<0) {
						cout <<"value of 'offset_first_layer' can not be negative. Value ignored. " << endl;
						offset_first_layer=0;
					}
				}

				options.clear();
				options.push_back("mirror");
				//options.push_back("surface");
				if (geometry=="planar")
					options.push_back("periodic");
				if (GetValue("lowerbound_z").size()>0) {
					cout << "lowerbound_z is not allowed in 2-gradient calculations" << endl;
					success=false;
				}
				if (GetValue("upperbound_z").size()>0) {
					cout << "upperbound_z is not allowed in 2-gradient calculations" << endl;
					success=false;
				}

				if (GetValue("lowerbound_x").size()==0) BC[0]="mirror";
				else if (!In[0]->Get_string(GetValue("lowerbound_x"),BC[0],options,"for 'lowerbound_x' boundary condition not recognized. Put 'mirror' or 'periodic' and put 'surface' inside system. ")) success=false;

				if (GetValue("upperbound").size()==0) BC[3]="mirror";
				else if (!In[0]->Get_string(GetValue("upperbound"),BC[3],options,"for 'upperbound_x' boundary condition not recognized. Put 'mirror' or 'periodic' and put 'surface' inside system.")) success = false;

				if (GetValue("lowerbound_y").size()==0) BC[1]="mirror";
				else if (!In[0]->Get_string(GetValue("lowerbound_y"),BC[1],options,"for 'lowerbound_y' boundary condition not recognized. Put 'mirror' or 'periodic' and put 'surface' inside system.")) success=false;

				if (GetValue("upperbound_y").size()==0) BC[4]="mirror";
				else if (!In[0]->Get_string(GetValue("upperbound_y"),BC[4],options,"for 'upperbound_Y' boundary condition not recognized. Put 'mirror' or 'periodic' and put 'surface' inside system.")) success = false;


				if (BC[0]=="periodic" || BC[1]=="periodic") {
					if (BC[0]!=BC[1]) {
						success=false;
						cout <<"For boundaries in x-direction: 'periodic' BC  should be set to upper and lower bounds " << endl;
					}
				}
				if (BC[2]=="periodic" || BC[3]=="periodic") {
					if (BC[2]!=BC[3]) {success=false;  cout <<"For boundaries in y-direction: 'periodic' BC should be set to upper and lower bounds " << endl;}
				}
				break;
			case 3:
				if (GetValue("upperbound").size()>0) {success=false; cout << "upperbound is only allowed in 1-gradient calculations" << endl;}
				if (GetValue("lowerbound").size()>0) {success=false; cout << "lowerbound is only allowed in 1-gradient calculations" << endl;}

				if (!In[0]->Get_int(GetValue("n_layers_x"),MX,1,1e6,"In 'lat' the parameter 'n_layers_x' is required"))
					success=false;
				if (!In[0]->Get_int(GetValue("n_layers_y"),MY,1,1e6,"In 'lat' the parameter 'n_layers_y' is required"))
					success=false;
				if (!In[0]->Get_int(GetValue("n_layers_z"),MZ,1,1e6,"In 'lat' the parameter 'n_layers_z' is required"))
					success=false;
				MX=fjc*(MX);
				MY=fjc*(MY);
				MZ=fjc*(MZ);
				options.clear();
				options.push_back("mirror"); //options.push_back("mirror_2");
				options.push_back("periodic");
				//options.push_back("surface");
				//options.push_back("shifted_mirror");


				if (GetValue("lowerbound_x").size()==0) BC[0]="mirror";
				else if (!In[0]->Get_string(GetValue("lowerbound_x"),BC[0],options,"for 'lowerbound_x' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success=false;

				if (GetValue("upperbound").size()==0) BC[3]="mirror";
				else if (!In[0]->Get_string(GetValue("upperbound_x"),BC[3],options,"for 'upperbound_x' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success = false;

				if (GetValue("lowerbound_y").size()==0) BC[1]="mirror";
				else if (!In[0]->Get_string(GetValue("lowerbound_y"),BC[1],options,"for 'lowerbound_y' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success=false;

				if (GetValue("upperbound_y").size()==0) BC[4]="mirror";
				else if (!In[0]->Get_string(GetValue("upperbound_y"),BC[4],options,"for 'upperbound_y' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success = false;

				if (GetValue("lowerbound_z").size()==0) BC[2]="mirror";
				else if (!In[0]->Get_string(GetValue("lowerbound_z"),BC[2],options,"for 'lowerbound_z' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success=false;

				if (GetValue("upperbound_z").size()==0) BC[5]="mirror";
				else if (!In[0]->Get_string(GetValue("upperbound_z"),BC[5],options,"for 'upperbound_z' boundary condition not recognized. Put 'mirror' or 'periodic' and put surface inside system. ")) success = false;

				if (BC[2]=="periodic" || BC[3]=="periodic") {
					if (BC[2] != BC[3]) {
					cout <<"In y-direction the boundary conditions do not match:" + BC[2] << " and " <<  BC[3] << endl;
					success=false;
					}
				}
				if (BC[4]=="periodic" || BC[5]=="periodic") {
					if (BC[4] != BC[5]) {
						cout <<"In z-direction the boundary conditions do not match:" + BC[4] << " and " <<  BC[5] << endl;
						success=false;
					}
				}

				break;
			default:
				cout << "gradients out of bounds " << endl;
				break;
		}
		if ((fjc>1) && (lattice_type ==hexagonal)) {success = false; cout << "For FJC-choices >3, we need lattice_type = 'hexagonal'." << endl; }
		if (gradients ==2 && fjc>2) {success = false; cout <<" When gradients is 2, FJC-choices are limited to 5 " << endl; }
		if (gradients ==3 && fjc>1) {success = false; cout <<" When gradients is 3, FJC-choices are limited to 3 " << endl; }

		if (GetValue("ignore_site_fraction").length()>0) {
			ignore_sites=In[0]->Get_bool(GetValue("ignore_sites"),false);
			if (!ignore_sites) cout <<"ignore_site_fraction is set to false. Full site fractions computed. " << endl;
		}

		if (GetValue("fcc_site_fraction").length()>0) {
			fcc_sites=In[0]->Get_bool(GetValue("fcc_site_fraction"),false);
			if (!fcc_sites) cout <<"fcc_site_fraction is set to false. Full site fractions computed. " << endl;
		}

		if (fcc_sites&&ignore_sites) {
			cout <<"can't combine 'fcc_site_fraction' with 'ignore_site_fraction'" <<endl; success=false;
		}

		if (GetValue("stencil_full").length()>0) {
			stencil_full=In[0]->Get_bool(GetValue("stencil_full"),false);
			if (gradients<3 && !stencil_full) cout << "In calculations with 'gradients' less than 3, stencil_full is set to true. " << endl;
		}
		//Initialize system size and indexing
		PutM();
		//if (lattice_type==simple_cubic) {
		//	lambda=1.0/6.0; //l0=4.0*l1;
		//} else {
		//	lambda=1.0/4.0; //l0=2.0*l1;
		//}
	Markov=1;
	//Markov=In[0]->Get_int(GetValue("Markov"),1);
	//if (Markov<1 || Markov>2) {
	//	cout <<" Integer value for 'Markov' is by default 1 and may be set to 2 for some mol_types and fjc-choices only. Markov value out of bounds. Proceed with caution. " << endl; success = false;
	//}
	k_stiff=0; //default value if in mol there is no k_stiff
	if (GetValue("k_stiff").size()>0) {
		k_stiff=In[0]->Get_Real(GetValue("k_stiff"),0);
		if (k_stiff<0 || k_stiff>10) {
			success =false;
			cout <<" Real value for 'k_stiff' out of bounds (0 < k_stiff < 10). " << endl;
			cout <<" For Markov == 2: u_bend (theta) = 0.5 k_stiff theta^2, where 'theta' is angle for bond direction deviating from the straight direction. " <<endl;
			cout <<" You may interpret 'k_stiff' as the molecular 'persistence length' " << endl;
			cout <<" k_stiff is a 'default value'. Use molecular specific values to overrule the default when appropriate (future implementation....) " << endl;
		}
		if (fjc>1 ) {
			success=false;
			cout <<" Work in progress.... Currently, Markov == 2 is only expected to work for fjc_choices < 5 " << endl;
		}
	}

	return success;
}

bool Lattice::PutVarInfo(string Var_type_, string Var_target_, Real Var_target_value_){
	bool success=true;
	Var_target = -1;
	Var_type=Var_type_;
	if (Var_type !="scan") {success=false; cout <<"Var type is not 'scan' and therefore var info rejected " << endl; }
	if (Var_target_=="n_layers") {
		Var_target=0;VarInitValue=MX;
		if (gradients>1) {success=false; cout <<"Var target 'n_layers' rejected because the number of gradients >1 " << endl; }
	}
	if (Var_target_=="n_layers_x") {
		Var_target=1; VarInitValue=MX;
		if (gradients==0) {success=false; cout <<"Var target 'n_layers_x' rejected because the number of gradients =1 " << endl; }
	}
	if (Var_target_=="n_layers_y") {
		Var_target=2; VarInitValue=MY;
		if (gradients==0) {success=false; cout <<"Var target 'n_layers_y' rejected because the number of gradients =1 " << endl; }
	}
	if (Var_target_=="n_layers_z") {
		Var_target=3; VarInitValue=MZ;
		if (gradients==0) {success=false; cout <<"Var target 'n_layers_z' rejected because the number of gradients =1 " << endl; }
	}
	if (Var_target_=="offset_firstlayer") {
		Var_target=4;VarInitValue=offset_first_layer;
		if (gradients>1) {success=false; cout <<"Var target 'offset_firstlayer' rejected because the number of gradients >1 " << endl; }
		if (geometry == "planar") {success=false; cout <<"Var target 'offset_firstlayer' rejected because 'geometry' is 'planar' " << endl; }
	}
	if (Var_target_=="offset_firstlayer_x") {
		Var_target=5; VarInitValue=offset_first_layer;
		if (gradients!=2) {success=false; cout <<"Var target 'offset_firstlayer_x' rejected because the number of gradients !=2 " << endl; }
		if (geometry != "cylindrical") {success=false; cout <<"Var target 'offset_firstlayer_x' rejected because 'geometry' not 'cylindrical' " << endl; }
	}
	return success;
}

bool Lattice::UpdateVarInfo(int step_nr) {
	bool success=true;
	switch(Var_target) {
		case 0:
			MX=VarInitValue+step_nr*Var_step;
			break;
		case 1:
			MX=VarInitValue+step_nr*Var_step;
			break;
		case 2:
			MY=VarInitValue+step_nr*Var_step;
			break;
		case 3:
			MZ=VarInitValue+step_nr*Var_step;
			break;
		case 4:
			offset_first_layer=VarInitValue+step_nr*Var_step;
			break;
		case 5:
			offset_first_layer=VarInitValue+step_nr*Var_step;
			break;
		default:
			cout <<"program error in UpdateVarInfo "<<endl;
			break;
	}
	return success;
}

int Lattice::PutVarScan(int step, int end_value) {
	int num_of_steps=-1;
	Var_step=step; if (step==0) cout <<"In var scan : of lattice variable, the value of step can not be negative" << endl;
	Var_end_value=end_value;
	if (Var_end_value <0) cout <<"In var: scan : of lattice variable, it is impossible to have a negative value for the end_value" << endl;
	if (step!=0) {
		num_of_steps = (end_value-VarInitValue)/step;
		if (num_of_steps<0) cout <<"in var : scan : of lattice variable, end_value and step are not consistent, try changing sign of step " << endl;
	}
	return num_of_steps;
}

bool Lattice::ResetInitValue() {
	bool success=true;
	switch(Var_target) {
		case 0:
			MX=VarInitValue;
			break;
		case 1:
			MX=VarInitValue;
			break;
		case 2:
			MY=VarInitValue;
			break;
		case 3:
			MZ=VarInitValue;
			break;
		case 4:
			offset_first_layer=VarInitValue;
			break;
		case 5:
			offset_first_layer=VarInitValue;
			break;
		default:
			cout <<"program error in ResetInitValue "<<endl;
			break;
	}
	return success;
}

void Lattice::PutParameter(string new_param) {
if (debug) cout <<"PutParameters in lattice " << endl;
	KEYS.push_back(new_param);
}

string Lattice::GetValue(string parameter){
if (debug) cout << "GetValue in lattice " << endl;
	int i=0;
	int length = PARAMETERS.size();
	while (i<length) {
		if (parameter==PARAMETERS[i]) {
			return VALUES[i];
		}
		i++;
	}
	return "" ;
}

Real Lattice::GetValue(Real* X,string s){ //need a check to find out if s contains 3 integers separated by ','
if (debug) cout << "GetValue in lattice " << endl;
if (X==NULL) cout << "pointer X is zero" << endl;
	int x=0,y=0,z=0;
	vector<string> sub;
	In[0]->split(s,',',sub);
	switch(gradients) {
		case 1:
			if (sub.size()==1) {
				x=In[0]->Get_int(sub[0],x);
				if (x<0||x>MX+1) {
					cout <<"Requested postition in 'kal' output out of bounds." << endl;
					return 0;
				} else return X[x];
			} else cout <<"Request for profile output does not contain the expected coordinate in 'kal' output" << endl;
			break;
		case 2:
			if (sub.size()==2) {
				x=In[0]->Get_int(sub[0],x);
				y=In[0]->Get_int(sub[1],y);
				if (x<0||x>MX+1||y<0||y>MY+1) {
					cout <<"Requested postition in 'kal' output out of bounds." << endl;
					return 0;
				} else return X[JX*x+y];
			} else cout <<"Request for profile output does not contain the expected coordinate in 'kal' output" << endl;
			break;
		case 3:
			if (sub.size()>2) {
				x=In[0]->Get_int(sub[0],x);
				y=In[0]->Get_int(sub[1],y);
				z=In[0]->Get_int(sub[2],y);
				if (x<0||x>MX+1||y<0||y>MY+1||z<0||z>MZ+1) {
					cout <<"Requested postition in 'kal' output out of bounds." << endl;
					return 0;
				} else return X[JX*x+JY*y+z];
			} else  cout <<"Request for profile output does not contain the coordinate in 'kal' output" << endl;
			return 0;
			break;
		default:
			break;
	}
	return 0;
}

bool Lattice:: PutMask(int* MASK,vector<int>px,vector<int>py,vector<int>pz,int R){
	bool success=true;
	int length =px.size();
	int X,Y,Z;

	for (int i =0; i<length; i++) {
		int xx,yy,zz;
		xx=px[i]; yy=py[i]; zz=pz[i];
		for (int x=xx-R; x<xx+R+1; x++)
		for (int y=yy-R; y<yy+R+1; y++)
		for (int z=zz-R; z<zz+R+1; z++) {
			if ((xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z) <=R*R) {
				X=x; Y=y; Z=z;
				if (x<1) X+=MX;
				if (y<1) Y+=MY;
				if (z<1) Z+=MZ;
				if (x>MX) X-=MX;
				if (y>MY) Y-=MY;
				if (z>MZ) Z-=MZ;
				MASK[P(X,Y,Z)]++;
			}
		}
		for (int x=1; x<MX; x++)
		for (int y=1; y<MY; y++)
		for (int z=1; z<MZ; z++)
		if (MASK[P(x,y,z)]>1) success=false;
	}

	return success;
}

bool Lattice::PrepareForCalculations(void) {
if (debug) cout <<"PrepareForCalculations in lattice" << endl;
	bool success=true;
	return success;
}

void Lattice::push(string s, Real X) {
if (debug) cout <<"push (Real) in lattice " << endl;
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void Lattice::push(string s, int X) {
if (debug) cout <<"push (int) in lattice " << endl;
	ints.push_back(s);
	ints_value.push_back(X);
}
void Lattice::push(string s, bool X) {
if (debug) cout <<"push (bool) in lattice " << endl;
	bools.push_back(s);
	bools_value.push_back(X);
}
void Lattice::push(string s, string X) {
if (debug) cout <<"push (string) in lattice " << endl;
	strings.push_back(s);
	strings_value.push_back(X);
}

Real* Lattice::GetPointer(string s,int &SIZE) {
if (debug) cout <<"GetPointer for lattice " + name << endl;
	vector<string> sub;
	SIZE=M;
	In[0]->split(s,';',sub);
	if (sub[0]=="profile" && sub[1]=="0") return L;
	if (sub[0]=="vector") {}
	return NULL;
}

int* Lattice::GetPointerInt(string s,int &SIZE) {
if (debug) cout <<"GetPointerInt for lattice " + name << endl;
	vector<string> sub;
	SIZE=M;
	In[0]->split(s,';',sub);
	if (sub[0]=="array"){//get with sub[1] the number and put the pointer to integer array in return.
	}

	return NULL;
}

void Lattice::PushOutput() {
if (debug) cout <<"PushOutput in lat " << endl;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	string mirror="mirror";
	string periodic="periodic";
	//string surface="surface";
	push("gradients",gradients);
	if (offset_first_layer>0) push("offset_first_layer",offset_first_layer);
	push("volume",volume);
	push("accessible volume",Accesible_volume);
	if (lattice_type == simple_cubic) push("lattice_type",simple_cubic); else push("lattice_type",hexagonal);
	push("bond_length",bond_length);
	push("FJC_choices",FJC);
	string s="profile;0"; push("L",s);

	switch (gradients) {
		case 3:
			push("n_layers_z",MZ);
			//if (BZ1==1) push("lowerbound_z",mirror); //should be fixed for fjc>1
			//if (BZM==MZ-1) push("upperbound_z",mirror);

			//if (BX1==MX) push("lowerbound_x",periodic);
			//if (BXM==1) push("upperbound_x",periodic);
			//if (BY1==MY) push("lowerbound_y",periodic);
			//if (BYM==1) push("upperbound_y",periodic);
			//if (BZ1==MZ) push("lowerbound_z",periodic);
			//if (BZM==1) push("upperbound_z",periodic);
			// Fall through
		case 2:
			push("n_layers_y",MY);
			//if (BY1==1) push("lowerbound_x",mirror);
			//if (BYM==MY-1) push("upperbound_x",mirror);
			// Fall through
		case 1:
			push("n_layers",MX);
			//if (BX1==1) push("lowerbound",mirror);
			//if (BXM==MX-1) push("upperbound",mirror);
			break;
		default:
			break;
	}
}


int Lattice::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
if (debug) cout <<"GetValue (long)  in lattice " << endl;

	for ( size_t i = 0 ; i<ints.size() ; ++i)
		if (prop==ints[i]) {
			int_result=ints_value[i];
			return 1;
		}

	for ( size_t i = 0 ; i<Reals.size() ; ++i)
		if (prop==Reals[i]) {
			Real_result=Reals_value[i];
			return 2;
		}

	for ( size_t i = 0 ; i<bools.size() ; ++i)
		if (prop==bools[i]) {
			if (bools_value[i]) string_result="true";
			else string_result="false";
			return 3;
		}

	for ( size_t i = 0 ; i<strings.size() ; ++i)
		if (prop==strings[i]) {
			string_result=strings_value[i];
			return 3;
		}

	return 0;
}




void Lattice::DistributeG1(Real *G1, Real *g1, int* Bx, int* By, int* Bz, int n_box) {
	int k=sub_box_on;
	tools::DistributeG1(G1, g1, Bx, By, Bz, M, m[k], n_box, mx[k], my[k], mz[k], MX, MY, MZ, jx[k], jy[k], JX, JY);
}

void Lattice::CollectPhi(Real* phi, Real* GN, Real* rho, int* Bx, int* By, int* Bz, int n_box) {
	int k=sub_box_on;
	tools::CollectPhi(phi, GN, rho, Bx, By, Bz, M, m[k], n_box, mx[k], my[k], mz[k], MX, MY, MZ, jx[k], jy[k], JX, JY);
}

void Lattice::ComputeGN(Real* GN, Real* Gg_f, int* H_Bx, int* H_By, int* H_Bz, int* H_Px2, int* H_Py2, int* H_Pz2, int N, int n_box) {
	int k=sub_box_on;
	for (int p=0; p<n_box; p++) Cp(GN+p,Gg_f+n_box*m[k]*N +p*m[k]+ jx[k]*(H_Px2[p]-H_Bx[p])+jy[k]*(H_Py2[p]-H_By[p])+(H_Pz2[p]-H_Bz[p]),1);

#ifdef CUDA //this transfer can go away when all is on GPU.
	//TransferDataToHost(H_GN,GN,n_box);
#endif
}


bool Lattice::ReadGuess(string filename, Real *x ,string &method, vector<string> &monlist, vector<string> &statelist, bool &charged, int &mx, int &my, int &mz, int &fjc, int readx) {
if (debug) cout <<"ReadGuess in output" << endl;
	bool success=true;
	ifstream in_file;
	string s_charge;
	int num_1,num_2;
	string name;
	in_file.open(filename.c_str());
	if (in_file.is_open()) {
		while (in_file && readx>-1) {
			in_file>>method;
			in_file>>mx>>my>>mz>>fjc;
			in_file>>s_charge;
			if (s_charge=="true") charged=true; else charged=false;
			in_file>>num_1;
			for (int i=0; i<num_1; i++) {
				in_file >> name;
				monlist.push_back(name);
			}
			in_file>>num_2;
			for (int i=0; i<num_2; i++) {
				in_file>> name;
				statelist.push_back(name);
			}
			if (readx==0) readx=-2; else {
				int m;
				if (my==0) {m=(mx+2*fjc);} else {if (mz==0) m=(mx+2*fjc)*(my+2*fjc); else {m=(mx+2*fjc)*(my+2*fjc)*(mz+2*fjc);}}
				int iv=(num_1+num_2)*m;
				if (charged) iv +=m;
				for (int i=0; i<iv; i++) {
					in_file >> x[i];
					//cout <<"R i " << i << " " << x[i] << endl;
				}
				readx=-3;
			}
		}
		in_file.close();
	} else {
		cout <<"inputfile " << filename << "is not found. Read guess for initial guess failed" << endl;
		success=false;
	}
	return success;
}

bool Lattice::StoreGuess(string Filename,Real *x,string method, vector<string> monlist,vector<string>statelist, bool charged, int start) {
if (debug) cout <<"StoreGuess in output" << endl;
	bool success=true;
	string s;
	int mon_length = monlist.size();
	int state_length = statelist.size();
	string filename;
	string outfilename;
	vector<string> sub;
	if (Filename == "") {
		outfilename=In[0]->name;
		In[0]->split(outfilename,'.',sub);
		char numc[2];
       	sprintf(numc,"%d",start);
		filename=sub[0].append("_").append(numc).append(".").append("outiv");
	} else {
		outfilename = Filename;
		In[0]->split(outfilename,'.',sub);
		char numc[2];
       	sprintf(numc,"%d",start);
		filename=sub[0].append("_").append(numc).append(".").append(sub[1]);
	}
	FILE *fp;
	fp=fopen(filename.c_str(),"w");
	fprintf(fp,"%s\n",method.c_str());
	fprintf(fp," %i\t%i\t%i\t%i\n" ,MX,MY,MZ, fjc);
	if (charged) fprintf(fp,"%s\n" ,"true"); else  fprintf(fp,"%s\n" ,"false");
	fprintf(fp,"%i\n",mon_length);
	for (int i=0; i<mon_length; i++) fprintf(fp,"%s\n",monlist[i].c_str());
	fprintf(fp,"%i\n",state_length);
	for (int i=0; i<state_length; i++) fprintf(fp,"%s\n",statelist[i].c_str());
	int iv=(mon_length+state_length)*M;
	if (charged) iv +=M;
	for (int i=0; i<iv; i++) fprintf(fp,"%e\n",x[i]);
	fclose(fp);
	return success;
}

bool Lattice::GenerateGuess(Real* x, string CalculationType, string GuessType, Real A_value, Real B_value) {
if (debug) cout <<"GenerateGuess in lat " << endl;
//GuessType: lamellae,Im3m,FCC,BCC,HEX,gyroid,Real_gyroid,Real_diamond,perforated_lamellae
//CalculationType: micro_emulsion,micro_phasesegregation
	bool success = true;
	int i,j,k;
	if (fjc>1) cout <<"GenerateGuess is not yet prepared to work for FJC-choices > 3 " << endl;
	switch (gradients) {
		case 1:
			if (GuessType=="lamellae") {
				for (i=0; i<MX+2; i++) {
					if (i<MX/2+1) {
						x[i]= A_value;
						x[i+M]= B_value;
					} else {
						x[i]=0;
						x[i+M]=0;
					}
				}
			}
			break;
		case 2:
			if (GuessType=="lamellae") {
				for (i=0; i<MX+2; i++)
					for (j=0; j<MY+2; j++) {
					if (j<MY/2+1) {
						x[JX*i+j]= A_value;
						x[JX*i+j+M]= B_value;
					} else {
						x[JX*i+j]=0;
						x[JX*i+j+M]=0;
					}
				}
			}
			break;
		case 3:
			if (GuessType=="lamellae") {
				for (i=0; i<MX+2; i++) for (j=0; j<MY+2; j++) for (k=0; k<MZ+2; k++) {
					if (k<MZ/2+1) {
						x[JX*i+JY*j+k]= A_value;
						x[JX*i+JY*j+k+M]= B_value;
					} else {
						x[JX*i+JY*j+k]=0;
						x[JX*i+JY*j+k+M]=0;
					}
				}
			}
			break;
		default:
			break;
	}
	return success;
}


bool Lattice::GuessVar(Real* x, Real theta,string GuessType, Real A_value, Real B_value){
if (debug) cout << "GuessVar in Lattice " << endl;
	bool success = true;
	int i;
	if (fjc>1) cout <<"GuessVar is not yet prepared to work for FJC-choices > 3 " << endl;
	int height = theta/1;
	switch (gradients) {
		case 1:
			if (GuessType=="lamellae"){
				for (i=0; i<MX+2; i++) {
                                        if (i<height+1) {
                                                x[i]= A_value;
                                                x[i+M]= B_value;
                                        } else {
                                                x[i]=0;
                                                x[i+M]=0;
                                        }

				}
			}
			break;
		case 2: cout << "var is not implemented in 2 or 3 gradient problems yet" << endl; success=false; break;
		case 3: cout << "var is not implemented in 2 or 3 gradient problems yet" << endl; success=false; break;

	}
	return success;
}



