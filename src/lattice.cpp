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
 	KEYS.push_back("bondlength");  KEYS.push_back("lattice_type");
	KEYS.push_back("lambda");
	KEYS.push_back("Z");
	KEYS.push_back("FJC_choices");
	sub_box_on = 0;
	all_lattice = false;
	fjc = 0;
	gradients=1;
	fjc=1;
	MX=MY=MZ=0;
}

Lattice::~Lattice() {
if (debug) cout <<"lattice destructor " << endl;
  //In this program, we will assume that the propagator will work on a simple cubic lattice.
			//Interactions will be treated either with simple cubic or FCC 'lambda's.
			DeAllocateMemory();
}

void Lattice::DeAllocateMemory(void) {
if (debug) cout <<"DeAllocateMemory in lattice " << endl;
	if (all_lattice) {
		if (fjc==1) {
			free(lambda_1);
			free(lambda1);
			free(lambda0);
			free(L);
		} else {
			free(L); free(LAMBDA);
		} 
	}
#ifdef CUDA
//	if (gradients==3) X=(Real*)AllOnDev(M);
//	if (gradients==3) cudaFree(X);
	all_lattice=false;
#endif
}

void Lattice::AllocateMemory(void) {
if (debug) cout <<"AllocateMemory in lattice " << endl;
	Real r,VL,LS;
	Real rlow,rhigh;
	Real one_over_lambda_min_two=1.0/lambda-2.0;
	int i {0}; int j {0}; int k {0};
	PutM();
	DeAllocateMemory();
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

	switch (gradients) {
		case 1:
			if (fjc==1) {
				if (geometry=="planar") {for (i=1; i<MX+1; i++) L[i]=1; }
				if (geometry=="cylindrical") {
					for (i=1; i<MX+1; i++) {
						r=offset_first_layer + i;
						L[i]=PIE*(pow(r,2)-pow(r-1,2));
						lambda1[i]=2.0*PIE*r/L[i];
						lambda_1[i]=2.0*PIE*(r-1)/L[i];
						lambda0[i]=one_over_lambda_min_two;
					}
				}
				if (geometry=="spherical") {
					for (i=1; i<MX+1; i++) {
						r=offset_first_layer + i;
						L[i]=4.0/3.0*PIE*(pow(r,3)-pow(r-1,3));
						lambda1[i]=4*PIE*pow(r,2)/L[i];
						lambda_1[i]=4*PIE*pow(r-1,2)/L[i];
						lambda0[i]=1.0/lambda-lambda1[i]-lambda_1[i];
					}
				}
			} else { //fjc != 1;
				if (geometry =="planar") {
					for (i=0; i<M; i++) {
						L[i]=1.0/fjc;
						LAMBDA[i]=1.0/(2*(FJC-1));
						LAMBDA[i+(FJC-1)*M]=1.0/(2*(FJC-1));
						//LAMBDA[i+(FJC-1)/2*M]=1.0-1.0/(FJC-1);
						LAMBDA[i+(FJC-1)/2*M]=1.0/(FJC-1);
						for (j=1; j<FJC/2; j++) {
							LAMBDA[i+j*M]=1.0/(FJC-1);
							LAMBDA[i+(FJC-j-1)*M]=1.0/(FJC-1);
							//LAMBDA[i+(FJC-1)/2*M]=LAMBDA[i+(FJC-1)/2*M]-LAMBDA[i+j*M]-LAMBDA[i+(FJC-j-1)*M];
						}
					}

				}
				if (geometry =="cylindrical") {

					for (i=fjc; i<M-fjc; i++) {
						r=offset_first_layer+1.0*(1.0*i-1.0*fjc+1.0)/fjc;
						rlow=r-0.5;
						rhigh=r+0.5;
						L[i]=PIE*(2.0*r)/fjc;
						VL=L[i]/PIE*fjc;
						if ((rlow-r)*2+r>0.0) LAMBDA[i]+=0.5/(1.0*FJC-1.0)*2.0*rlow/VL;
						if ((rhigh-r)*2+r<MX) LAMBDA[i+(FJC-1)*M]+=0.5/(1.0*FJC-1.0)*2.0*rhigh/VL;
						else {
							if (2*rhigh-r-MX>-0.001 && 2*rhigh-r-MX < 0.001) {
								LAMBDA[i+(FJC-1)*M]+= 0.5/(1.0*FJC-1.0)*2.0*rhigh/VL;
							}
							for (j=1; j<=fjc; j++) {
								if (2*rhigh-r-MX > 0.99*j/fjc && 2*rhigh-r-MX < 1.01*j/fjc) {
									LAMBDA[i+(FJC-1)*M]+= 0.5/(1.0*FJC-1.0)*2.0*(rhigh-1.0*j/fjc)/VL;
								}
							}
						}
						for (j=1; j<fjc; j++) {
							rlow += 0.5/(fjc);
							rhigh -= 0.5/(fjc);
							if ((rlow-r)*2+r>0.0) LAMBDA[i+j*M]+=1.0/(1.0*FJC-1.0)*2.0*rlow/VL;
							if ((rhigh-r)*2+r < offset_first_layer+MX) LAMBDA[i+(FJC-1-j)*M]+=1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
							else {
								if (2*rhigh-r-MX>-0.001 && 2*rhigh-r-MX < 0.001 ) {
									LAMBDA[i+(FJC-1-j)*M]+= 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
								}
								for (k=1; k<=fjc; k++) {
									if (2*rhigh-r-MX > 0.99*k/fjc && 2*rhigh-r-MX < 1.01*k/fjc) {
										LAMBDA[i+(FJC-1-j)*M]+= 1.0/(1.0*FJC-1.0)*2.0*(rhigh-1.0*k/fjc)/VL;
									}
								}
							}
						}
						LS=0;
						for (j=0; j<FJC; j++) LS+=LAMBDA[i+j*M];
						LAMBDA[i+(FJC/2)*M]+=1.0-LS;
					}
				}


				if (geometry =="spherical") {//todo : test
					for (i=fjc; i<M-fjc; i++) {
						r=offset_first_layer+1.0*(1.0*i-1.0*fjc+1.0)/fjc;
						rlow=r-0.5;
						rhigh=r+0.5;
						L[i]=PIE*4.0/3.0*(rhigh*rhigh*rhigh-rlow*rlow*rlow)/fjc;
						VL=L[i]/PIE*fjc;
						if ((rlow-r)*2+r>0.0) LAMBDA[i]+=0.5/(1.0*FJC-1.0)*4.0*rlow*rlow/VL;
						if ((rhigh-r)*2+r<MX) LAMBDA[i+(FJC-1)*M]+=0.5/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
						else {
							if (2*rhigh-r-MX>-0.001 && 2*rhigh-r-MX < 0.001) {
								LAMBDA[i+(FJC-1)*M]+= 0.5/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
							}
							for (j=1; j<=fjc; j++) {
								if (2*rhigh-r-MX > 0.99*j/fjc && 2*rhigh-r-MX < 1.01*j/fjc) {
									LAMBDA[i+(FJC-1)*M]+= 0.5/(1.0*FJC-1.0)*4.0*(rhigh-1.0*j/fjc)*(rhigh-1.0*j/fjc)/VL;
								}
							}
						}
						for (j=1; j<fjc; j++) {
							rlow += 0.5/(fjc);
							rhigh -= 0.5/(fjc);
							if ((rlow-r)*2+r>0.0) LAMBDA[i+j*M]+=1.0/(1.0*FJC-1.0)*4.0*rlow*rlow/VL;
							if ((rhigh-r)*2+r < offset_first_layer+MX) LAMBDA[i+(FJC-1-j)*M]+=1.0/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
							else {
								if (2*rhigh-r-MX>-0.001 && 2*rhigh-r-MX < 0.001 ) {
									LAMBDA[i+(FJC-1-j)*M]+= 1.0/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
								}
								for (k=1; k<=fjc; k++) {
									if (2*rhigh-r-MX > 0.99*k/fjc && 2*rhigh-r-MX < 1.01*k/fjc) {
										LAMBDA[i+(FJC-1-j)*M]+= 1.0/(1.0*FJC-1.0)*4.0*(rhigh-1.0*k/fjc)*(rhigh-1.0*k/fjc)/VL;
									}
								}
							}
						}
						LS=0;
						for (j=0; j<FJC; j++) LS+=LAMBDA[i+j*M];
						LAMBDA[i+(FJC/2)*M]+=1.0-LS;
					}
				}
			}
			break;
		case 2:
			if (geometry=="planar") {
				for (i=0; i<M; i++) L[i]=1;
			}
			if (geometry=="cylindrical") {
				for (i=1; i<MX+1; i++)
				for (j=1; j<MY+1; j++) {
					r=offset_first_layer + i;
					L[i*JX+j]=PIE*(pow(r,2)-pow(r-1,2));
					lambda1[i*JX+j]=2.0*PIE*r/L[i*JX+j];
					lambda_1[i*JX+j]=2.0*PIE*(r-1)/L[i*JX+j];
					lambda0[i*JX+j]=one_over_lambda_min_two;
				}

			}
			break;
		case 3:
			break;
		default:
			break;
	}
#ifdef CUDA
		if (gradients==3)
			X=(Real*)AllOnDev(M);
#endif
		all_lattice=(gradients<2 && geometry!="planar");
}

bool Lattice::PutM() {
	bool success=true;
	switch (gradients) {
		case 3:
			if (BC[4]=="mirror") BZ1=1;
			if (BC[4]=="periodic") BZ1=MZ;
			//if (BC[5]=="surface") BZM=MZ+1;
			if (BC[5]=="mirror") BZM=MZ;
			if (BC[5]=="periodic") BZM=1;
			//Fall through
		case 2:
			//if (BC[2]=="surface") BY1=0;
			if (BC[2]=="mirror") BY1=1;
			if (BC[2]=="periodic") BY1=MY;
			//if (BC[3]=="surface") BYM=MY+1;
			if (BC[3]=="mirror") BYM=MY;
			if (BC[3]=="periodic") BYM=1;
			//Fall through
		case 1:
			//if (BC[0]=="surface") BX1=0;
			if (BC[0]=="mirror") BX1=1;
			if (BC[0]=="periodic") BX1=MX;
			//if (BC[1]=="surface") BXM=MX+1;
			if (BC[1]=="mirror") BXM=MX;
			if (BC[1]=="periodic") BXM=1;
	}
	switch(gradients) {
		case 1:
			JX=1; JY=0; JZ=0; M=MX+2;
			if (fjc>1) M *=fjc;
			if (geometry=="planar") {volume = MX;}
			if (geometry=="spherical") {volume = 4/3*PIE*(pow(MX+offset_first_layer,3)-pow(offset_first_layer,3));}
			if (geometry=="cylindrical") {volume = PIE*(pow(MX+offset_first_layer,2)-pow(offset_first_layer,2));}
			break;
		case 2:
			if (geometry=="cylindrical")
				volume = MY*PIE*(pow(MX+offset_first_layer,2)-pow(offset_first_layer,2));
			else volume = MX*MY;
			JX=(MY+2); JY=1; JZ=0; M=(MX+2)*(MY+2);
			break;
		case 3:
			volume = MX*MY*MZ;
			JX=(MZ+2)*(MY+2); JY=(MZ+2); JZ=1; M = (MX+2)*(MY+2)*(MZ+2);
			break;
		default:
			break;
	}
	return success;
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
	bool success;
	mx.push_back(0); my.push_back(0); mz.push_back(0); jx.push_back(0); jy.push_back(0); m.push_back(0); n_box.push_back(0);
	string Value;

	success = In[0]->CheckParameters("lat",name,start,KEYS,PARAMETERS,VALUES);
	if (!success) return success;

	vector<string> options;

		bond_length=0;
		if (GetValue("bondlength").size()>0) {
			bond_length =  In[0]->Get_Real(GetValue("bondlength"),5e-10);
			if (bond_length < 1e-10 || bond_length > 1e-8) {cout <<" bondlength out of range 1e-10..1e-8 " << endl; success=false;}
		}

		lambda=0;
		lattice_type="";
		options.push_back("simple_cubic"); options.push_back("FCC"); options.push_back("hexagonal");
		Value=GetValue("lattice_type");
		if (Value.length()>0) {
			if (!In[0]->Get_string(Value,lattice_type,options,"Input for 'lattice_type' not recognized. 'Simple_cubic': lambda=1/6; Z=6; FCC: lambda=1/3; Z=3; 'hexagonal': lambda=1/4; Z=4;")) success = false; else {
				if (lattice_type == "simple_cubic") {lambda=1.0/6.0; Z=6;}
				if (lattice_type == "FCC") {lambda = 1.0/3.0; Z=3;}
				if (lattice_type == "hexagonal") {lambda=1.0/4.0; Z=4;}
			}
		} else {
			if (lambda>0) {
				if (GetValue("lambda").size()>0) cout <<"'lattice_type' has been set: 'lambda'-value is ignored" << endl;
				if (GetValue("Z").size()>0) cout <<"'lattice_type' has been set: 'Z' (coordination number) is ignored" << endl;
			} else {
				if (GetValue("lambda").size()>0) {
					lambda=In[0]->Get_Real(GetValue("lambda"),1.0/6.0);
					if (lambda < 0.0 || lambda > 1.0/3.0) {lattice_type="simple_cubic"; lambda = 1.0/6.0; cout <<"'lambda'-value out of bounds 0..1/3. Default 1/6 is used instead" << endl; }
					else {
						if (lambda>0.33 && lambda <0.34) {lattice_type="FCC"; lambda=1.0/3.0;}
						if (lambda>0.16 && lambda < 0.167) {lattice_type="simple_cubic"; lambda = 1.0/6.0;}
						if (lambda>0.24 && lambda < 0.26) {lattice_type="hexagonal"; lambda = 1.0/4.0;}
					}
					Z=In[0]->Get_int(GetValue("Z"),6);
					if (Z<3 || Z>12) {lattice_type="simple_cubic"; lambda = 1.0/6.0; cout <<"coordination number 'Z' is out of bounds 3 .. 12. Default Z=6 is used instead"<<endl; }
				}
				if (lambda>0) {
					if (GetValue("Z").size()>0) {
						Z=In[0]->Get_int(GetValue("Z"),6);
						if (Z<3 || Z>12) {lattice_type="simple_cubic"; lambda = 1.0/6.0; cout <<"coordination number 'Z' is out of bounds 3 .. 12. Default Z=6 is used instead"<<endl; }
						else {
							lambda = 1.0/Z;
							if (Z==3) lattice_type="FCC";
							if (Z==6) lattice_type="simple_cubic";
							if (Z==4) lattice_type="hexagonal";
						}
					}
				}
				if (lambda==0) {
					cout<<"'lattice_type'-info is missing. 'lambda'-info is missing. 'Z'-info is missing. Default values; lambda=1/6, Z=6, lattice_type='simple cubic' used instead. " << endl;
					lattice_type  = "simple_cubic"; lambda=1.0/6.0; Z=0;
				}
			}
		}
		offset_first_layer =0;
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

				options.clear();
				options.push_back("mirror"); //options.push_back("mirror_2");
				//options.push_back("surface");
				if (GetValue("lowerbound_x").size()>0) {success=false; cout << "lowerbound_x is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("lowerbound_y").size()>0) {success=false; cout << "lowerbound_y is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("lowerbound_z").size()>0) {success=false; cout << "lowerbound_z is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("upperbound_x").size()>0) {success=false; cout << "upperbound_x is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("upperbound_y").size()>0) {success=false; cout << "upperbound_y is not allowed in 1-gradient calculations" << endl;}
				if (GetValue("upperbound_z").size()>0) {success=false; cout << "upperbound_z is not allowed in 1-gradient calculations" << endl;}

				Value.clear();
				Value=GetValue("lowerbound");
				if (Value.length() > 0 && !In[0]->Get_string(Value,BC[0],options,"for 'lowerbound' boundary condition not recognized."))
					success = false;
				else //default to
					BC[0]="mirror";

				Value.clear();
				Value=GetValue("upperbound");
				if (Value.length()>0 && !In[0]->Get_string(Value,BC[1],options,"for 'upperbound' boundary condition not recognized."))
					success = false;
				else //default to
					BC[1] = "mirror";
				break;
			case 2:
				if (lattice_type == "") {
					success=false;
					cout <<" in two gradient calculations, you should set lattice type to either 'simple_cubic' or 'FCC'"<<endl;
				}
				if (Z==0) {
					lattice_type  = "FCC";
					lambda=1.0/3.0; Z=3;
					cout <<"Correction: in two-gradient case the default of the lattice is 'FCC'. " << endl;
				}
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
				options.push_back("mirror"); //options.push_back("mirror_2");
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
				
				Value.clear();
				Value=GetValue("lowerbound_x");
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,BC[0],options,"for 'lowerbound_y' boundary condition not recognized."))
						success=false;
				} else
					BC[0]="mirror";
	

				Value.clear();
				Value=GetValue("upperbound_x");
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,BC[1],options,"for 'upperbound_x' boundary condition not recognized."))
						success=false;
				} else
					BC[1]="mirror";

				Value.clear();
				Value=GetValue("lowerbound_y");
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,BC[2],options,"for 'lowerbound_y' boundary condition not recognized."))
						success=false;
				} else
					BC[2]="mirror";

				Value.clear();
				Value=GetValue("upperbound_y");
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,BC[3],options,"for 'upperbound_y' boundary condition not recognized."))
						success=false;
				} else
					BC[3]="mirror";


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
				if (!In[0]->Get_int(GetValue("n_layers_x"),MX,1,1e6,"In 'lat' the parameter 'n_layers_x' is required"))
					success=false;
				if (!In[0]->Get_int(GetValue("n_layers_y"),MY,1,1e6,"In 'lat' the parameter 'n_layers_y' is required"))
					success=false;
				if (!In[0]->Get_int(GetValue("n_layers_z"),MZ,1,1e6,"In 'lat' the parameter 'n_layers_z' is required"))
					success=false;

				options.clear();
				options.push_back("mirror"); //options.push_back("mirror_2");
				options.push_back("periodic");
				//options.push_back("shifted_mirror");

				Value.clear();
				Value=GetValue("lowerbound_x");
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,BC[0],options,"for 'lowerbound_y' boundary condition not recognized."))
						success=false;
				} else
					BC[0]="mirror";
	

				Value.clear();
				Value=GetValue("upperbound_x");
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,BC[1],options,"for 'upperbound_x' boundary condition not recognized."))
						success=false;
				} else
					BC[1]="mirror";

				Value.clear();
				Value=GetValue("lowerbound_y");
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,BC[2],options,"for 'lowerbound_y' boundary condition not recognized."))
						success=false;
				} else
					BC[2]="mirror";

				Value.clear();
				Value=GetValue("upperbound_y");
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,BC[3],options,"for 'upperbound_y' boundary condition not recognized."))
						success=false;
				} else
					BC[3]="mirror";

				Value.clear();
				Value=GetValue("lowerbound_z");
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,BC[4],options,"for 'lowerbound_z' boundary condition not recognized."))
						success=false;
				} else
					BC[4]="mirror";

				Value.clear();
				Value=GetValue("upperbound_z");
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,BC[5],options,"for 'upperbound_z' boundary condition not recognized."))
						success=false;
				} else
					BC[5]="mirror";

				if (BC[0]=="periodic" || BC[1]=="periodic") {
					if (BC[0] != BC[1]) {
						cout <<"In x-direction the boundary conditions do not match:" + BC[0] << " and " <<  BC[1] << endl;
						success=false;
					}
				}
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
		FJC=3;	fjc=1;
		if (success && GetValue("FJC_choices").length()>0) {
			if (!In[0]->Get_int(GetValue("FJC_choices"),FJC,"FJC_choices can adopt only few integer values: 3 + i*4, with i = 1, 2, 3, 4, ..."))
				success=false;
			else if (FJC %4 != 3) {
				cout << "FJC_choices can adopt only few integer values: 3 + i*4, with i = 1, 2, 3, 4, ...." <<endl;
				success=false;
			}
			if (success && (gradients>1)) {
				cout << "FJC_choices can only be used in combination with 'gradients == planar' " << endl;
				success=false;
			}
			if (success && (lattice_type != "hexagonal")) {
				cout << "FJC_choices can only be used in combination with 'lattice_type == hexagonal' "<<endl;
				success=false;
			}
			fjc=(FJC-1)/2;
		}

		//Initialize system size and indexing
		PutM();
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
		num_of_steps = (end_value-VarInitValue)/step+1;
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

void Lattice::TimesL(Real* X){
	switch(gradients) {
		case 1:
			//fall-through
		case 2:
			if (geometry!="planar") Times(X,X,L,M);
			break;
		case 3:
			break;
		default:
			break;
	}
}

void Lattice::DivL(Real* X){
	switch(gradients) {
		case 1:
			//fall-through
		case 2:
			if (geometry!="planar") Div(X,L,M);
			break;
		case 3:
			break;
		default:
			break;
	}
}

Real Lattice:: Moment(Real* X,int n) {
	Real Result=0;
	Real cor; 
	if (gradients !=1 || geometry!="planar" ) {cout << "Moments analysis can only be done on one-gradient, planar system " << endl;  return Result; }
	remove_bounds(X);
	for (int i = fjc; i<M; i++) {
		cor = 1.0*(i-fjc+1)/fjc; Result += pow(cor,n)*X[i]; 
	}	
	return Result/fjc;
}

Real Lattice::WeightedSum(Real* X){
	Real sum{0};
	remove_bounds(X);
	switch(gradients) {
		case 1:
			if (fjc==1 && geometry=="planar")
				Sum(sum,X,M);
			else
				Dot(sum,X,L,M);
			break;
		case 2:
			if (geometry=="planar")
				Sum(sum,X,M);
			else
				Dot(sum,X,L,M);
			break;
		case 3:
			Sum(sum,X,M);
			break;
		default:
			return 0;
			break;
	}
	return sum;
}

void Lattice::vtk(string filename, Real* X, string id) {
if (debug) cout << "vtk in lattice " << endl;
	FILE *fp;
	int i,j,k;
	fp = fopen(filename.c_str(),"w+");
	switch(gradients) {
		case 1:
			cout << "for system with one gradient there is no VTK output available " << endl;
			break;
		case 2:
			fprintf(fp,"# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\n",MX,MY,1);
			fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",MX*MY);
			fprintf(fp,"SCALARS %s double\nLOOKUP_TABLE default\n",id.c_str());
			for (i=1; i<MX+1; i++)
			for (j=1; j<MY+1; j++)
			fprintf(fp,"%f\n",X[i*JX+j]);
			break;
		case 3:
			fprintf(fp,"# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\n",MX,MY,MZ);
			fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",MX*MY*MZ);
			fprintf(fp,"SCALARS %s double\nLOOKUP_TABLE default\n",id.c_str());
			for (i=1; i<MX+1; i++)
			for (j=1; j<MY+1; j++)
			for (k=1; k<MZ+1; k++)
			fprintf(fp,"%f\n",X[i*JX+j*JY+k]);
			break;
		default:
			break;
	}
	fclose(fp);
}
void Lattice::PutProfiles(FILE* pf,vector<Real*> X){
if (debug) cout <<"PutProfiles in lattice " << endl;
	int x,y,z,i;
	int length=X.size();
	switch(gradients) {
		case 1:
			for (x=0; x<M; x++){
				if (fjc==1) fprintf(pf,"%i\t",x); else fprintf(pf,"%e\t",1.0*(x+1)/fjc-1.0);
				for (i=0; i<length; i++) fprintf(pf,"%.20g\t",X[i][x]);
				fprintf(pf,"\n");
			}
			break;
		case 2:
			for (x=1; x<MX+1; x++)
			for (y=1; y<MY+1; y++){
				fprintf(pf,"%i\t%i\t",x,y);
				for (i=0; i<length; i++) fprintf(pf,"%.20g\t",X[i][x*JX+y]);
				fprintf(pf,"\n");
			}
			break;
		case 3:
			for (x=1; x<MX+1; x++)
			for (y=1; y<MY+1; y++)
			for (z=1; z<MZ+1; z++) {
				fprintf(pf,"%i\t%i\t%i\t",x,y,z);
				for (i=0; i<length; i++) fprintf(pf,"%.20g\t",X[i][x*JX+y*JY+z]);
				fprintf(pf,"\n");
			}
			break;
		default:
			break;
	}
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
void Lattice::PushOutput() {
if (debug) cout <<"PushOutput in lattice " << endl;
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
	push("lattice_type",lattice_type);
	push("bond_length",bond_length);
	push("FJC_choices",FJC);
	string s="profile;0"; push("L",s);

	switch (gradients) {
		case 3:
			push("n_layers_z",MZ);
			if (BZ1==1) push("lowerbound_z",mirror);
			if (BZM==MZ-1) push("upperbound_z",mirror);

			if (BX1==MX) push("lowerbound_x",periodic);
			if (BXM==1) push("upperbound_x",periodic);
			if (BY1==MY) push("lowerbound_y",periodic);
			if (BYM==1) push("upperbound_y",periodic);
			if (BZ1==MZ) push("lowerbound_z",periodic);
			if (BZM==1) push("upperbound_z",periodic);
			// Fall through
		case 2:
			push("n_layers_y",MY);
			if (BY1==1) push("lowerbound_x",mirror);
			if (BYM==MY-1) push("upperbound_x",mirror);
			// Fall through
		case 1:
			push("n_layers",MX);
			if (BX1==1) push("lowerbound",mirror);
			if (BXM==MX-1) push("upperbound",mirror);
			break;
		default:
			break;
	}
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



void Lattice::Side(Real *X_side, Real *X, int M) { //this procedure should use the lambda's according to 'lattice_type'-, 'lambda'- or 'Z'-info;
if (debug) cout <<" Side in lattice " << endl;
	Real* SIDE;
	int j,k;
	Zero(X_side,M);
	switch(gradients) {
		case 1:
			set_bounds(X);
			if (fjc==1) {
				if (geometry=="planar") {
					Add(X_side+1,X,M-1); Add(X_side,X+1,M-1);
					if (lattice_type=="simple_cubic")  {
						YplusisCtimesX(X_side,X,4.0,M);
						Norm(X_side,lambda,M);
					} else {
						Add(X_side,X,M);
						Norm(X_side,lambda,M);
					}
				} else {
					AddTimes(X_side,X,lambda0,M);
					AddTimes(X_side+1,X,lambda_1+1,M-1);
					AddTimes(X_side,X+1,lambda1,M-1);
					Norm(X_side,lambda,M);
				}
			} else {
				for (j=0; j<FJC/2; j++) {
					k=(FJC-1)/2-j;
					AddTimes(X_side+k,X,LAMBDA+j*M+k,M-k);
					AddTimes(X_side,X+k,LAMBDA+(FJC-j-1)*M,M-k);
				}
				AddTimes(X_side,X,LAMBDA+((FJC-1)/2)*M,M);
			}
			break;
		case 2:
			//set_bounds(X);
			if (geometry=="planar") {
				if (lattice_type =="simple_cubic" ) {
					Add(X_side,X+1,M-1);
					Add(X_side+1,X,M-1);
					Add(X_side+JX,X,M-JX);
					Add(X_side,X+JX,M-JX);
					YplusisCtimesX(X_side,X,2.0,M);
					Norm(X_side,lambda,M);
				} else {
					Add(X_side+1,X,M-1);
					Add(X_side,X+1,M-1);
					Add(X_side+JX,X,M-JX);
					Add(X_side,X+JX,M-JX);
					Add(X_side+1+JX,X,M-1-JX);
					Add(X_side+1,X+JX,M-1-JX);
					Add(X_side,X+1+JX,M-1-JX);
					Add(X_side+JX,X+1,M-1-JX);
					Add(X_side,X,M);
					Norm(X_side,lambda*lambda,M);
				}
			} else {
				if (lattice_type =="simple_cubic") {
					Add(X_side,X+1,M-1);
					Add(X_side+1,X,M-1);
					YplusisCtimesX(X_side,X,2.0,M);
					AddTimes(X_side+JX,X,lambda_1+JX,M-JX);
					AddTimes(X_side,X+JX,lambda1,M-JX);
					Norm(X_side,lambda,M);

				} else { cout <<"Warning: In FCC the function side is not correctly implemented" << endl;
					SIDE=new Real[M];
					Zero(SIDE,M);
					Add(SIDE,X+1,M-1);
					Add(SIDE+1,X,M-1);
					Add(SIDE,X,M);
					AddTimes(X_side,SIDE,lambda0,M);
					AddTimes(X_side+JX,SIDE,lambda_1+JX,M-JX);
					AddTimes(X_side,SIDE+JX,lambda1,M-JX);
					Norm(X_side,lambda*lambda,M);
					delete [] SIDE;
				}
			}
			break;
		case 3:
			//set_bounds(X);
			Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
			Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
			Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
			Norm(X_side,1.0/6.0,M);
			break;
		default:
			break;
	}
}

void Lattice::propagate(Real *G, Real *G1, int s_from, int s_to,int M) { //this procedure should function on simple cubic lattice.
if (debug) cout <<" propagate in lattice " << endl;
	Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	int JX_=JX, JY_=JY;
	int k=sub_box_on;
	int kk;
	int j;
	switch(gradients) {
		case 1:
			Zero(gs,M); set_bounds(gs_1);

			if (fjc==1) {
				if (geometry=="planar") {
					Add(gs+1,gs_1,M-1); Add(gs,gs_1+1,M-1);
					YplusisCtimesX(gs,gs_1,4.0,M);
					Norm(gs,lambda,M); Times(gs,gs,G1,M);
				} else {
					AddTimes(gs,gs_1,lambda0,M);
					AddTimes(gs+1,gs_1,lambda_1+1,M-1); //is this correct?
					AddTimes(gs,gs_1+1,lambda1,M-1);
					Norm(gs,lambda,M); Times(gs,gs,G1,M);
				}
			} else {
				for (j=0; j<FJC/2; j++) {
					kk=(FJC-1)/2-j;
					AddTimes(gs+kk,gs_1,LAMBDA+j*M+kk,M-kk);
					AddTimes(gs,gs_1+kk,LAMBDA+(FJC-j-1)*M,M-kk);
				}
				AddTimes(gs,gs_1,LAMBDA+(FJC-1)/2*M,M);
				Times(gs,gs,G1,M);
			}
			break;
		case 2:
			Zero(gs,M); set_bounds(gs_1);
			if (geometry=="planar") {
				Add(gs,gs_1+1,M-1);
				Add(gs+1,gs_1,M-1);
				Add(gs+JX,gs_1,M-JX);
				Add(gs,gs_1+JX,M-JX);
				YplusisCtimesX(gs,gs_1,2.0,M);
				Norm(gs,lambda,M); Times(gs,gs,G1,M);
			} else {
				Add(gs,gs_1+1,M-1);
				Add(gs+1,gs_1,M-1);
				YplusisCtimesX(gs,gs_1,2.0,M);
				AddTimes(gs+JX,gs_1,lambda_1+JX,M-JX);
				AddTimes(gs,gs_1+JX,lambda1,M-JX);
				Norm(gs,lambda,M);Times(gs,gs,G1,M);
			}
			break;
		case 3:
			if (k>0) {JX_=jx[k]; JY_=jy[k];}
			Zero(gs,M); set_bounds(gs_1);
			Add(gs+JX_,gs_1,M-JX_); Add(gs,gs_1+JX_,M-JX_);
			Add(gs+JY_,gs_1,M-JY_); Add(gs,gs_1+JY_,M-JY_);
			Add(gs+1,gs_1,M-1);  Add(gs,gs_1+1, M-1);
			Norm(gs,1.0/6.0,M); Times(gs,gs,G1,M);
			break;
		default:
			break;
	}
}

//Specify which types can be used by remove_bounds
template void Lattice::remove_bounds<int>(int*);
template void Lattice::remove_bounds<Real>(Real*);

template <typename T>
void Lattice::remove_bounds(T *X){
if (debug) cout <<" remove_bounds (Real) in lattice " << endl;
	int x,y;
	int j;
	switch(gradients) {
		case 1:
			if (fjc==1) {
				X[0]=0; X[MX+1]=0;
			} else {
				for (j=0; j<fjc; j++) {
					X[j]=0;
					X[(MX+2)*fjc-j-1]=0;
				}
			}
			break;
		case 2:
			for (x=1; x<MX+1; x++) {
				X[x*JX+0] = 0;
				X[x*JX+MY+1]=0;
			}
			for (y=1; y<MY+1; y++) {
				X[0+y] = 0;
				X[(MX+1)*JX+y]=0;
			}
			X[0]             =0;
			X[(MX+1)*JX]     =0;
			X[0+MY+1]        =0;
			X[(MX+1)*JX+MY+1]=0;
			break;
		case 3:
			if (sub_box_on!=0) {int k=sub_box_on;
				for (int i=0; i<n_box[k]; i++)
				RemoveBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
			} else {
			RemoveBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
			}
			break;
		default:

			break;
	}
}

//Specify which types can be used by set_bounds
template void Lattice::set_bounds<int>(int*);
template void Lattice::set_bounds<Real>(Real*);

template <typename T>
void Lattice::set_bounds(T *X){
if (debug) cout <<"set_bounds (Reals) in lattice " << endl;
	int x,y;
	int j;
	switch(gradients) {
		case 1:
			if (fjc==1) {
				X[0]=X[BX1]; X[MX+1]=X[BXM];
			} else {
				for (j=0; j<fjc; j++) {
					X[j]=X[(BX1+1)*fjc-j-1];
					X[(MX+2)*fjc-j-1]=X[BXM*fjc+j];
				}
			}
			break;
		case 2:
			for (x=1; x<MX+1; x++) {
				X[x*JX+0] = X[x*JX+BY1];
				X[x*JX+MY+1]=X[x*JX+BYM];
			}
			for (y=1; y<MY+1; y++) {
				X[0+y] = X[BX1*JX+y];
				X[(MX+1)*JX+y]=X[BXM*JX+y];
			}
			X[0]             =X[1]*X[JX]; if (X[0]>0) X[0]=sqrt(X[0]);
			X[(MX+1)*JX]     =X[MX*JX]*X[(MX+1)*JX+1]; if (X[(MX+1)*JX]>0) X[(MX+1)*JX]=sqrt(X[(MX+1)*JX]);
			X[0+MY+1]        =X[JX+MY+1]*X[0+MY]; if (X[MY+1]>0) X[MY+1]=sqrt(X[MY+1]);
			X[(MX+1)*JX+MY+1]=X[MX*JX+MY+1]*X[(MX+1)*JX+MY]; if (X[(MX+1)*JX+MY+1]>0) X[(MX+1)*JX+MY+1]=sqrt(X[(MX+1)*JX+MY+1]);
			break;
		case 3:
			if (sub_box_on!=0) {
				int k=sub_box_on;
				for (int i=0; i<n_box[k]; i++)
					SetBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
			} else
				SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
			break;
		default:
			break;
	}
}


bool Lattice::ReadRange(int* r, int* H_p, int &n_pos, bool &block, string range, string seg_name, string range_type) {
if (debug) cout <<"ReadRange in lattice " << endl;
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
		switch(gradients) {
			case 1:
				if (coor.size()!=1) {cout << "In mon " + seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes do not come as a single coordinate 'x'" << endl; success=false;}
				else {
					diggit=coor[0].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[0]=In[0]->Get_int(coor[0],0); else {
						recognize_keyword=false;
						if (coor[0]=="firstlayer") {recognize_keyword=true; r[0] = 1;}
						//if (coor[0]=="lowerbound") {recognize_keyword=true; r[0] = 0;}
						//if (coor[0]=="upperbound") {recognize_keyword=true; r[0] = MX+1;}
						if (coor[0]=="lastlayer")  {recognize_keyword=true; r[0] = MX;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
							success=false;
						}
					}
					if (r[0] < 1-a || r[0] > MX+a) {cout << "In mon " + seg_name + ", for 'pos 1', the x-coordinate in '" + range_type + "' is out of bounds: "<< 1-a <<" .." << MX+a << endl; success =false;}
				}
				coor.clear(); In[0]->split(set[1],',',coor);

				if (coor.size()!=1) {cout << "In mon " + seg_name+ ", for 'pos 2', in '" + range_type + "' the coordinates do not come as a single coordinate 'x'" << endl; success=false;}
				else {
					diggit=coor[0].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[3]=In[0]->Get_int(coor[0],0); else {
						recognize_keyword=false;
						if (coor[0]=="firstlayer") {recognize_keyword=true; r[3] = 1;}
						//if (coor[0]=="lowerbound") {recognize_keyword=true; r[3] = 0;}
						//if (coor[0]=="upperbound") {recognize_keyword=true; r[3] = MX+1;}
						if (coor[0]=="lastlayer")  {recognize_keyword=true; r[3] = MX;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer'" << endl;
							success=false;
						}
					}
					if (r[3] < 1-a || r[3] > MX+a) {cout << "In mon " + seg_name+ ", for 'pos 2', the x-coordinate in '" + range_type + "' is out of bounds; "<< 1-a <<" .." << MX+a << endl; success =false;}
					if (r[0] > r[3]) {cout << "In mon " + seg_name+ ", for 'pos 1', the x-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
				}
				break;
			case 2:
				if (coor.size()!=2) {cout << "In mon " + 	seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes do not come in set of two: 'x,y'" << endl; success=false;}
				else {
					diggit=coor[0].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[0]=In[0]->Get_int(coor[0],0); else {
						recognize_keyword=false;
						if (coor[0]=="firstlayer") {recognize_keyword=true; r[0] = 1;}
						//if (coor[0]=="lowerbound") {recognize_keyword=true; r[0] = 0;}
						//if (coor[0]=="upperbound") {recognize_keyword=true; r[0] = MX+1;}
						if (coor[0]=="lastlayer")  {recognize_keyword=true; r[0] = MX;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
							success=false;
						}
					}
					if (r[0] < 1-a || r[0] > MX+a) {cout << "In mon " + seg_name + ", for 'pos 1', the x-coordinate in '" + range_type + "' is out of bounds: "<< 1-a <<" .." << MX+a << endl; success =false;}
					diggit=coor[1].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[1]=In[0]->Get_int(coor[1],0); else {
						recognize_keyword=false;
						if (coor[1]=="firstlayer") {recognize_keyword=true; r[1] = 1;}
						//if (coor[1]=="lowerbound") {recognize_keyword=true; r[1] = 0;}
						//if (coor[1]=="upperbound") {recognize_keyword=true; r[1] = MY+1;}
						if (coor[1]=="lastlayer")  {recognize_keyword=true; r[1] = MY;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
							success=false;
						}
					}
					if (r[1] < 1-a || r[1] > MY+a) {cout << "In mon " + seg_name+ ", for 'pos 1', the y-coordinate in '" + range_type + "' is out of bounds: "<< 1-a <<" .."<< MY+a << endl; success =false;}
				}
				coor.clear(); In[0]->split(set[1],',',coor);

				if (coor.size()!=2) {cout << "In mon " + seg_name+ ", for 'pos 2', in '" + range_type + "', the coordinates do not come in set of two: 'x,y'" << endl; success=false;}
				else {
					diggit=coor[0].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[3]=In[0]->Get_int(coor[0],0); else {
						recognize_keyword=false;
						if (coor[0]=="firstlayer") {recognize_keyword=true; r[3] = 1;}
						//if (coor[0]=="lowerbound") {recognize_keyword=true; r[3] = 0;}
						//if (coor[0]=="upperbound") {recognize_keyword=true; r[3] = MX+1;}
						if (coor[0]=="lastlayer")  {recognize_keyword=true; r[3] = MX;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
							success=false;
						}
					}
					if (r[3] <1-a || r[3] > MX+a) {cout << "In mon " + seg_name+ ", for 'pos 2', the x-coordinate in '" + range_type + "' is out of bounds: " << 1-a << " .. " << MX+a<< endl; success=false;}
					diggit=coor[1].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[4]=In[0]->Get_int(coor[1],0); else {
						recognize_keyword=false;
						if (coor[1]=="firstlayer") {recognize_keyword=true; r[4] = 1;}
						//if (coor[1]=="lowerbound") {recognize_keyword=true; r[4] = 0;}
						//if (coor[1]=="upperbound") {recognize_keyword=true; r[4] = MY+1;}
						if (coor[1]=="lastlayer")  {recognize_keyword=true; r[4] = MY;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' " << endl;
							success=false;
						}
					}
					if (r[4] < 1-a || r[4] > MY+a) {cout << "In mon " + seg_name+ ", for 'pos 2', the y-coordinate in '" + range_type + "' is out of bounds; "<< 1-a <<" .."<< MY+a << endl; success =false;}
					if (r[0] > r[3]) {cout << "In mon " + seg_name+ ", for 'pos 1', the x-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
					if (r[1] > r[4]) {cout << "In mon " + seg_name+ ", for 'pos 1', the y-coordinate in '" + range_type + "' should be less than that of 'pos 2'" << endl; success =false;}
				}
				break;
			case 3:
				if (coor.size()!=3) {cout << "In mon " + 	seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes do not come in set of three: 'x,y,z'" << endl; success=false;}
				else {
					diggit=coor[0].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[0]=In[0]->Get_int(coor[0],0); else {
						recognize_keyword=false;
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
				break;
			default:
				break;
		}
	} else {
		string s;
		In[0]->split(set[0],')',coor);
		s=coor[0].substr(0,1);
		if (s!="(") { //now expect one of the keywords
				block=true;
			if (gradients==1) {
				recognize_keyword=false;
				if (coor[0]=="firstlayer")         {recognize_keyword=true; r[0] = r[3] = 1; }
				//if (coor[0]=="lowerbound" && a==1) {recognize_keyword=true; r[0] = r[3] = 0; }
				//if (coor[0]=="upperbound" && a==1) {recognize_keyword=true; r[0] = r[3] = MX+1; }
				if (coor[0]=="lastlayer")          {recognize_keyword=true; r[0] = r[3] = MX;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the input: 'firstlayer' or 'lastlayer' ";
					if (a==1) cout << "or 'upperbound' or 'lowerbound'";
					cout <<" is expected. " << endl;
					success=false;
				}
			} else {
				cout << "In mon " + seg_name + " and  range_type " + range_type + ", the info was not recognised because when 'gradients>1' the lonely keywords 'firstlayer' 'lastlayers' do not work." << endl;
				success=false;
			}
		} else {
			int px{0},py{0},pz{0};
			string s;

			if (coor.size()==0)
				block=false;
			else {
				for (size_t i = 0 ; i < coor.size() ; ++i) {
					s=coor[i].substr(1,coor[i].size()-1);
					In[0]->split(s,',',xyz);
					switch(gradients) {
					case 3:
								if (xyz.size()!=3) {
									cout << "In mon " + seg_name+ " pinned_range  the expected 'triple coordinate' -with brackets- structure '(x,y,z)' was not found. " << endl;  success = false;
									break;
								}
									px=In[0]->Get_int(xyz[0],0);
									if (px < 1 || px > MX) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the x-coordinate in pinned_range out of bounds: 1.." << MX << endl; success =false;}
									py=In[0]->Get_int(xyz[1],0);
									if (py < 1 || py > MY) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 1.." << MY << endl; success =false;}
									pz=In[0]->Get_int(xyz[2],0);
									if (pz < 1 || pz > MZ) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 1.." << MZ << endl; success =false;}
								H_p[i]=px*JX+py*JY+pz;
						break;
					case 2:
								if (xyz.size()!=2) {
									cout << "In mon " + seg_name+ " pinned_range  the expected 'pair of coordinate' -with brackets- structure '(x,y)' was not found. " << endl;  success = false;
								} else {
									px=In[0]->Get_int(xyz[0],0);
									if (px < 0 || px > MX+1) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the x-coordinate in pinned_range out of bounds: 0.." << MX+1 << endl; success =false;}
									py=In[0]->Get_int(xyz[1],0);
									if (py < 0 || py > MY+1) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 0.." << MY+1 << endl; success =false;}
								}
								H_p[i]=px*JX+py;
						break;
					case 1:
							if (xyz.size()!=1) {
								cout << "In mon " + seg_name+ " pinned_range  the expected 'single coordinate' -with brackets- structure '(x)' was not found. " << endl;  success = false;
							} else {
								px=In[0]->Get_int(xyz[0],0);
								if (px < 0 || px > MX+1) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the x-coordinate in pinned_range out of bounds: 0.." << MX+1 << endl; success =false;}
							}
							H_p[i]=px;
					break;
				default:
					break;
				}
			}
		}

		}
	}
	return success;
}

bool Lattice::ReadRangeFile(string filename,int* H_p, int &n_pos, string seg_name, string range_type) {
if (debug) cout <<"ReadRangeFile in lattice " << endl;
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

	switch(gradients) {
		case 1:
				if (length == MX) { //expect to read 'mask file';
					if (n_pos==0) {
						for (i = 0 ; i < length ; ++i) {
							if (In[0]->Get_int(lines[i],0)==1) n_pos++;
						}
						if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any unities." << endl;}
					} else {
						p_i=0;
						for (x=1; x<MX+1; x++) {
							if (In[0]->Get_int(lines[x-1],0)==1) {H_p[p_i]=x; p_i++;}
						}
					}
				} else { //expect to read x only
					px=0,py=0; i=0;
					if (n_pos==0) n_pos=length;
					else {
						while (i<length) {
							xyz.clear();
							In[0]->split(lines[i],',',xyz);
							length_xyz=xyz.size();
							if (length_xyz!=1) {
								cout << "In mon " + seg_name + " " +range_type+"_filename  the expected 'single coordinate' 'x' was not found. " << endl;  success = false;
							} else {
								px=In[0]->Get_int(xyz[0],0);
								if (px < 1 || px > MX) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the x-coordinate in "+range_type+"_filename out of bounds: 1.." << MX << endl; success =false;}
							}
							H_p[i]=px;
							i++;
						}
					}
				}
			break;
		case 2:
				if (length == MX*MY) { //expect to read 'mask file';
					i=0;
					if (n_pos==0) {
						while (i<length){
							if (In[0]->Get_int(lines[i],0)==1) n_pos++;
							i++;
						};
						if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any unities." << endl;}
					} else {
						i=0; p_i=0;
						for (x=1; x<MX+1; x++) for (y=1; y<MY+1; y++)  {
							if (In[0]->Get_int(lines[i],0)==1) {H_p[p_i]=x*JX+y; p_i++;}
							i++;
						}
					}
				} else { //expect to read x,y
					px=0,py=0; i=0;
					if (n_pos==0) n_pos=length;
					else {
						while (i<length) {
							xyz.clear();
							In[0]->split(lines[i],',',xyz);
							length_xyz=xyz.size();
							if (length_xyz!=2) {
								cout << "In mon " + seg_name + " " +range_type+"_filename  the expected 'pair of coordinates' 'x,y' was not found. " << endl;  success = false;
							} else {
								px=In[0]->Get_int(xyz[0],0);
								if (px < 1 || px > MX) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the x-coordinate in "+range_type+"_filename out of bounds: 1.." << MX << endl; success =false;}
								py=In[0]->Get_int(xyz[1],0);
								if (py < 1 || py > MY) {cout << "In mon " + seg_name + ", for 'pos' "<< i << ", the y-coordinate in "+range_type+"_filename out of bounds: 1.." << MY << endl; success =false;}
							}
							H_p[i]=px*JX+py;
							i++;
						}
					}
				}
			break;
		case 3:
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
							if (In[0]->Get_int(lines[i],0)==1) {H_p[p_i]=x*JX+y*JY+z; p_i++;}
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
							H_p[i]=px*JX+py*JY+pz;
							i++;
						}
					}
				}
			break;
		default:
			break;
	}
	return success;
}

bool Lattice::CreateMASK(int* H_MASK, int* r, int* H_P, int n_pos, bool block) {
if (debug) cout <<"CreateMask for lattice " + name << endl;
	bool success=true;
	H_Zero(H_MASK,M);
	if (block) {
		switch(gradients) {
			case 1:
				if (fjc==1) {
					for (int x=r[0]; x<r[3]+1; x++)
						H_MASK[x]=1;
				} else {
					for (int x=r[0]*fjc; x<(r[3]+1)*fjc; x++)
						H_MASK[x]=1;
				}
			break;
			case 2:
				for (int x=r[0]; x<r[3]+1; x++)
					for (int y=r[1]; y<r[4]+1; y++)
						H_MASK[x*JX+y]=1;
			break;
			case 3:
				for (int x=r[0]; x<r[3]+1; x++)
					for (int y=r[1]; y<r[4]+1; y++)
						for (int z=r[2]; z<r[5]+1; z++)
							H_MASK[x*JX+y*JY+z]=1;
			break;
			default:
			break;
		}
	} else {
		for (int i = 0; i<n_pos; i++) H_MASK[H_P[i]]=1;
	}
	return success;
}

bool Lattice::GenerateGuess(Real* x, string CalculationType, string GuessType, Real A_value, Real B_value) {
if (debug) cout <<"GenerateGuess in lattice " << endl;
//GuessType: lamellae,Im3m,FCC,BCC,HEX,gyroid,Real_gyroid,Real_diamond,perforated_lamellae
//CalculationType: micro_emulsion,micro_phasesegregation
	bool success = true;
	int i,j,k;
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

Real Lattice::ComputeTheta(Real* phi) {
	Real result=0; remove_bounds(phi);
	if (gradients<3 && geometry !="planar") Dot(result,phi,L,M); else {if (fjc==1) Sum(result,phi,M); else  Dot(result,phi,L,M);}
	return result;
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
			for (int i=0; i<num_1; i++) {in_file >> name;
				 monlist.push_back(name);
			}
			in_file>>num_2;
			for (int i=0; i<num_2; i++) {in_file>> name;
				statelist.push_back(name);
			}
//if (readx==1) cout <<"again..."<< endl;
//cout <<"method " << method << endl;
//cout <<"mx , my , mz " << mx <<" , " << my << " , " << mz << endl;
//if (charged) cout <<"charged" << endl; else cout <<"not charged" << endl;
//for (int i=0; i<num; i++) cout << monlist[i] << endl;
			if (readx==0) readx=-2; else {
				int m;
				if (my==0) {m=(mx+2);} else {if (mz==0) m=(mx+2)*(my+2); else {m=(mx+2)*(my+2)*(mz+2);}}
				int iv=(num_1+num_2)*m;
				if (charged) iv +=m;
				for (int i=0; i<iv; i++) {
					in_file >> x[i];  //cout <<"R i " << i << " " << x[i] << endl;
				}
				readx=-3;
			}
		}
		in_file.close();
//if (readx==-3) for (int i=0; i<M; i++) cout <<"lat i " << i << " " << x[3*M+i] << endl;
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



void Lattice::UpdateEE(Real* EE, Real* psi) {
	int x,y;
	Real Exmin,Explus,Eymin,Eyplus;
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_bT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.

	switch (gradients) {
		case 1:
			Explus = psi[0]-psi[1];
			Explus *= Explus;
			for (x=1; x<M-1; x++) {
				Exmin = Explus;
				Explus = psi[x]-psi[x+1];
				Explus *= Explus;
				if (geometry=="planar")
					EE[x] = pf*(Exmin + Explus);
				else
					EE[x] = pf*(lambda_1[x]*Exmin+lambda1[x]*Explus);
			}
			if (fjc>1) cout <<"Error in updateEE " << endl;
			break;
		case 2:
			for (x=1; x<MX+1; x++) {
				Eyplus = psi[x*JX]-psi[x*JX+1];
				Eyplus *=Eyplus;
				for (y=1; y<MY+1; y++) {
					Eymin=Eyplus;
					Eyplus = psi[x*JX+y]-psi[x*JX+y+1];
					Eyplus *=Eyplus;
					Exmin = psi[x*JX+y]-psi[(x-1)*JX+y];
					Exmin *=Exmin;
					Explus=psi[x*JX+y]-psi[(x+1)*JX+y];
					Explus *=Explus;
					if (geometry=="flat") {
						EE[x*JX+y]=pf*(Exmin+Explus+Eymin+Eyplus);
					} else {
						EE[x*JX+y]=pf*(lambda_1[x*JX+y]*Exmin+lambda1[x*JX+y]*Explus+Eymin+Eyplus);
					}
				}
			}
			break;
		case 3:
			Zero(EE,M);
			AddGradSquare(EE+1,psi,psi+1,psi+2,M-2);
			AddGradSquare(EE+JX,psi,psi+JX,psi+2*JX,M-2*JX);
			AddGradSquare(EE+JY,psi,psi+JY,psi+2*JY,M-2*JY);
			Norm(EE,pf,M);
			break;
		default:
			cout <<"GetEE wrongly implemented in lattice: " << endl;
			break;
	}

}

void Lattice::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, int* Mask) { //not only update psi but also g (from newton).
	int x,y;
	#ifndef CUDA
	int z;
	#endif

	Real epsXplus,epsXmin,epsYplus,epsYmin;
	#ifndef CUDA
	Real epsZplus, epsZmin;
	#endif

	Real C = e*e/(eps0*k_BT*bond_length);
	switch (gradients) {
		case 1:
			C*=2;
			epsXplus=(eps[0]+eps[1]);
			for (x=1; x<M-1; x++) {
				epsXmin=epsXplus;
				epsXplus=eps[x]+eps[x+1];
				if (Mask[x]==0) { //use updated psi at subsequent coordinates = Jacobian trick. Upwind to larger x.
					if (geometry == "planar" ) {
						psi[x] = (epsXmin*psi[x-1]+epsXplus*psi[x+1]+C*q[x])/(epsXplus+epsXmin);
					} else {
						psi[x] = (lambda_1[x]*epsXmin*psi[x-1]+lambda1[x]*epsXplus*psi[x+1]+C*q[x])/(lambda_1[x]*epsXplus+lambda1[x]*epsXmin);
					}
				}
			}
			epsXmin=eps[M-1]+eps[M-2];
			for (x=M-2; x>0; x--) { //perform backward update as well in case potentials need to 'diffuse' to lower x: downwind.
				epsXplus=epsXmin;
				epsXmin=eps[x]+eps[x-1];
				if (Mask[x]==0) {
					if (geometry=="planar") {
						psi[x] = (epsXmin*psi[x-1]+epsXplus*psi[x+1]+C*q[x])/(epsXplus+epsXmin);
					} else {
						psi[x] = (lambda_1[x]*epsXmin*psi[x-1]+lambda1[x]*epsXplus*psi[x+1]+C*q[x])/(lambda_1[x]*epsXplus+lambda1[x]*epsXmin);
					}
					g[x]-=psi[x];
				}
			}
			if (fjc>1) cout <<"error in UpdatePsi" << endl;
			break;
		case 2:
			C*=4;
			for (x=1; x<MX+1; x++) {
				epsYplus = eps[x*JX]+eps[x*JX+1];
				for (y=1; y<MY+1; y++) {
					epsYmin=epsYplus;
					epsYplus= eps[x*JX+y]+eps[x*JX+y+1];
					epsXmin = eps[x*JX+y]+eps[(x-1)*JX+y];
					epsXplus= eps[x*JX+y]+eps[(x+1)*JX+y];
					if (Mask[x*JX+y]==0) {
						if (geometry=="flat") {
							psi[x*JX+y]= (epsXmin*psi[(x-1)*JX+y]+epsXplus*psi[(x+1)*JX+y]+
								      epsYmin*psi[x*JX+y-1]+epsYplus*psi[x*JX+y+1]+C*q[x*JX+y])/
								      (epsXmin+epsXplus+epsYmin+epsYplus);
						} else {
							psi[x*JX+y]= (lambda_1[x*JX+y]*epsXmin*psi[(x-1)*JX+y]+lambda1[x*JX+y]*epsXplus*psi[(x+1)*JX+y]+
								      epsYmin*psi[x*JX+y-1]+epsYplus*psi[x*JX+y+1]+C*q[x*JX+y])/
								     (lambda_1[x*JX+y]*epsXmin+lambda1[x*JX+y]*epsXplus+epsYmin+epsYplus);
						}
					}
				}
			}
			for (x=MX; x>0; x--) {
				epsYmin = eps[x*JX+MY]+eps[x*JX+MY+1];
				for (y=MY; y>0; y--) {
					epsYplus=epsYmin;
					epsYmin= eps[x*JX+y]+eps[x*JX+y+1];
					epsXmin = eps[x*JX+y]+eps[(x-1)*JX+y];
					epsXplus= eps[x*JX+y]+eps[(x+1)*JX+y];
					if (Mask[x*JX+y]==0) {
						if (geometry=="flat") {
							psi[x*JX+y]= (epsXmin*psi[(x-1)*JX+y]+epsXplus*psi[(x+1)*JX+y]+
								      epsYmin*psi[x*JX+y-1]+epsYplus*psi[x*JX+y+1]+C*q[x*JX+y])/
								     (epsXmin+epsXplus+epsYmin+epsYplus);
						} else {
							psi[x*JX+y]= (lambda_1[x*JX+y]*epsXmin*psi[(x-1)*JX+y]+lambda1[x*JX+y]*epsXplus*psi[(x+1)*JX+y]+
								      epsYmin*psi[x*JX+y-1]+epsYplus*psi[x*JX+y+1]+C*q[x*JX+y])/
							             (lambda_1[x*JX+y]*epsXmin+lambda1[x*JX+y]*epsXplus+epsYmin+epsYplus);
						}
						g[x*JX+y] -=psi[x*JX+y];
					}
				}
			}
			break;
		case 3: //Do not know how to do Jacobi trick in CUDA.
			C*=6;
#ifdef CUDA
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

#endif
			break;
		default:
			cout <<"UpdatePsi wrongly implemented in lattice: " << endl;
			break;
	}
}

void Lattice::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, int* Mask) {//Not only update q (charge), but also g (from newton).
	int x,y;
	#ifndef CUDA
	int z;
	#endif

	Real epsXplus,epsXmin,epsYplus,epsYmin;
	#ifndef CUDA
	Real epsZplus, epsZmin;
	#endif
	Real C = -e*e/(eps0*k_BT*bond_length);
	switch (gradients) {
		case 1:
			C *=2.0;
			epsXplus=eps[0]+eps[1];
			for (x=1; x<M-1; x++) {
				epsXmin=epsXplus;
				epsXplus=eps[x]+eps[x+1];
				if (Mask[x]==1) {
					if (geometry=="planar") {
						q[x] = (epsXmin*psi[x-1]+epsXplus*psi[x+1]-(epsXplus+epsXmin)*psi[x])/C;
					} else {
						q[x] = (lambda_1[x]*epsXmin*psi[x-1]+lambda1[x]*epsXplus*psi[x+1]-(lambda1[x]*epsXplus+lambda_1[x]*epsXmin)*psi[x])/C;

					}
					g[x]=-q[x];
				}
			}
			if (fjc>1) cout <<"error in updateQ " << endl;
			break;
		case 2:
			C *=4.0;
			for (x=1; x<MX+1; x++) {
				epsYplus = eps[x*JX]+eps[x*JX+1];
				for (y=1; y<MY+1; y++) {
					epsYmin=epsYplus;
					epsYplus= eps[x*JX+y]+eps[x*JX+y+1];
					epsXmin = eps[x*JX+y]+eps[(x-1)*JX+y];
					epsXplus= eps[x*JX+y]+eps[(x+1)*JX+y];
					if (Mask[x*JX+y]==1) {
						if (geometry=="flat") {
							q[x*JX+y]= (epsXmin*psi[(x-1)*JX+y]+epsXplus*psi[(x+1)*JX+y]+
								epsYmin*psi[x*JX+y-1]+epsYplus*psi[x*JX+y+1]-
								(epsXmin+epsXplus+epsYmin+epsYplus)*psi[x*JX+y])/C;
						} else {
							q[x*JX+y]= (lambda_1[x*JX+y]*epsXmin*psi[(x-1)*JX+y]+lambda1[x*JX+y]*epsXplus*psi[(x+1)*JX+y]+
							 	epsYmin*psi[x*JX+y-1]+epsYplus*psi[x*JX+y+1]-
								(lambda_1[x*JX+y]*epsXmin+lambda1[x*JX+y]*epsXplus+epsYmin+epsYplus)*psi[x*JX+y])/C;
						}
						g[x*JX+y]=-q[x*JX+y];
					}
				}
			}
			break;
		case 3:
			C *=6.0;
#ifdef CUDA
			UpQ(g+JX+JY+1,q+JX+JY+1,psi+JX+JY+1,eps+JX+JY+1,JX,JY,C,Mask+JX+JY+1,M-2*(JX+JY+1));
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
#endif
			break;
		default:
			cout <<"UpdateQ wrongly implemented in lattice: " << endl;
			break;
	}
}

/*

All below is commented out. This is stored here to recover from earlier program.

void Sideh(Real *X_side, Real *X, int M) {
	Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);

	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Add(X_side+JX,X+JY,MM-JX-JY); Add(X_side+JY,X+JX,MM-JY-JX);
	Add(X_side+1,X+JY,MM-JY-1); Add(X_side+JY,X+1,MM-JY-1);
	Add(X_side+1,X+JX,MM-JX-1); Add(X_side+JX,X+1,MM-JX-1);
	Norm(X_side,1.0/12.0,M);
}

void Side(Real *X_side, Real *X, int M) {
	Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Norm(X_side,1.0/6.0,M);
}

void advanced_average(Real *X_side, Real *X, int M){
        Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);

	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);

	Add(X_side+JX+JY,X,M-JX-JY); Add(X_side, X+JX+JY, M-JX-JY);
	Add(X_side+JY+1,X,M-JY-1); Add(X_side, X+JY+1, M-JY-1);
	Add(X_side+JX+1,X,M-JX-1); Add(X_side, X+JX+1, M-JX-1);
        Add(X_side+JX-JY,X,M-JX+JY); Add(X_side, X+JX-JY, M-JX+JY);
	Add(X_side+JY-1,X,M-JY+1); Add(X_side, X+JY-1, M-JY+1);
	Add(X_side+JX-1,X,M-JX+1); Add(X_side, X+JX-1, M-JX+1);

	Add(X_side+JX+JY+1,X,M-JX-JY-1); Add(X_side, X+JX+JY+1, M-JX-JY-1);
	Add(X_side+JX+JY-1,X,M-JX-JY+1); Add(X_side, X+JX+JY-1, M-JX-JY+1);
	Add(X_side+JX-JY+1,X,M-JX+JY-1); Add(X_side, X+JX-JY+1, M-JX+JY-1);
	Add(X_side-JX+JY+1,X,M+JX-JY-1); Add(X_side, X-JX+JY+1, M+JX-JY-1);

	Norm(X_side,1.0/26.0,M);
}

void Propagateh(Real* G, Real* G1, int s_from, int s_to) { //on small boxes
	int MMM=M*n_box;
	Real *gs = G+MMM*s_to, *gs_1 = G+MMM*s_from, *g = G1;
	Zero(gs,MMM);
	for (int p=0; p<n_box; p++) SetBoundaries(gs_1+M*p,jx,jy,bx1,bxm,by1,bym,bz1,bzm,Mx,My,Mz);

	Add(gs+jx,gs_1,MMM-jx); Add(gs,gs_1+jx,MMM-jx);
	Add(gs+jy,gs_1,MMM-jy); Add(gs,gs_1+jy,MMM-jy);
	Add(gs+1,gs_1,MMM-1);  Add(gs,gs_1+1, MMM-1);
	Add(gs+jx,gs_1+jy,MMM-jx-jy); Add(gs+jy,gs_1+jx,MMM-jx-jy);
	Add(gs+1,gs_1+jy,MMM-jy-1); Add(gs+jy,gs_1+1,MMM-jy-1);
	Add(gs+1,gs_1+jx,MMM-jx-1); Add(gs+jx,gs_1+1,MMM-jx-1);
	Norm(gs,1.0/12.0,MMM); Times(gs,gs,g,MMM);
}

void PROPAGATEh(Real *G, Real *G1, int s_from, int s_to) { //on big box
	Real *gs = G+MM*(s_to), *gs_1 = G+MM*(s_from), *g = G1;
	Zero(gs,MM);
	SetBoundaries(gs_1,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(gs+JX,gs_1,MM-JX); Add(gs,gs_1+JX,MM-JX);
	Add(gs+JY,gs_1,MM-JY); Add(gs,gs_1+JY,MM-JY);
	Add(gs+1,gs_1,MM-1);  Add(gs,gs_1+1, MM-1);
	Add(gs+JX,gs_1+JY,MM-JX-JY); Add(gs+JY,gs_1+JX,MM-JX-JY);
	Add(gs+1,gs_1+JY,MM-JY-1); Add(gs+JY,gs_1+1,MM-JY-1);
	Add(gs+1,gs_1+JX,MM-JX-1); Add(gs+JX,gs_1+1,MM-JX-1);
	Norm(gs,1.0/12.0,MM); Times(gs,gs,g,MM);
}

void Propagate(Real* G, Real* G1, int s_from, int s_to) { //on small boxes
	int MMM=M*n_box;
	Real *gs = G+MMM*s_to, *gs_1 = G+MMM*s_from, *g = G1;
	Zero(gs,MMM);
	for (int p=0; p<n_box; p++) SetBoundaries(gs_1+M*p,jx,jy,bx1,bxm,by1,bym,bz1,bzm,Mx,My,Mz);
	Add(gs+jx,gs_1,MMM-jx); Add(gs,gs_1+jx,MMM-jx);
	Add(gs+jy,gs_1,MMM-jy); Add(gs,gs_1+jy,MMM-jy);
	Add(gs+1,gs_1,MMM-1);  Add(gs,gs_1+1, MMM-1);
	Norm(gs,1.0/6.0,MMM); Times(gs,gs,g,MMM);
}

void PROPAGATE(Real *G, Real *G1, int s_from, int s_to) { //on big box
	Real *gs = G+MM*(s_to), *gs_1 = G+MM*(s_from), *g = G1;
	Zero(gs,MM);
	SetBoundaries(gs_1,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(gs+JX,gs_1,MM-JX); Add(gs,gs_1+JX,MM-JX);
	Add(gs+JY,gs_1,MM-JY); Add(gs,gs_1+JY,MM-JY);
	Add(gs+1,gs_1,MM-1);  Add(gs,gs_1+1, MM-1);
	Norm(gs,1.0/6.0,MM); Times(gs,gs,g,MM);
}


void Lattice::Edis(Real *X_side, Real *X, int M) { //operation opposite to Side ;-) ....careful it distroys info in X_side.
if (debug) cout <<"Edis in lattice " << endl;
	int j;
	Zero(X_side,M); set_bounds(X);
	if (fjc==1) {
	} else {
		for (j=0; j<fjc/2; j++) {
			AddTimes(X_side+j,X,VL+j*M+j,M-j);
			AddTimes(X_side,X+j,VL+(fjc-j-1)*M,M-j);
		}
		AddTimes(X_side,X,VL+((fjc-1)/2)*M,M);
	}
	Cp(X,X_side,M);
}
*/
