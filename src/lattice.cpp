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
  	KEYS.push_back("lattice_type");
	KEYS.push_back("stencil_full");
	//KEYS.push_back("lambda");
	//KEYS.push_back("Z");
	KEYS.push_back("FJC_choices");
	sub_box_on = 0;
	all_lattice = false;
	ignore_sites=false;
	stencil_full=false;
	gradients=1;
	fjc=1;
	MX=MY=MZ=0;
	offset_first_layer=0;
}

Lattice::~Lattice() {
if (debug) cout <<"lattice destructor " << endl;
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
			free(B_X1);free(B_Y1);free(B_Z1);free(B_XM);free(B_YM);free(B_ZM);
			free(L); free(LAMBDA);
			#ifdef CUDA
			cudaFree(X);
			#else
			free(X);
			#endif
		}
	}
all_lattice=false;
}

void Lattice::AllocateMemory(void) {
if (debug) cout <<"AllocateMemory in lattice " << endl;
	Real r, VL, LS;
	Real rlow, rhigh;
	int i {0}; int j {0}; int k {0}; //int kk {0};

	DeAllocateMemory();
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
			if (BC[4]=="mirror") {
				if (fjc==1) BZ1=1; else {
					for (k=0; k<fjc; k++) B_Z1[k]=fjc+k;
				}
			}
			if (BC[4]=="periodic") {
				if (fjc==1) BZ1=MZ; else {
					for (k=0; k<fjc; k++) B_Z1[k]=MZ+fjc-k-1;
				}
			}
			if (BC[5]=="mirror") {
				if (fjc==1) BZM=MZ; else {
					for (k=0; k<fjc; k++) B_ZM[k]=MZ+fjc-k-1;
				}
			}
			if (BC[5]=="periodic") {
				if (fjc==1) BZM=1; else {
					for (k=0; k<fjc; k++) B_ZM[k]=fjc+k;
				}
			}
			//Fall through
		case 2:
			if (BC[2]=="mirror") {
				if (fjc==1) BY1=1; else {
					for (k=0; k<fjc; k++) B_Y1[k]=fjc+k-1;
				}
			}
			if (BC[2]=="periodic") {
				if (fjc==1) BY1=MY; else {
					for (k=0; k<fjc; k++) B_Y1[k]=MY+fjc-k-1;
				}
			}
			if (BC[3]=="mirror") {
				if (fjc==1) BYM=MY; else {
					for (k=0; k<fjc; k++) B_YM[k]=MY+fjc-k;
				}
			}
			if (BC[3]=="periodic") {
				if (fjc==1) BYM=1; else {
					for (k=0; k<fjc; k++) B_YM[k]=fjc+k;
				}
			}
			//Fall through
		case 1:
			if (BC[0]=="mirror") {
				if (fjc==1) BX1=1; else {
					for (k=0; k<fjc; k++) B_X1[k]=fjc+k;
				}
			}
			if (BC[0]=="periodic") {
				if (fjc==1) BX1=MX;else {
					for (k=0; k<fjc; k++) B_X1[k]=MX+fjc-k-1;
				}
			}
			if (BC[1]=="mirror") {
				if (fjc==1) BXM=MX;else {
					for (k=0; k<fjc; k++) B_XM[k]=MX+fjc-k-1;
				}
			}
			if (BC[1]=="periodic") {
				if (fjc==1) BXM=1; else {
					for (k=0; k<fjc; k++) B_XM[k]=fjc+k;
				}
			}
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

#ifdef CUDA
	X=(Real*)AllOnDev(M);
#else
	X=(Real*)malloc(M*sizeof(Real));
#endif

	switch (gradients) {
		case 1:
			if (geometry=="planar") {
				for (i=1; i<MX+1; i++) L[i]=1;
			}

			if (fjc==1) {

				if (geometry=="cylindrical") {
					for (i=1; i<MX+1; i++) {
						r=offset_first_layer + i;
						L[i]=PIE*(pow(r,2)-pow(r-1,2));
						lambda1[i]=2.0*PIE*r/L[i]*lambda;
						lambda_1[i]=2.0*PIE*(r-1)/L[i]*lambda;
						lambda0[i]=1.0-2.0*lambda;
					}
				}
				if (geometry=="spherical") {
					for (i=1; i<MX+1; i++) {
						r=offset_first_layer + i;
						L[i]=4.0/3.0*PIE*(pow(r,3)-pow(r-1,3));
						lambda1[i]=4.0*PIE*pow(r,2)/L[i]*lambda;
						lambda_1[i]=4.0*PIE*pow(r-1,2)/L[i]*lambda;
						lambda0[i]=1.0-lambda1[i]-lambda_1[i];
					}
				}
			}

			if (fjc>1) {
				if (geometry=="planar") {
					for (i = 0; i < M; i++) {
						L[i] = 1.0/fjc;
						LAMBDA[i] = 1.0/(2*(FJC-1));
						LAMBDA[i + (FJC-1)*M] = 1.0/(2*(FJC-1));
						LAMBDA[i+(FJC-1)/2*M] = 1.0/(FJC-1);
						for (j = 1; j < FJC/2; j++) {
							LAMBDA[i+j*M] = 1.0/(FJC-1);
							LAMBDA[i+(FJC-j-1)*M] = 1.0/(FJC-1);
						}
					}
				}
				if (geometry == "cylindrical") {
					for (i = fjc; i < M - fjc; i++) {
						r = offset_first_layer+1.0*(i-fjc+1.0)/fjc;
						rlow = r - 0.5;
						rhigh = r + 0.5;
						L[i] = PIE * (2.0 * r) / fjc;
						VL = L[i] / PIE * fjc;
						if ((rlow - r) * 2 + r > 0.0)
							LAMBDA[i] += 1.0/(1.0*FJC-1.0)*rlow/VL;
						if ((rhigh - r) * 2 + r < 1.0*MX/fjc)
							LAMBDA[i+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
						else {
							if (2*rhigh-r-1.0*MX/fjc > -0.001 && 2 * rhigh-r-1.0*MX/fjc < 0.001) {
								LAMBDA[i+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
							}
							for (j = 1; j <= fjc; j++) {
								if (2*rhigh-r-1.0*MX/fjc > 0.99*j/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*j/fjc) {
									LAMBDA[i+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*(rhigh-1.0*j/fjc)/VL;
								}
							}
						}
						for (j = 1; j < fjc; j++) {
							rlow += 0.5/(fjc);
							rhigh -= 0.5/(fjc);
							if ((rlow-r)*2+r > 0.0)
								LAMBDA[i+j*M] += 1.0/(1.0*FJC-1.0)*2.0*rlow/VL;
							if ((rhigh-r)*2+r < offset_first_layer+1.0*MX/fjc)
								LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
							else {
								if (2 * rhigh-r-1.0*MX/fjc > -0.001 && 2*rhigh-r-1.0*MX/fjc < 0.001) {
									LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
								}
								for (k = 1; k <= fjc; k++) {
									if (2 * rhigh-r-1.0*MX/fjc > 0.99*k/fjc && 2*rhigh-r-1.0*MX/fjc<1.01*k/fjc) {
										LAMBDA[i + (FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*(rhigh-1.0*k/fjc)/VL;
									}
								}
							}
						}
						LS = 0;
						for (j = 0; j < FJC; j++)
							LS += LAMBDA[i+j*M];
						LAMBDA[i+(FJC/2)*M] += 1.0 - LS;
					}
				}

				if (geometry == "spherical") {
					for (i = fjc; i < M - fjc; i++) {
						r = offset_first_layer+1.0*(1.0*i-1.0*fjc+1.0)/fjc;
						rlow = r-0.5;
						rhigh = r+0.5;
						L[i] = PIE*4.0/3.0*(rhigh*rhigh*rhigh-rlow*rlow*rlow)/fjc;
						VL = L[i] / PIE * fjc;
						if ((rlow-r)*2+r > 0.0)
							LAMBDA[i] += 0.5/(1.0*FJC-1.0)*4.0*rlow*rlow/VL;
						if ((rhigh -r)*2+r < 1.0*MX/fjc)
							LAMBDA[i+(FJC-1)*M] += 0.5/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
						else {
							if (2*rhigh-r-1.0*MX/fjc>-0.001 && 2*rhigh-r-1.0*MX/fjc<0.001) {
								LAMBDA[i+(FJC-1)*M] += 0.5/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
							}
							for (j = 1; j <= fjc; j++) {
								if (2*rhigh-r-1.0*MX/fjc > 0.99*j/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*j/fjc) {
									LAMBDA[i+(FJC-1)*M] += 0.5/(1.0*FJC-1.0)*4.0*(rhigh-1.0*j/fjc)*(rhigh-1.0*j/fjc)/VL;
								}
							}
						}
						for (j = 1; j < fjc; j++) {
							rlow += 0.5/(fjc);
							rhigh -= 0.5/(fjc);
							if ((rlow-r)*2+r > 0.0)
								LAMBDA[i+j*M] += 1.0/(1.0*FJC-1.0)*4.0*rlow*rlow/VL;
							if ((rhigh - r) * 2 + r < offset_first_layer + 1.0*MX/fjc)
								LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
							else {
								if (2*rhigh-r-1.0*MX/fjc > -0.001 && 2*rhigh-r-1.0*MX/fjc < 0.001) {
									LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*4.0*rhigh*rhigh/VL;
								}
								for (k = 1; k <= fjc; k++) {
									if (2*rhigh-r-1.0*MX/fjc > 0.99*k/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*k/fjc) {
										LAMBDA[i+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*4.0*(rhigh-1.0*k/fjc)*(rhigh-1.0*k/fjc)/VL;
									}
								}
							}
						}
						LS = 0;
						for (j = 0; j < FJC; j++)
							LS += LAMBDA[i+j*M];
						LAMBDA[i+(FJC/2)*M] += 1.0-LS;
					}
				}
			}
			break;
		case 2:
			if (fjc==1) {
				if (geometry=="planar") {
					for (i=0; i<M; i++) L[i]=1;
				}
				if (geometry=="cylindrical") {
					for (int x=1; x<MX+1; x++)
					for (int y=1; y<MY+1; y++) {
						r=offset_first_layer + 1.0*x;
						L[P(x,y)]=PIE*(pow(r,2)-pow(r-1,2));
						lambda1[P(x,y)]=2.0*PIE*r/L[P(x,y)]*lambda;
						lambda_1[P(x,y)]=2.0*PIE*(r-1)/L[P(x,y)]*lambda;
						lambda0[P(x,y)]=1.0-2.0*lambda;
					}
				}
			}
			if (fjc==2) {
				if (geometry=="planar") {
					for (i=0; i<M; i++) L[i]=1.0/fjc;
				}
				if (geometry == "cylindrical") {
					for (int y=1-fjc; y<MY+fjc; y++) {
						for (int x = fjc; x < MX+fjc; x++) {
							r = offset_first_layer+1.0*(x-fjc+1.0)/fjc;
							rlow = r - 0.5;
							rhigh = r + 0.5;
							L[P(x,y)] = PIE * (2.0 * r) / fjc;
							VL = L[P(x,y)] / PIE * fjc;
							if ((rlow - r) * 2 + r > 0.0) {
								LAMBDA[P(x,y)] += 1.0/(1.0*FJC-1.0)*rlow/VL;
							}
							if ((rhigh - r) * 2 + r < 1.0*MX/fjc) {
								LAMBDA[P(x,y)+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
							} else {
								if (2*rhigh-r-1.0*MX/fjc > -0.001 && 2 * rhigh-r-1.0*MX/fjc < 0.001) {
									LAMBDA[P(x,y)+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*rhigh/VL;
								}
								for (j = 1; j <= fjc; j++) {
									if (2*rhigh-r-1.0*MX/fjc > 0.99*j/fjc && 2*rhigh-r-1.0*MX/fjc < 1.01*j/fjc) {
										LAMBDA[P(x,y)+(FJC-1)*M] += 1.0/(1.0*FJC-1.0)*(rhigh-1.0*j/fjc)/VL;
									}
								}								}
							for (j = 1; j < fjc; j++) {
								rlow += 0.5/(fjc);
								rhigh -= 0.5/(fjc);
								if ((rlow-r)*2+r > 0.0)
								LAMBDA[P(x,y)+j*M] += 1.0/(1.0*FJC-1.0)*2.0*rlow/VL;
								if ((rhigh-r)*2+r < offset_first_layer+1.0*MX/fjc)
								LAMBDA[P(x,y)+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
								else {
									if (2 * rhigh-r-1.0*MX/fjc > -0.001 && 2*rhigh-r-1.0*MX/fjc < 0.001) {
										LAMBDA[P(x,y)+(FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*rhigh/VL;
									}
									for (k = 1; k <= fjc; k++) {
										if (2 * rhigh-r-1.0*MX/fjc > 0.99*k/fjc && 2*rhigh-r-1.0*MX/fjc<1.01*k/fjc) {
											LAMBDA[P(x,y) + (FJC-1-j)*M] += 1.0/(1.0*FJC-1.0)*2.0*(rhigh-1.0*k/fjc)/VL;
										}
									}
								}
							}
							LS = 0;
							for (j = 0; j < FJC; j++)
							LS += LAMBDA[P(x,y)+j*M];
							LAMBDA[P(x,y)+(FJC/2)*M] += 1.0 - LS;
						}
					}
				}
			}
			for (k=0; k<FJC; k++) for (int y=1-fjc; y<MY+fjc; y++) for (int x=MX+1; x< MX+fjc; x++)
				LAMBDA[P(x,y)+k*M]=LAMBDA[P(2*MX-x+1,y)+(FJC-k-1)*M];
			break;
		case 3:
			break;
		default:
			break;
	}

	all_lattice=(gradients<3 && geometry!="planar");
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

bool Lattice::PutM() {
	bool success=true;

	switch(gradients) {
		case 1:
			JX=fjc; JY=0; JZ=0; M=MX+2*fjc;
			if (geometry=="planar") {volume = MX; }
			if (geometry=="spherical") {volume = 4/3*PIE*(pow(MX+offset_first_layer,3)-pow(offset_first_layer,3));}
			if (geometry=="cylindrical") {volume = PIE*(pow(MX+offset_first_layer,2)-pow(offset_first_layer,2));}
			break;
		case 2:
			if (geometry=="cylindrical")
				volume = MY*PIE*(pow(MX+offset_first_layer,2)-pow(offset_first_layer,2));
			else volume = MX*MY;
			JX=MY+2*fjc; JY=1; JZ=0; M=(MX+2*fjc)*(MY+2*fjc);
			break;
		case 3:
			volume = MX*MY*MZ;
			JX=(MZ+2*fjc)*(MY+2*fjc); JY=MZ+2*fjc; JZ=1; M = (MX+2*fjc)*(MY+2*fjc)*(MZ+2*fjc);
			break;
		default:
			break;
	}

	Accesible_volume=volume;
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


		lattice_type="simple_cubic"; lambda=1.0/6.0; Z=6;

		options.push_back("simple_cubic"); options.push_back("hexagonal");
		Value=GetValue("lattice_type");
		if (Value.length()>0) {
			if (!In[0]->Get_string(Value,lattice_type,options,"Input for 'lattice_type' not recognized. 'simple_cubic', 'FCC' or 'hexagonal'.")) success = false; else {
				if (lattice_type == "simple_cubic") {lambda=1.0/6.0; Z=6;}
				//if (lattice_type == "FCC") {lambda = 1.0/3.0; Z=3;}
				if (lattice_type == "hexagonal") {lambda=1.0/4.0; Z=4;}
			}
		} else {
			success=false; cout <<"Namics can not run without input for 'lattice_type'" << endl;
		}
/*			else {
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
*/
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
				//if (lattice_type == "") {
				//	success=false;
				//	cout <<" in two gradient calculations, you should set lattice type to either 'simple_cubic' or 'FCC'"<<endl;
				//}
				//if (Z==0) {
				//	lattice_type  = "FCC";
				//	lambda=1.0/3.0; Z=3;
				//	cout <<"Correction: in two-gradient case the default of the lattice is 'FCC'. " << endl;
				//}
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
			if (!In[0]->Get_int(GetValue("FJC_choices"),FJC,"FJC_choices can adopt only few integer values: 3 + i*2, with i = 0, 1, 2, 3, ..."))
				success=false;
			else if ((FJC-3) %2 != 0) {
				cout << "FJC_choices can adopt only few integer values: 3 + i*2, with i = 0, 1, 2, 3, ...." <<endl;
				success=false;
			}

			fjc=(FJC-3)/2+1;
		}
		if ((fjc>1) && (lattice_type !="hexagonal")) {success = false; cout << "For FJC-choices >3, we need lattice_type = 'hexagonal'." << endl; }
		if (gradients ==2 && fjc>2) {success = false; cout <<" When gradients is 2, FJC-choices are limited to 5 " << endl; }
		if (gradients ==3 && fjc>1) {success = false; cout <<" When graidents is 3, FJC-choices are limited to 3 " << endl; }

		if (GetValue("ignore_site_fraction").length()>0) ignore_sites=In[0]->Get_bool(GetValue("ignore_sites"),false);
		if (GetValue("stencil_full").length()>0) {
			stencil_full=In[0]->Get_bool(GetValue("stencil_full"),false);
			if (gradients<3 && !stencil_full) cout << "In calculations with 'gradients' less than 3, stencil_full is set to true. " << endl;
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
if (debug) cout << "TimesL in lattice " << endl;
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
if (debug) cout << "DivL in lattice " << endl;
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
if (debug) cout << "Moment in lattice " << endl;
	Real Result=0;
	Real cor;
	if (gradients !=1 || geometry!="planar" ) {cout << "Moments analysis is only implemented for one-gradient, planar system " << endl;  return Result; }
	remove_bounds(X);
	for (int i = fjc; i<M; i++) {
		cor = (i-fjc+0.5)/fjc; Result += pow(cor,n)*X[i];
	}
	return Result/fjc;
}

Real Lattice::WeightedSum(Real* X){
if (debug) cout << "weighted sum in lattice " << endl;
	Real sum{0};
	remove_bounds(X);
	switch(gradients) {
		case 1:
			if (geometry=="planar") {
				Sum(sum,X,M); sum/=fjc;
			} else Dot(sum,X,L,M);
			break;
		case 2:
			if (geometry=="planar") {
				Sum(sum,X,M); sum = sum/(fjc*fjc);
			} else	{
				Dot(sum,X,L,M);
			}
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

void Lattice::vtk(string filename, Real* X, string id,bool writebounds) {
if (debug) cout << "vtk in lattice " << endl;
	FILE *fp;
	int i;
	fp = fopen(filename.c_str(),"w+");
	switch(gradients) {
		case 1:
			cout << "for system with one gradient there is no VTK output available " << endl;
			break;
		case 2:
			fprintf(fp,"# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i %i\n",MX,MY,1);
			if (writebounds) {
				fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",(MX+2*fjc)*(MY+2*fjc));
			} else {
				fprintf(fp,"SPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA %i\n",MX*MY);
			}
			fprintf(fp,"SCALARS %s double\nLOOKUP_TABLE default\n",id.c_str());
			if (writebounds) for (i=0; i<M; i++) fprintf(fp,"%f\n",X[i]);
			for (int x=1; x<MX+1; x++)
			for (int y=1; y<MY+1; y++)
			fprintf(fp,"%f\n",X[P(x,y)]);
			break;
		case 3:
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
				fprintf(fp,"%f\n",X[P(x,y,z)]);
			}
			break;
		default:
			break;
	}
	fclose(fp);
}
void Lattice::PutProfiles(FILE* pf,vector<Real*> X,bool writebounds){
if (debug) cout <<"PutProfiles in lattice " << endl;
	int x,y,z,i;
	int length=X.size();
	int a;
	if (writebounds) a=fjc; else a = 0;
	switch(gradients) {
		case 1:
			for (x=1-a; x<MX+1+a; x++){
				fprintf(pf,"%e\t",offset_first_layer+1.0*x/fjc-1/(2.0*fjc));
				for (i=0; i<length; i++) fprintf(pf,"%.20g\t",X[i][fjc-1+x]);
				fprintf(pf,"\n");
			}
	    break;
		case 2:
			for (x=1-a; x<MX+1+a; x++)
			for (y=1-a; y<MY+1+a; y++){
				fprintf(pf,"%e\t%e\t",offset_first_layer+1.0*x/fjc-1/(2.0*fjc),1.0*y/fjc-1/(2.0*fjc));
			  for (i=0; i<length; i++) fprintf(pf,"%.20g\t",X[i][P(x,y)]);
				fprintf(pf,"\n");
			}
			break;
		case 3:
			for (x=1-a; x<MX+1+a; x++)
			for (y=1-a; y<MY+1+a; y++)
			for (z=1-a; z<MZ+1+a; z++) {
				fprintf(pf,"%e\t%e\t%e\t",1.0*x*fjc-1/(2.0*fjc),1.0*y*fjc-1/(2.0*fjc),1.0*z*fjc-1/(2.0*fjc));
				for (i=0; i<length; i++) fprintf(pf,"%.20g\t",X[i][x*JX+y*JY+fjc-1+z]);
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
	if (ignore_sites) {
		Cp(X_side,X,M); return;
	}
	Zero(X_side,M);set_bounds(X);
	int j, kk;
	switch(gradients) {
		case 1:
			if (fjc==1) {
				if (geometry=="planar") {
					if (lattice_type=="simple_cubic") {
						YplusisCtimesX(X_side+1,X,1.0/6.0,M-1);
						YplusisCtimesX(X_side,X+1,1.0/6.0,M-1);
						YplusisCtimesX(X_side,X,4.0/6.0,M);
					} else {
						YplusisCtimesX(X_side+1,X,1.0/4.0,M-1);
						YplusisCtimesX(X_side,X+1,1.0/4.0,M-1);
						YplusisCtimesX(X_side,X,2.0/4.0,M);
					}
				} else {
					AddTimes(X_side,X,lambda0,M);
					AddTimes(X_side+1,X,lambda_1+1,M-1);
					AddTimes(X_side,X+1,lambda1,M-1);
				}
			} else {

				for (j = 0; j < FJC/2; j++) {
					kk = (FJC-1)/2-j;
					AddTimes(X_side+kk, X, LAMBDA+j*M+kk, M-kk);
					AddTimes(X_side, X+kk, LAMBDA+(FJC-j-1)*M, M-kk);
				}
				AddTimes(X_side, X, LAMBDA+(FJC-1)/2*M, M);

			}
			break;

		case 2:
			if (fjc==1) {
				if (geometry=="planar") {
					if (lattice_type=="simple_cubic") {
						YplusisCtimesX(X_side,X,    16.0/36.0,M);
						YplusisCtimesX(X_side+1,X,   4.0/36.0,M-1);
						YplusisCtimesX(X_side,X+1,   4.0/36.0,M-1);
						YplusisCtimesX(X_side+JX,X,  4.0/36.0,M-JX);
						YplusisCtimesX(X_side,X+JX,  4.0/36.0,M-JX);
						YplusisCtimesX(X_side+JX+1,X,1.0/36.0,M-JX-1);
						YplusisCtimesX(X_side+JX,X+1,1.0/36.0,M-JX);
						YplusisCtimesX(X_side+1,X+JX,1.0/36.0,M-JX);
						YplusisCtimesX(X_side,X+JX+1,1.0/36.0,M-JX-1);
					} else {
						YplusisCtimesX(X_side,X,    12.0/8.0,M);
						YplusisCtimesX(X_side+1,X,   6.0/8.0,M-1);
						YplusisCtimesX(X_side,X+1,   6.0/8.0,M-1);
						YplusisCtimesX(X_side+JX,X,  6.0/8.0,M-JX);
						YplusisCtimesX(X_side,X+JX,  6.0/8.0,M-JX);
						YplusisCtimesX(X_side+JX+1,X,3.0/8.0,M-JX-1);
						YplusisCtimesX(X_side+JX,X+1,3.0/8.0,M-JX);
						YplusisCtimesX(X_side+1,X+JX,3.0/8.0,M-JX);
						YplusisCtimesX(X_side,X+JX+1,3.0/8.0,M-JX-1);
					}
				} else {
					if (lattice_type =="simple_cubic") {
						YplusisCtimesX(X_side,X,4.0/6.0,M);
						AddTimes(X_side+JX,X,lambda_1+JX,M-JX);
						AddTimes(X_side,X+JX,lambda1,M-JX);
						YplusisCtimesX(X_side+1,X,1.0/6.0,M-1);
						YplusisCtimesX(X_side,X+1,1.0/6.0,M-1);
						Norm(X_side,4.0,M);
						AddTimes(X_side+JX+1,X,lambda_1+JX+1,M-JX-1);
						AddTimes(X_side+JX,X+1,lambda_1+JX,M-JX);
						AddTimes(X_side+1,X+JX,lambda1+1,M-JX);
						AddTimes(X_side,X+JX+1,lambda1,M-JX-1);
					} else {
						YplusisCtimesX(X_side,X,2.0/4.0,M);
						AddTimes(X_side+JX,X,lambda_1+JX,M-JX);
						AddTimes(X_side,X+JX,lambda1,M-JX);
						YplusisCtimesX(X_side+1,X,1.0/4.0,M-1);
						YplusisCtimesX(X_side,X+1,1.0/4.0,M-1);
						Norm(X_side,2.0,M);
						AddTimes(X_side+JX+1,X,lambda_1+JX+1,M-JX-1);
						AddTimes(X_side+JX,X+1,lambda_1+JX,M-JX);
						AddTimes(X_side+1,X+JX,lambda1+1,M-JX);
						AddTimes(X_side,X+JX+1,lambda1,M-JX-1);
						Norm(X_side,6.0,M);
					}
				}
			}
			if (fjc==2) {
				if (geometry=="planar") {
					Add(X_side,X,M);
					Add(X_side+JX,X,M-JX);
					Add(X_side,X+JX,M-JX);
					Add(X_side+1,X,M-1);
					Add(X_side,X+1,M-1);
					Add(X_side+JX+1,X,M-JX-1);
					Add(X_side,X+JX+1,M-JX-1);
					Add(X_side+1,X+JX,M-JX);
					Add(X_side+JX,X+1,M-JX);
					Norm(X_side,2.0,M);
					Add(X_side+2*JX,X,M-2*JX);
					Add(X_side,X+2*JX,M-2*JX);
					Add(X_side+2*JX,X+1,M-2*JX);
					Add(X_side+2*JX+1,X,M-2*JX-1);
					Add(X_side+1,X+2*JX,M-2*JX);
					Add(X_side,X+2*JX+1,M-2*JX-1);
					Add(X_side+JX+2,X,M-JX-2);
					Add(X_side+JX,X+2,M-JX);
					Add(X_side,X+JX+2,M-JX-2);
					Add(X_side+2,X+JX,M-JX);
					Add(X_side+2,X,M-2);
					Add(X_side,X+2,M-2);
					Norm(X_side,2.0,M);
					Add(X_side+2*JX+2,X,M-2*JX-2);
					Add(X_side,X+2*JX+2,M-2*JX-2);
					Add(X_side+2*JX,X+2,M-2*JX);
					Add(X_side+2,X+2*JX,M-2*JX);
					Norm(X_side,1.0/64.0,M);
				} else {
					Add(X_side,X,M);
					Add(X_side+1,X,M-1);
					Add(X_side,X+1,M-1);
					Norm(X_side,1.0/4.0,M);
					AddTimes(X_side+JX,X,LAMBDA+M+JX,M-JX);
					AddTimes(X_side+JX+1,X,LAMBDA+M+JX+1,M-JX-1);
					AddTimes(X_side+JX,X+1,LAMBDA+M+JX+1,M-JX);
					AddTimes(X_side,X+JX,LAMBDA+3*M,M-JX);
					AddTimes(X_side,X+JX+1,LAMBDA+3*M,M-JX-1);
					AddTimes(X_side+1,X+JX,LAMBDA+3*M,M-JX);
					Norm(X_side,2.0,M);
					AddTimes(X_side+2*JX,X,LAMBDA+2*JX,M-2*JX);
					AddTimes(X_side,X+2*JX,LAMBDA+4*M,M-2*JX);
					AddTimes(X_side+2*JX,X+1,LAMBDA+2*JX,M-2*JX);
					AddTimes(X_side+2*JX+1,X,LAMBDA+2*JX+1,M-2*JX-1);
					AddTimes(X_side+1,X+2*JX,LAMBDA+4*M+1,M-2*JX);
					AddTimes(X_side,X+2*JX+1,LAMBDA+4*M,M-2*JX-1);
					AddTimes(X_side+JX+2,X,LAMBDA+M+JX+2,M-JX-2);
					AddTimes(X_side+JX,X+2,LAMBDA+M+JX,M-JX);
					AddTimes(X_side,X+JX+2,LAMBDA+3*M,M-JX-2);
					AddTimes(X_side+2,X+JX,LAMBDA+3*M+2,M-JX);
					AddTimes(X_side+2,X,LAMBDA+2*M+2,M-2);
					AddTimes(X_side,X+2,LAMBDA+2*M,M-2);
					Norm(X_side,2.0,M);
					AddTimes(X_side+2*JX+2,X,LAMBDA+2*JX+2,M-2*JX-2);
					AddTimes(X_side,X+2*JX+2,LAMBDA+4*M,M-2*JX-2);
					AddTimes(X_side+2*JX,X+2,LAMBDA+2*JX,M-2*JX);
					AddTimes(X_side+2,X+2*JX,LAMBDA+4*M+2,M-2*JX);
					Norm(X_side,1.0/16.0,M);
				}
			}
			break;
		case 3:
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

			break;
		default:
			break;
	}
}


void Lattice::propagate(Real *G, Real *G1, int s_from, int s_to,int M) { //this procedure should function on simple cubic lattice.
if (debug) cout <<" propagate in lattice " << endl;
	Real *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	int JX_=JX, JY_=JY; int kk;
	int k=sub_box_on;
	int j;

	Zero(gs,M); set_bounds(gs_1);
	switch(gradients) {
		case 1:
			if (fjc==1) {
				if (geometry=="planar") {
					if (lattice_type=="simple_cubic") {
						YplusisCtimesX(gs+1,gs_1,1.0/6.0,M-1);
						YplusisCtimesX(gs,gs_1,4.0/6.0,M);
						YplusisCtimesX(gs,gs_1+1,1.0/6.0,M-1);
						Times(gs,gs,G1,M);
					} else {
						YplusisCtimesX(gs+1,gs_1,1.0/4.0,M-1);
						YplusisCtimesX(gs,gs_1,2.0/4.0,M);
						YplusisCtimesX(gs,gs_1+1,1.0/4.0,M-1);
						Times(gs,gs,G1,M);
					}
				} else { //cylindrical or spherical
					AddTimes(gs,gs_1,lambda0,M);
					AddTimes(gs+1,gs_1,lambda_1+1,M-1);
					AddTimes(gs,gs_1+1,lambda1,M-1);
					Times(gs,gs,G1,M);
				}
			} else {
				for (j = 0; j < FJC/2; j++) {
					kk = (FJC-1)/2-j;
					AddTimes(gs+kk, gs_1, LAMBDA+j*M+kk, M-kk);
					AddTimes(gs, gs_1+kk, LAMBDA+(FJC-j-1)*M, M-kk);
				}
				AddTimes(gs, gs_1, LAMBDA+(FJC-1)/2*M, M);
				Times(gs, gs, G1, M);
			}
			break;
		case 2:
			if (fjc==1) {
				if (geometry=="planar") {
					if (lattice_type=="simple_cubic") { //9 point stencil
						YplusisCtimesX(gs,gs_1,    16.0/36.0,M);
						YplusisCtimesX(gs+1,gs_1,   4.0/36.0,M-1);
						YplusisCtimesX(gs,gs_1+1,   4.0/36.0,M-1);
						YplusisCtimesX(gs+JX,gs_1,  4.0/36.0,M-JX);
						YplusisCtimesX(gs,gs_1+JX,  4.0/36.0,M-JX);
						YplusisCtimesX(gs+JX+1,gs_1,1.0/36.0,M-JX-1);
						YplusisCtimesX(gs+JX,gs_1+1,1.0/36.0,M-JX);
						YplusisCtimesX(gs+1,gs_1+JX,1.0/36.0,M-JX);
						YplusisCtimesX(gs,gs_1+JX+1,1.0/36.0,M-JX-1);
						Times(gs,gs,G1,M);
					} else { //hexagonal //9 point stencil
						YplusisCtimesX(gs,gs_1,    12.0/8.0,M);
						YplusisCtimesX(gs+1,gs_1,   6.0/8.0,M-1);
						YplusisCtimesX(gs,gs_1+1,   6.0/8.0,M-1);
						YplusisCtimesX(gs+JX,gs_1,  6.0/8.0,M-JX);
						YplusisCtimesX(gs,gs_1+JX,  6.0/8.0,M-JX);
						YplusisCtimesX(gs+JX+1,gs_1,3.0/8.0,M-JX-1);
						YplusisCtimesX(gs+JX,gs_1+1,3.0/8.0,M-JX);
						YplusisCtimesX(gs+1,gs_1+JX,3.0/8.0,M-JX);
						YplusisCtimesX(gs,gs_1+JX+1,3.0/8.0,M-JX-1);
						Times(gs,gs,G1,M);
					}
				} else { //geometry==cylindrical. use \lambda's.
					if (lattice_type=="simple cubic") {
						YplusisCtimesX(gs,gs_1,4.0/6.0,M);
						AddTimes(gs+JX,gs_1,lambda_1+JX,M-JX);
						AddTimes(gs,gs_1+JX,lambda1,M-JX);
						YplusisCtimesX(gs+1,gs_1,1.0/6.0,M-1);
						YplusisCtimesX(gs,gs_1+1,1.0/6.0,M-1);
						Norm(gs,4.0,M);
						AddTimes(gs+JX+1,gs_1,lambda_1+JX+1,M-JX-1);
						AddTimes(gs+JX,gs_1+1,lambda_1+JX,M-JX);
						AddTimes(gs+1,gs_1+JX,lambda1+1,M-JX);
						AddTimes(gs,gs_1+JX+1,lambda1,M-JX-1);
						Times(gs,gs,G1,M);
					} else {
						YplusisCtimesX(gs,gs_1,2.0/4.0,M);
						AddTimes(gs+JX,gs_1,lambda_1+JX,M-JX);
						AddTimes(gs,gs_1+JX,lambda1,M-JX);
						YplusisCtimesX(gs+1,gs_1,1.0/4.0,M-1);
						YplusisCtimesX(gs,gs_1+1,1.0/4.0,M-1);
						Norm(gs,2.0,M);
						AddTimes(gs+JX+1,gs_1,lambda_1+JX+1,M-JX-1);
						AddTimes(gs+JX,gs_1+1,lambda_1+JX,M-JX);
						AddTimes(gs+1,gs_1+JX,lambda1+1,M-JX);
						AddTimes(gs,gs_1+JX+1,lambda1,M-JX-1);
						Norm(gs,6.0,M);
						Times(gs,gs,G1,M);
					}
				}
			}
			if (fjc==2) { //25 point stencil only fjc==2 implemented....
				if (geometry=="planar") { //lattice_type = hexagonal
					Add(gs,gs_1,M);
					Add(gs+JX,gs_1,M-JX);
					Add(gs,gs_1+JX,M-JX);
					Add(gs+1,gs_1,M-1);
					Add(gs,gs_1+1,M-1);
					Add(gs+JX+1,gs_1,M-JX-1);
					Add(gs,gs_1+JX+1,M-JX-1);
					Add(gs+1,gs_1+JX,M-JX);
					Add(gs+JX,gs_1+1,M-JX);
					Norm(gs,2.0,M);
					Add(gs+2*JX,gs_1,M-2*JX);
					Add(gs,gs_1+2*JX,M-2*JX);
					Add(gs+2*JX,gs_1+1,M-2*JX);
					Add(gs+2*JX+1,gs_1,M-2*JX-1);
					Add(gs+1,gs_1+2*JX,M-2*JX);
					Add(gs,gs_1+2*JX+1,M-2*JX-1);
					Add(gs+JX+2,gs_1,M-JX-2);
					Add(gs+JX,gs_1+2,M-JX);
					Add(gs,gs_1+JX+2,M-JX-2);
					Add(gs+2,gs_1+JX,M-JX);
					Add(gs+2,gs_1,M-2);
					Add(gs,gs_1+2,M-2);
					Norm(gs,2.0,M);
					Add(gs+2*JX+2,gs_1,M-2*JX-2);
					Add(gs,gs_1+2*JX+2,M-2*JX-2);
					Add(gs+2*JX,gs_1+2,M-2*JX);
					Add(gs+2,gs_1+2*JX,M-2*JX);
					Norm(gs,1.0/64.0,M);
					Times(gs,gs,G1,M);
				} else { //lattice_type hexagonal.; cylindrical geometry.
					Add(gs,gs_1,M);
					Add(gs+1,gs_1,M-1);
					Add(gs,gs_1+1,M-1);
					Norm(gs,1.0/4.0,M);
					AddTimes(gs+JX,gs_1,LAMBDA+M+JX,M-JX);
					AddTimes(gs+JX+1,gs_1,LAMBDA+M+JX+1,M-JX-1);
					AddTimes(gs+JX,gs_1+1,LAMBDA+M+JX,M-JX);
					AddTimes(gs,gs_1+JX,LAMBDA+3*M,M-JX);
					AddTimes(gs,gs_1+JX+1,LAMBDA+3*M,M-JX-1);
					AddTimes(gs+1,gs_1+JX,LAMBDA+3*M+1,M-JX);
					Norm(gs,2.0,M);
					AddTimes(gs+2*JX,gs_1,LAMBDA+2*JX,M-2*JX);
					AddTimes(gs,gs_1+2*JX,LAMBDA+4*M,M-2*JX);
					AddTimes(gs+2*JX,gs_1+1,LAMBDA+2*JX,M-2*JX);
					AddTimes(gs+2*JX+1,gs_1,LAMBDA+2*JX+1,M-2*JX-1);
					AddTimes(gs+1,gs_1+2*JX,LAMBDA+4*M+1,M-2*JX);
					AddTimes(gs,gs_1+2*JX+1,LAMBDA+4*M,M-2*JX-1);
					AddTimes(gs+JX+2,gs_1,LAMBDA+1*M+JX+2,M-JX-2);
					AddTimes(gs+JX,gs_1+2,LAMBDA+1*M+JX,M-JX);
					AddTimes(gs,gs_1+JX+2,LAMBDA+3*M,M-JX-2);
					AddTimes(gs+2,gs_1+JX,LAMBDA+3*M+2,M-JX);
					AddTimes(gs+2,gs_1,LAMBDA+2*M+2,M-2);
					AddTimes(gs,gs_1+2,LAMBDA+2*M,M-2);
					Norm(gs,2.0,M);
					AddTimes(gs+2*JX+2,gs_1,LAMBDA+2*JX+2,M-2*JX-2);
					AddTimes(gs,gs_1+2*JX+2,LAMBDA+4*M,M-2*JX-2);
					AddTimes(gs+2*JX,gs_1+2,LAMBDA+2*JX,M-2*JX);
					AddTimes(gs+2,gs_1+2*JX,LAMBDA+4*M+2,M-2*JX);
					Norm(gs,1.0/16.0,M);
					Times(gs,gs,G1,M);
				}
			}
			break;
		case 3:
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
if (debug) cout <<"remove_bounds in lattice " << endl;
	int x,y,z;
	int k;
	switch(gradients) {
		case 1:
			if (fjc==1) {
				X[0]=0;
				X[MX+1]=0;
			} else {
				for (k=0; k<fjc; k++) {
					X[(fjc-1)-k]=0;
					X[MX+fjc+k]=0;
				}
			}
			break;
		case 2:
			if (fjc==1) {
				for (x=0; x<MX+1; x++) {
					//if (P(x,0) <0 || P(x,0)>M) cout <<"P(x,0) " <<P(x,0) << endl;
					//if (P(x,MY+1) <0 || P(x,MY+1)>M) cout <<"P(x,MY+1) " <<P(x,MY+1) << endl;
					X[P(x,0)] = 0;
					X[P(x,MY+1)]=0;
				}
				for (y=0; y<MY+1; y++) {
					//if (P(0,y) <0 || P(0,y)>M) cout <<"P(0,y) " <<P(0,y) << endl;
					//if (P(MX+1,y) <0 ||P(MX+1,y)>M) cout <<"P(MX+1,y) " <<P(MX+1,y) << endl;
					X[P(0,y)] = 0;
					X[P(MX+1,y)]=0;
				}
			} else {
				for (x=1-fjc; x<MX+fjc; x++) {
					for (k=0; k<fjc; k++) X[P(x,-k)] =0;
					for (k=0; k<fjc; k++) X[P(x,MY+1+k)]=0;
				}
				for (y=1-fjc; y<MY+fjc; y++) {
					for (k=0; k<fjc; k++) X[P(-k,y)] = 0;
					for (k=0; k<fjc; k++) X[P(MX+1+k,y)]=0;
				}
			}
			break;
		case 3:
			if (sub_box_on!=0) {
				int k=sub_box_on;
				for (int i=0; i<n_box[k]; i++)
					RemoveBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
			} else
				if (fjc==1) RemoveBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ); else {
					for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++){
						for (k=0; k<fjc; k++) X[x*JX+y*JY+(fjc-1)-k] = 0;
						for (k=0; k<fjc; k++) X[x*JX+y*JY+MZ+fjc-k]  = 0;
					}
					for (y=fjc; y<MY+fjc; y++) for (z=fjc; z<MZ+fjc; z++)  {
						for (k=0; k<fjc; k++) X[(fjc-k-1)*JX+y*JY+z*JZ] = 0;
						for (k=0; k<fjc; k++) X[(MX+fjc-k)*JX+y*JY+z*JZ] = 0;
					}
					for (z=fjc; z<MZ+fjc; z++) for (x=fjc; x<MX+fjc; x++){
						for (k=0; k<fjc; k++) X[x*JX+(fjc-k-1)*JY+z*JZ] = 0;
						for (k=0; k<fjc; k++) X[x*JX+(MY+fjc-k)*JY+z*JZ] = 0;
					}
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
if (debug) cout <<"set_bounds in lattice " << endl;
	int x,y,z;
	int k=0;
	switch(gradients) {
		case 1:
			if (fjc==1) {
				X[0]=X[BX1];
				X[MX+1]=X[BXM];
			} else {
				for (k=0; k<fjc; k++) {
					X[(fjc-1)-k]=X[B_X1[k]];
					X[MX+fjc+k]=X[B_XM[k]];
				}
			}
			break;
		case 2:
			if (fjc==1) {
				for (x=1; x<MX+1; x++) {
					X[x*JX+0] = X[x*JX+BY1];
					X[x*JX+MY+1]=X[x*JX+BYM];
				}
				for (y=1; y<MY+1; y++) {
					X[0+y] = X[BX1*JX+y];
					X[(MX+1)*JX+y]=X[BXM*JX+y];
				}
				//corners
				for (x=0; x<1; x++) {
					X[x*JX+0] = X[x*JX+1];
					X[x*JX+MY+1]=X[x*JX+MY];
				}
				for (x=MX+1; x<MX+2; x++) {
					X[x*JX+0] = X[x*JX+1];
					X[x*JX+MY+1]=X[x*JX+MY];
				}
			} else {
				for (x=1; x<MX+1; x++) {
					for (k=0; k<fjc; k++) {
						//if (P(x,-k)<0 || P(x,-k)>M) cout << "P(x,-k) =" << P(x,-k) << endl;
						//if (P(x,1+k)<0 || P(x,1+k)>M) cout << "P(x,1+k) =" << P(x,1+k) << endl;
						X[P(x,-k)] = X[P(x,1+k)];
					}
					for (k=0; k<fjc; k++) {
						//if (P(x,MY+k)<0 || P(x,MY+k)>M) cout << "P(x,MY+k) =" << P(x,MY+k) << endl;
						//if (P(x,MY+1-k)<0 || P(x,MY+1-k)>M) cout << "P(x,MY+1-k) =" << P(x,MY+1-k) << endl;
						X[P(x,MY+1+k)]=X[P(x,MY-k)];
					}
				}
				for (y=1; y<MY+1; y++) {
					for (k=0; k<fjc; k++) {
						//if (P(-k,y)<0 || P(-k,y)>M) cout << "P(-k,y) =" << P(-k,y) << endl;
						//if (P(1+k,y)<0 || P(1+k,y)>M) cout << "P(1+k,y) =" << P(1+k,y) << endl;
						X[P(-k,y)] = X[P(1+k,y)];
					}
					for (k=0; k<fjc; k++) {
						//if (P(MX+1+k,y)<0 || P(MX+1+k,y)>M) cout << "P(MX+1+k,y) =" << P(MX+1+k,y) << endl;
						//if (P(MX-k,y)<0 || P(MX-k,y)>M) cout << "P(MX-k,y) =" << P(MX-k,y) << endl;
						X[P(MX+1+k,y)]=X[P(MX-k,y)];
					}
				}
				//corners
				for (x=1-fjc; x<1; x++) {
					for (k=0; k<fjc; k++) {
						//if (P(x,-k)<0 || P(x,-k)>M) cout << "1. P(x,-k) =" << P(x,-k) << "x" << x << "k" << k << endl;//
						//if (P(x,1+k)<0 || P(x,1+k)>M) cout << "P(x,1+k) =" << P(x,1+k) << endl;
						X[P(x,-k)] = X[P(x,1+k)];
					}
					for (k=0; k<fjc; k++) {
						//if (P(x,MY+1+k)<0 || P(x,MY+1+k)>M) cout << "P(x,MY+1+k) =" << P(x,MY+1+k) << endl;
						//if (P(x,MY-k)<0 || P(x,MY-k)>M) cout << "P(x,MY-k) =" << P(x,MY-k) << endl;
						X[P(x,MY+1+k)]=X[P(x,MY-k)];
					}
				}
				for (x=MX+1; x<MX+1+fjc; x++) {
					for (k=0; k<fjc; k++) {
						//if (P(x,-k)<0 || P(x,-k)>M) cout << "2. P(x,-k) =" << P(x,-k) << endl;
						//if (P(x,1+k)<0 || P(x,1+k)>M) cout << "P(x,1+k) =" << P(x,1+k) << endl;
						X[P(x,-k)] = X[P(x,1+k)];
					}
					for (k=0; k<fjc; k++) {
						//if (P(x,MY+1+k)<0 || P(x,MY+1+k)>M) cout << "P(x,MY+1+k)) =" << P(x,MY+1+k) << endl;
						//if (P(x,MY-k)<0 || P(x,MY-k)>M) cout << "P(x,MY-k) =" << P(x,MY-k) << endl;
						X[P(x,MY+1+k)]=X[P(x,MY-k)];
					}
				}
			}
			break;
		case 3:
			if (sub_box_on!=0) {
				int k=sub_box_on;
				for (int i=0; i<n_box[k]; i++)
					SetBoundaries(X+i*m[k],jx[k],jy[k],1,mx[k],1,my[k],1,mz[k],mx[k],my[k],mz[k]);
			} else
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
						H_p[i]=px*JX+py*JY+fjc-1+pz;
						break;
					case 2:
						if (xyz.size()!=2) {
							cout << "In mon " + seg_name+ " pinned_range  the expected 'pair of coordinate' -with brackets- structure '(x,y)' was not found. " << endl;  success = false;
						} else {
							px=In[0]->Get_int(xyz[0],0);
							if (px < 1 || px > MX) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the x-coordinate in pinned_range out of bounds: 0.." << MX+1 << endl; success =false;}
							py=In[0]->Get_int(xyz[1],0);
							if (py < 1 || py > MY) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 0.." << MY+1 << endl; success =false;}
						}
						H_p[i]=px*JX+fjc-1+py;
						break;
					case 1:
						if (xyz.size()!=1) {
							cout << "In mon " + seg_name+ " pinned_range  the expected 'single coordinate' -with brackets- structure '(x)' was not found. " << endl;  success = false;
						} else {
							px=In[0]->Get_int(xyz[0],0);
							if (px < 1 || px > MX) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the x-coordinate in pinned_range out of bounds: 0.." << MX+1 << endl; success =false;}
						}
						H_p[i]=fjc-1+px;
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
						if (In[0]->Get_int(lines[x-1],0)==1) {H_p[p_i]=fjc-1+x; p_i++;}
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
						H_p[i]=fjc-1+px;
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
						if (In[0]->Get_int(lines[i],0)==1) {H_p[p_i]=P(x,y); p_i++;}
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
						cout <<"reading px " << px << " and py " << py << endl;
						H_p[i]=P(px,py);
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
			break;
		default:
			break;
	}
	return success;
}

bool Lattice::FillMask(int* Mask, vector<int>px, vector<int>py, vector<int>pz, string filename) {
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
	switch (gradients) {
		case 1:
			if (readfile) {
				if (MX!=length) {success=false; cout <<"inputfile for filling delta_range has not expected length in x-direction" << endl;
				} else {
					for (int x=1; x<MX+1; x++) Mask[x]=In[0]->Get_int(lines[x],-1);
				}
			} else  {
				for (int i=0; i<length_px; i++) {
					p=px[i]; if (p<1 || p>MX) {success=false; cout<<" x-value in delta_range out of bounds; " << endl; }
					else Mask[fjc-1+px[i]]=1;
				}
			}
			break;
		case 2:
			if (readfile) {
				if (MX*MY!=length) {success=false; cout <<"inputfile for filling delta_range has not expected length in x,y-directions" << endl;
				} else {
					for (int x=1; x<MX+1; x++)
					for (int y=1; y<MY+1; y++) Mask[x*JX + y]=In[0]->Get_int(lines[x*JX+y],-1);
				}
			} else  {
				for (int i=0; i<length_px; i++) {
					p=px[i]; if (p<1 || p>MX) {success=false; cout<<" x-value in delta_range out of bounds; " << endl; }
					p=py[i]; if (p<1 || p>MY) {success=false; cout<<" y-value in delta_range out of bounds; " << endl; }
					if (success) Mask[P(px[i],py[i])]=1; //Mask[px[i]*JX + fjc-1+ py[i]]=1;
				}
			}

			break;
		case 3:
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
			break;
		default:
			break;
	}
	for (int i=0; i<M; i++) if (!(Mask[i]==0 || Mask[i]==1)) {success =false; cout <<"Delta_range does not contain '0' or '1' values. Check delta_inputfile values"<<endl; }
	return success;
}

bool Lattice::CreateMASK(int* H_MASK, int* r, int* H_P, int n_pos, bool block) {
if (debug) cout <<"CreateMask for lattice " + name << endl;
	bool success=true;
	H_Zero(H_MASK,M);
	if (block) {
		switch(gradients) {
			case 1:
				for (int x=r[0]; x<r[3]+1; x++)
						H_MASK[fjc-1+x]=1;
			break;
			case 2:
				for (int x=r[0]; x<r[3]+1; x++)
				for (int y=r[1]; y<r[4]+1; y++){
						H_MASK[P(x,y)]=1;
					}
						//H_MASK[x*JX+fjc-1+y]=1;
			break;
			case 3:
				for (int x=r[0]; x<r[3]+1; x++)
				for (int y=r[1]; y<r[4]+1; y++)
				for (int z=r[2]; z<r[5]+1; z++)
							H_MASK[x*JX+y*JY+fjc-1+z]=1;
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
 	if (fjc>1) cout <<"ReadGuess is not  yet prepared for working with FJC-choices > 3 " << endl;
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
 	if (fjc>1) cout <<"StoreGuess is not  yet prepared for working with FJC-choices > 3 " << endl;
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



void Lattice::UpdateEE(Real* EE, Real* psi, Real* E) {
	Real pf=0.5*eps0*bond_length/k_BT*(k_BT/e)*(k_BT/e); //(k_BT/e) is to convert dimensionless psi to real psi; 0.5 is needed in weighting factor.
	set_bounds(psi);
	Zero(EE,M);
	Real Exmin,Explus,Eymin,Eyplus;
	int x,y,z;
	int r;

	switch (gradients) {
		case 1:
			if (geometry =="planar" ) {
				pf =pf/2.0*fjc*fjc;
				Explus=psi[fjc-1]-psi[fjc]; //van Male approach.
				Explus *=Explus;

				for (x=fjc; x<MX+fjc; x++) {
					Exmin=Explus;
					Explus=psi[x]-psi[x+1];
					Explus *=Explus;
					EE[x]=pf*(Exmin+Explus);
				}
			}
			if (geometry=="cylindrical" ) {
				r=offset_first_layer*fjc;
				pf=pf*PIE;
				for (x=fjc; x<MX+fjc; x++) {
					r++;
					Exmin=psi[x]-psi[x-1];
					Exmin*=(r-1)*Exmin;
					Explus=psi[x]-psi[x+1];
					Explus*=(r)*Explus;
					EE[x]=pf*(Exmin+Explus)/L[x];
				}


/*
				pf=pf*PIE;
				r=offset_first_layer*fjc-1.5;
				if (offset_first_layer <=0) {
					r++;
					Explus = psi[fjc+1]-psi[fjc+2];
					Explus *=(r+0.5)*Explus;
					EE[fjc+1]=pf*Explus/L[fjc+1];
					x=fjc+2;
				} else {
					Explus=psi[fjc]-psi[fjc+1];
					Explus *=(r+0.5)*Explus;
					x=fjc+1;
				}

				for (; x<MX+fjc; x++) {
					r +=1.0;
					Exmin=Explus;
					Explus=(psi[x]-psi[x+1]);
					Explus *= (r+0.5)*Explus;
					EE[x]=pf*(Exmin+Explus)/L[x];
				}
*/

			}
			if (geometry=="spherical" ) {
				pf=pf*PIE*2/fjc;
				r=offset_first_layer*fjc +1.0;
				Explus=r*(psi[fjc]-psi[fjc+1]);
				Explus *=Explus;
				EE[fjc]=pf*Explus/(L[fjc]);
				for (x=fjc+1; x<MX+fjc; x++) {
					r +=1.0;
					Exmin=Explus;
					Explus=r*(psi[x]-psi[x+1]);
					Explus *=Explus;
					EE[x]=pf*(Exmin+Explus)/L[x];
				}

			}
			break;
		case 2:
			if (geometry=="planar") {
				pf = pf/2.0;
				for (x=fjc; x<MX+fjc; x++) {
					for (y=fjc; y<MY+fjc; y++) {
						z=x*JX+y;
						Exmin=psi[z]-psi[z-JX];
						Exmin*=Exmin;
						Explus=psi[z]-psi[z+JX];
						Explus*=Explus;
						Eymin=psi[z]-psi[z-1];
						Eymin*=Eymin;
						Eyplus=psi[z]-psi[z+1];
						Eyplus*=Eyplus;
						EE[x*MX+y]=pf*(Exmin+Explus+Eymin+Eyplus);
					}
				}
			} else {
				pf = pf/2;
				r=offset_first_layer*fjc;
				for (x=fjc; x<MX+fjc; x++) {
					r++;
					for (y=fjc; y<MY+fjc; y++) {
						z=x*JX+y;
						Exmin=psi[z]-psi[z-JX];
						Exmin*=(r-0.5)*Exmin*2*PIE;
						Explus=psi[z]-psi[z+JX];
						Explus*=(r+0.5)*Explus*2*PIE;
						Eymin=(psi[z]-psi[z-1])*fjc;
						Eymin*=Eymin;
						Eyplus=(psi[z]-psi[z+1])*fjc;
						Eyplus*=Eyplus;
						EE[z]=pf*((Exmin+Explus)/L[z]+Eymin+Eyplus);
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
			break;
*/
		default:
			break;
	}
}


void Lattice::UpdatePsi(Real* g, Real* psi ,Real* q, Real* eps, int* Mask, bool grad_epsilon, bool fixedPsi0) { //not only update psi but also g (from newton).
	int x,y,i;
#ifndef CUDA
	int z;
#endif

#ifndef CUDA
	Real epsZplus, epsZmin;
#endif
	Real a,b,c,a_,b_,c_;
	Real r;
	Real epsXplus, epsXmin, epsYplus,epsYmin;
	set_bounds(eps);
	Real C =e*e/(eps0*k_BT*bond_length);

   if (!fixedPsi0) {
	switch (gradients) {
		case 1:

			if (geometry=="planar") { //van Male approach;
				C=C*2.0/fjc/fjc;
				epsXplus=eps[fjc-1]+eps[fjc];
				a=0; b=psi[fjc-1]; c=psi[fjc];
				for (x=fjc; x<MX+fjc; x++) {
					epsXmin=epsXplus;
					epsXplus=eps[x]+eps[x+1];
					a=b; b=c; c=psi[x+1];
					X[x]=(epsXmin*a  +C*q[x] + epsXplus*c)/(epsXmin+epsXplus);
				 }
			}
			if (geometry=="cylindrical") {
				C=C/PIE;
				r=offset_first_layer*fjc;
				epsXplus=r*(eps[fjc-1]+eps[fjc]);
				a=0; b=psi[fjc-1]; c=psi[fjc];
				for (x=fjc; x<MX+fjc; x++) {
					r++;
					epsXmin=epsXplus;
					epsXplus=r*(eps[x]+eps[x+1]);
					a=b; b=c; c=psi[x+1];
					X[x]=(epsXmin*a + C*q[x]*L[x] + epsXplus*c)/(epsXmin+epsXplus);
				 }
			}
			if (geometry=="spherical") {
				C=C/(2.0*PIE)*fjc;
				r=offset_first_layer*fjc;
				epsXplus=r*r*(eps[fjc-1]+eps[fjc]);
				a=0; b=psi[fjc-1]; c=psi[fjc];
				for (x=fjc; x<MX+fjc; x++) {
					epsXmin=epsXplus;
					r++;
					epsXplus=r*r*(eps[x]+eps[x+1]);
					a=b; b=c; c=psi[x+1];
					X[x]=(epsXmin*a + C*q[x]*L[x] + epsXplus*c)/(epsXmin+epsXplus);
				 }
			}
			//Cp(psi,X,M);
			YisAminB(g,g,X,M);
			break;
		case 2:
			if (geometry=="planar") {
				C=C*2.0/fjc/fjc;
				for (x=fjc; x<MX+fjc; x++) {
					for (y=fjc; y<MY+fjc; y++) {
						i=x*JX+y;
						epsXmin=eps[i]+eps[i-JX];
						epsXplus=eps[i]+eps[i+JX];
						epsYmin=eps[i]+eps[i-1];
						epsYplus=eps[i]+eps[i+1];
						X[i]= (C*q[i]+epsXmin*psi[i-JX]+epsXplus*psi[i+JX]+epsYmin*psi[i-1]+epsYplus*psi[i+1])/
						(epsXmin+epsXplus+epsYmin+epsYplus);
					}
				}
			} else {
				C=C/2/fjc/fjc;
				r=offset_first_layer*fjc;
				for (x=fjc; x<MX+fjc; x++) {
					r++;
					for (y=fjc; y<MY+fjc; y++) {
						i=x*JX+y;
						epsXmin=2*PIE*(r-1)*(eps[i]+eps[i-JX])/L[i]*fjc*fjc;
						epsXplus=2*PIE*r*(eps[i]+eps[i+JX])/L[i]*fjc*fjc;
						epsYmin=eps[i]+eps[i-1];
						epsYplus=eps[i]+eps[i+1];
						X[i]= (C*q[i]+epsXmin*psi[i-JX]+epsXplus*psi[i+JX]+epsYmin*psi[i-1]+epsYplus*psi[i+1])/
						(epsXmin+epsXplus+epsYmin+epsYplus);
					}
				}

			}
			//Cp(psi,X,M);
			YisAminB(g,g,X,M);
			break;
		case 3:
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
			break;
		default:
			break;
	} //end of switch
   } else { //fixedPsi0 is true
	switch (gradients) {
		case 1:
			a=0; b=psi[fjc-1]; c=psi[fjc];
			for (x=fjc; x<MX+fjc; x++) {
				a=b; b=c; c=psi[x+1];
				if (Mask[x] == 0) psi[x]=0.5*(a+c)+q[x]*C/eps[x];
			}

			if (geometry=="cylindrical") {
				a=0; b=psi[fjc-1]; c=psi[fjc];
				for (x=fjc; x<MX+fjc; x++) {
					a=b; b=c; c=psi[x+1];
					if (Mask[x] == 0) psi[x]+=(c-a)/(2.0*(offset_first_layer*fjc+x-fjc+0.5))*fjc;
				}
			}
			if (geometry=="spherial") {
				a=0; b=psi[fjc-1]; c=psi[fjc];
				for (x=fjc; x<MX+fjc; x++) {
					a=b; b=c; c=psi[x+1];
					if (Mask[x] == 0) psi[x]+=(c-a)/(offset_first_layer*fjc+x-fjc+0.5)*fjc;
				}
			}
			if (grad_epsilon) {
				a=0; b=psi[fjc-1]; c=psi[fjc];a_=0; b_=eps[fjc-1]; c_=eps[fjc];
				for (x=fjc; x<MX+fjc; x++) {//for all geometries
					a=b; b=c; c=psi[x+1]; a_=b_; b_=c_; c_=eps[x+1];
					if (Mask[x] == 0) {
						psi[x]+=0.25*(c_-a_)*(c-a)/eps[x]*fjc*fjc;

					}
				}
			}
			for (x=fjc; x<MX+fjc; x++)
			if (Mask[x] == 0) {
				g[x]-=psi[x];
			}
			break;
		case 2:
			for (x=fjc; x<MX+fjc; x++) {
				for (y=fjc; y<MY+fjc; y++){
					if (Mask[x*JX+y] == 0)
					X[x*JX+y]=0.25*(psi[(x-1)*JX+y]+psi[(x+1)*JX+y]
						        +psi[x*JX+y-1]  +psi[x*JX+y+1])
							 +0.5*q[x*JX+y]*C/eps[x*JX+y];
				}
			}

			if (geometry=="cylindrical") { //radial geometry in x no radial geometry in y
				for (x=fjc; x<MX+fjc; x++) {
					for (y=fjc; y<MY+fjc; y++){
						if (Mask[x*JX+y] == 0)
						X[x*JX+y]+=(psi[(x+1)*JX+y]-psi[(x-1)*JX+y])/(2.0*(offset_first_layer+x-fjc+0.5))*fjc;
					}
				}
			}
			if (grad_epsilon) {
				for (x=fjc; x<MX+fjc; x++) {
					for (y=fjc; y<MY+fjc; y++) {
						if (Mask[x*JX+y] == 0) {
							X[x*JX+y]+=0.25*(eps[(x+1)*JX+y]-eps[(x-1)*JX+y])*(psi[(x+1)*JX+y]-psi[(x-1)*JX+y]+
                                              		   	   eps[x*JX+y+1]  -eps[x*JX+y-1])  *(psi[x*JX+y+1]  -psi[x*JX+y-1])/
					           		   	   eps[x*JX+y]*fjc*fjc;
						}
					}
				}
			}
			for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)
			if (Mask[x*JX+y] == 0) {
				psi[x*JX+y]=X[x*JX+y];
				g[x*JX+y]-=psi[x*JX+y];
			}
			break;
		case 3:
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
			break;
		default:
			break;
	} //end of switch

   } //end of if then else

}


void Lattice::UpdateQ(Real* g, Real* psi, Real* q, Real* eps, int* Mask,bool grad_epsilon) {//Not only update q (charge), but also g (from newton).
	int x,y;
	#ifndef CUDA
	int z;
	#endif
	Real a,b,c,a_,b_,c_;
	Real epsXplus,epsXmin,epsYplus,epsYmin,epsZplus,epsZmin;

	Real C = -e*e/(eps0*k_BT*bond_length);
	switch (gradients) {
		case 1:
			a=0; b=psi[fjc-1]; c=psi[fjc];
			for (x=fjc; x<MX+fjc; x++) { //for all geometries
				a=b; b=c; c=psi[x+1];
				if (Mask[x] == 1) q[x] = -0.5*(a-2*b+c)*fjc*fjc*eps[x]/C;
			}

			if (geometry=="cylindrical") {
				a=0; b=psi[fjc-1]; c=psi[fjc];
				for (x=fjc; x<MX+fjc; x++) {
					a=b; b=c; c=psi[x+1];
					if (Mask[x] == 1) q[x]-=(c-a)/(2.0*(offset_first_layer*fjc+x-fjc+0.5))*fjc*eps[x]/C;
				}
			}
			if (geometry=="spherial") {
				a=0; b=psi[fjc-1]; c=psi[fjc];
				for (x=fjc; x<MX+fjc; x++) {
					a=b; b=c; c=psi[x+1];
					if (Mask[x] == 1) q[x]-=(c-a)/(offset_first_layer*fjc+x-fjc+0.5)*fjc*eps[x]/C;
				}
			}
			if (grad_epsilon) {
				a=0; b=psi[fjc-1]; c=psi[fjc]; a_=0; b_=eps[fjc-1]; c_=eps[fjc];
				for (x=fjc; x<MX+fjc; x++) {//for all geometries
					a=b; b=c; c=psi[x+1]; a_=b_; b_=c_; c_=eps[x+1];
					if (Mask[x] == 1) q[x]-=0.25*(c_-a_)*(c-a)*fjc*fjc/C;
				}
			}
			for (x=fjc; x<MX+fjc; x++)
			if (Mask[x] == 1) {
				g[x]=-q[x];
			}
			break;
		case 2:
			for (x=fjc; x<MX+fjc; x++) {
				for (y=fjc; y<MY+fjc; y++){ //for all geometries
					if (Mask[x*JX+y] == 1)
					q[x*JX+y]=-0.5*(psi[(x-1)*JX+y]+psi[(x+1)*JX+y]
						        +psi[x*JX+y-1]  +psi[x*JX+y+1]
							 -4*psi[x*JX+y])*fjc*fjc*eps[x*JX+y]/C;
				}
			}

			if (geometry=="cylindrical") { //radial geometry in x no radial geometry in y
				for (x=fjc; x<MX+fjc; x++) {
					for (y=fjc; y<MY+fjc; y++){
						if (Mask[x*JX+y] == 1)
						q[x*JX+y]-=(psi[(x+1)*JX+y]-psi[(x-1)*JX+y])/(2.0*(offset_first_layer+x-fjc+0.5))*fjc*eps[x]/C;
					}
				}
			}
			if (grad_epsilon) {
				for (x=fjc; x<MX+fjc; x++) {
					for (y=fjc; y<MY+fjc; y++) {
						if (Mask[x*JX+y] == 1)
						q[x*JX+y]-=0.25*(eps[(x+1)*JX+y]-eps[(x-1)*JX+y])*(psi[(x+1)*JX+y]-psi[(x-1)*JX+y]+
                                              		   eps[x*JX+y+1]  -eps[x*JX+y-1])  *(psi[x*JX+y+1]  -psi[x*JX+y-1])*fjc*fjc/C;
						}
				}
			}
			for (x=fjc; x<MX+fjc; x++) for (y=fjc; y<MY+fjc; y++)
			if (Mask[x*JX+y] == 1) {
				g[x*JX+y]=-q[x*JX+y];
			}
			break;
		case 3:
#ifdef CUDA
			C *=6;
			Cp(X,psi,M);
		UpQ(g + JX + JY + 1, q + JX + JY + 1, psi + JX + JY + 1, eps + JX + JY + 1, JX, JY, C, Mask + JX + JY + 1, M - 2 * (JX + JY + 1));
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
			break;
		default:
			break;

	}
}
