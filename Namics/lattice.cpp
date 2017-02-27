#include "lattice.h" 
Lattice::Lattice(vector<Input*> In_,string name_) { //this file contains switch (gradients). In this way we keep all the lattice issues in one file!
if (debug) cout <<"Lattice constructor" << endl; 
	In=In_; name=name_; 
	KEYS.push_back("gradients"); KEYS.push_back("n_layers"); KEYS.push_back("offset_first_layer");
	KEYS.push_back("n_layers_x");   KEYS.push_back("n_layers_y"); KEYS.push_back("n_layers_z");
	KEYS.push_back("lowerbound"); KEYS.push_back("upperbound");
	KEYS.push_back("lowerbound_x"); KEYS.push_back("upperbound_x");
	KEYS.push_back("lowerbound_y"); KEYS.push_back("upperbound_y");
	KEYS.push_back("lowerbound_z"); KEYS.push_back("upperbound_z");
 	KEYS.push_back("bond_length");  KEYS.push_back("lattice_type");  
}

Lattice::~Lattice() {
if (debug) cout <<"lattice destructor " << endl;
  //In this program, we will assume that the propagator will work on a simple cubic lattice. 
			//Interactions will be treated either with simple cubic or FCC 'lambda's.
if (lambda_1) delete [] lambda_1;
if (lambda1) delete [] lambda1;
if (lambda0) delete [] lambda0; 
if (L) delete [] L; 
}

void Lattice::AllocateMemory(void) {
if (debug) cout <<"AllocateMemory in lattice " << endl; 
	double r;
	double one_over_lambda_min_two=1.0/lambda-2.0;
	int i,j;
	switch (gradients) {
		case 1:
			if (geometry!="planar") {
				L = new double[M]; Zero(L,M);
				lambda_1 = new double[M]; Zero(lambda_1,M);
				lambda1 = new double[M]; Zero(lambda1,M);
				lambda0 = new double[M]; Zero(lambda0,M);
			}
			if (geometry=="cylindrical") {
				for (i=1; i<MX+1; i++) {
					r=offset_first_layer + i; 
					L[i]=PIE*(pow(r,2)-pow(r-1,2));
					lambda1[i]=2*PIE*r/L[i];
					lambda_1[i]=2*PIE*(r-1)/L[i];
					lambda0[i]=one_over_lambda_min_two;
				}
			}
			if (geometry=="spherical") {
				for (i=1; i<MX+1; i++) {
					r=offset_first_layer + i; 
					L[i]=4/3*PIE*(pow(r,3)-pow(r-1,3));
					lambda1[i]=4*PIE*pow(r,2)/L[i];
					lambda_1[i]=4*PIE*pow((r-1),2)/L[i];
					lambda0[i]=1.0/lambda-lambda1[i]-lambda_1[i];
				}
			}
			break;
		case 2:
			if (geometry!="planar") {
				L=new double[M]; Zero(L,M); 
				lambda_1 = new double[M]; Zero(lambda_1,M);
				lambda1 = new double[M]; Zero(lambda1,M);
				lambda0 = new double[M]; Zero(lambda0,M);				
			}
			if (geometry=="cylindrical") {
				for (i=1; i<MX+1; i++) 
				for (j=1; j<MY+1; j++) {
					r=offset_first_layer + i; 
					L[i*JX+j]=PIE*(pow(r,2)-pow(r-1,2));
					lambda1[i*JX+j]=2*PIE*r/L[i*JX+j];
					lambda_1[i*JX+j]=2*PIE*(r-1)/L[i*JX+j];
					lambda0[i*JX+j]=one_over_lambda_min_two;
				}
			
			}
			break;
		case 3:
			break;
		default:
			break;
	}
}

bool Lattice::CheckInput(int start) {
if (debug) cout <<"CheckInput in lattice " << endl; 
	bool success;
	string Value;
	string VALUE1,VALUE2,VALUE3,VALUE4,VALUE5,VALUE6;

	vector<string> options;
	success = In[0]->CheckParameters("lat",name,start,KEYS,PARAMETERS,VALUES);
	if (success){
		bond_length =  In[0]->Get_int(GetValue("bond_length"),5e-10);
		if (bond_length < 0 || bond_length > 1e-8) {cout <<" bond_length out of range 0..1e-8 " << endl; success=false;}
		options.push_back("simple_cubic"); options.push_back("FCC");
		Value=GetValue("lattice_type"); 
		if (Value.length()>0) {
			if (!In[0]->Get_string(Value,lattice_type,options,"Input for 'lattice_type' not recognized.")) success = false; else {
				if (lattice_type == "simple_cubic") lambda=1.0/6.0;
				if (lattice_type == "FCC") lambda = 1.0/3.0; 
			} 
		} else {
			lattice_type  = "simple_cubic"; lambda=1.0/6.0;
		} 
		offset_first_layer =0;
		gradients=In[0]->Get_int(GetValue("gradients"),3);
		if (gradients<0||gradients>3) {cout << "value of gradients out of bounds 1..3; default value '1' is used instead " << endl; gradients=1;}	
		switch(gradients) {
			case 1: 
				MX = In[0]->Get_int(GetValue("n_layers"),-123);
				if (MX==-123) {success=false; cout <<"In 'lat' the parameter 'n_layers' is required. Problem terminated" << endl;} 
				else {
					if (MX<0 || MX >1e6) {success = false; cout <<"n_layers out of bounds, currently: 0..1e6; Problem terminated" << endl; }
				}
				options.clear(); 
				options.push_back("spherical");
				options.push_back("cylindrical");
				options.push_back("flat");options.push_back("planar");
				if (GetValue("geometry").size()>0) {
					if (!In[0]->Get_string(GetValue("geometry"),geometry,options,"In lattice input for 'geometry' not recognized.")) success=false; 
				} else geometry = "planar";
				if (geometry=="flat") geometry="planar"; 	
				if (geometry!="planar") {
					offset_first_layer=In[0]->Get_double(GetValue("offset_first_layer"),0);
					if (offset_first_layer<0) {
						cout <<"value of 'offset_first_layer' can not be negative. Value ignored. " << endl; 
						offset_first_layer=0;
					}
				}
				
				if (geometry=="planar") {volume = MX;}
				if (geometry=="spherical") {volume = 4/3*PIE*(pow(MX+offset_first_layer,3)-pow(offset_first_layer,3));}
				if (geometry=="cylindrical") {volume = PIE*(pow(MX+offset_first_layer,2)-pow(offset_first_layer,2));}
				
				JX=1; JY=0; M=MX+2; 
				options.clear(); 
				options.push_back("mirror"); //options.push_back("mirror_2");
				options.push_back("surface"); 
				if (GetValue("lowerbound_x").size()>0) {cout << "lowerbond_x is not allowed in 1-gradient calculations" << endl; success=false;}
				if (GetValue("lowerbound_y").size()>0) {cout << "lowerbond_y is not allowed in 1-gradient calculations" << endl; success=false;}
				if (GetValue("lowerbound_z").size()>0) {cout << "lowerbond_z is not allowed in 1-gradient calculations" << endl; success=false;}
				if (GetValue("upperbound_x").size()>0) {cout << "upperbond_x is not allowed in 1-gradient calculations" << endl; success=false;}
				if (GetValue("upperbound_y").size()>0) {cout << "upperbond_y is not allowed in 1-gradient calculations" << endl; success=false;}
				if (GetValue("upperbound_z").size()>0) {cout << "upperbond_z is not allowed in 1-gradient calculations" << endl; success=false;}
				
				Value.clear();
				Value=GetValue("lowerbound"); 
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,VALUE1,options,"for 'lowerbound' boundary condition not recognized.")) success=false; 
					if (VALUE1=="mirror") BX1=1; 
					//if (VALUE1=="mirror_2") BX1=2;
					if (VALUE1=="surface")  BX1=0;
				} else {
					VALUE1="mirror";
					BX1=1;
				} 
				Value.clear();
				Value=GetValue("upperbound"); 
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,VALUE2,options,"for 'upperbound' boundary condition not recognized.")) success=false; 
					if (VALUE2=="mirror") BXM=MX; 
					//if (VALUE2=="mirror_2") BXM=MX-1;
					if (VALUE2=="surface")  BXM=MX+1;
				} else {
					VALUE2="mirror";
					BXM=MX-1;
				} 
				
				break;
			case 2: 
				MX = In[0]->Get_int(GetValue("n_layers_x"),-123);
				if (MX==-123) {success=false; cout <<"In 'lat' the parameter 'n_layers_x' is required. Problem terminated" << endl;} 
				else {
					if (MX<0 || MX >1e6) {success = false; cout <<"n_layers_x out of bounds, currently: 0.. 1e6; Problem terminated" << endl; }
				}
				MY = In[0]->Get_int(GetValue("n_layers_y"),-123);
				if (MY==-123) {success=false; cout <<"In 'lat' the parameter 'n_layers_y' is required. Problem terminated" << endl;} 
				else {
					if (MY<0 || MY >1e6) {success = false; cout <<"n_layers_y out of bounds, currently: 0.. 1e6; Problem terminated" << endl; }
				}
				options.clear(); 
				options.push_back("cylindrical");
				options.push_back("flat");options.push_back("planar");
				if (GetValue("geometry").size()>0) {
					if (!In[0]->Get_string(GetValue("geometry"),geometry,options,"In lattice input for 'geometry' not recognized.")) success=false; 
				} else geometry = "planar";
				if (geometry=="flat") geometry="planar"; 
				if (geometry=="planar") {volume = MX*MY;}

				if (geometry!="planar") {
					offset_first_layer=In[0]->Get_double(GetValue("offset_first_layer"),0);
					if (offset_first_layer<0) {
						cout <<"value of 'offset_first_layer' can not be negative. Value ignored. " << endl; 
						offset_first_layer=0;
					}
				}

				if (geometry=="cylindrical") {volume = MY*PIE*(pow(MX+offset_first_layer,2)-pow(offset_first_layer,2));}
				JX=(MX+2); JY=1; M=(MX+2)*(MY+2); 

				options.clear(); 
				options.push_back("mirror"); //options.push_back("mirror_2"); 
				options.push_back("surface");
				if (geometry=="planar") options.push_back("periodic"); 
				if (GetValue("lowerbound_z").size()>0) {cout << "lowerbond_z is not allowed in 2-gradient calculations" << endl; success=false;}
				if (GetValue("upperbound_z").size()>0) {cout << "upperbond_z is not allowed in 2-gradient calculations" << endl; success=false;}
				
				Value.clear();
				Value=GetValue("lowerbound_x"); 
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,VALUE1,options,"for 'lowerbound_x' boundary condition not recognized.")) success=false; 
					if (VALUE1=="mirror") BX1=1; 
					//if (VALUE1=="mirror_2") BX1=2;
					if (VALUE1=="surface")  BX1=0;
					if (VALUE1=="periodic") BX1=MX;
				} else {
					VALUE1="mirror";
					BX1=1;
				} 
				Value.clear();
				Value=GetValue("upperbound_x"); 
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,VALUE2,options,"for 'upperbound_x' boundary condition not recognized.")) success=false; 
					if (VALUE2=="mirror") BXM=MX; 
					//if (VALUE2=="mirror_2") BXM=MX-1;
					if (VALUE2=="surface")  BXM=MX+1;
					if (VALUE2=="periodic") BXM=1;
				} else {
					VALUE2="mirror";
					BXM=MX-1;
				}
 
				if (geometry !="planar") options.push_back("periodic"); 
				Value.clear();
				Value=GetValue("lowerbound_y"); 
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,VALUE3,options,"for 'lowerbound_x' boundary condition not recognized.")) success=false; 
					if (VALUE3=="mirror") BY1=1; 
					//if (VALUE3=="mirror_2") BY1=2;
					if (VALUE3=="surface")  BY1=0;
					if (VALUE3=="periodic") BY1=MY; 
				} else {
					VALUE3="mirror";
					BY1=1;
				} 
				Value.clear();
				Value=GetValue("upperbound_x"); 
				if (Value.length()>0) {
					if (!In[0]->Get_string(Value,VALUE4,options,"for 'upperbound_x' boundary condition not recognized.")) success=false; 
					if (VALUE4=="mirror") BYM=MY; 
					//if (VALUE4=="mirror_2") BYM=MY-1;
					if (VALUE4=="surface")  BYM=MY+1;
					if (VALUE4=="periodic") BYM=1;
				} else {
					VALUE4="mirror";
					BYM=MY-1;
				}
				if (VALUE1=="periodic" || VALUE2=="periodic") {
					if (VALUE1!=VALUE2) {success=false;  cout <<"For boundaries in x-direction: 'periodic' BC  should be set to upper and lower bounds " << endl;} 
				}
				if (VALUE3=="periodic" || VALUE4=="periodic") {
					if (VALUE3!=VALUE4) {success=false;  cout <<"For boundaries in y-direction: 'periodic' BC should be set to upper and lower bounds " << endl;} 
				}
				break;
			case 3:
				if (!In[0]->Get_int(GetValue("n_layers_x"),MX,1,1e6,"In 'lat' the parameter 'n_layers_x' is required")) {success=false;}
				if (!In[0]->Get_int(GetValue("n_layers_y"),MY,1,1e6,"In 'lat' the parameter 'n_layers_y' is required")) {success=false;}
				if (!In[0]->Get_int(GetValue("n_layers_z"),MZ,1,1e6,"In 'lat' the parameter 'n_layers_z' is required")) {success=false;}
				volume=MX*MY*MZ;
				JX=(MX+2)*(MY+2); JY=(MY+2); M = (MX+2)*(MY+2)*(MZ+2);   
				Mx=In[0]->Get_int(GetValue("sub_box_size"),MX/2);
				if (Mx>MX) {cout << "'sub_box_size' can not exceed the size of the main box." << endl; success=false;} else My=Mz=Mx;
				
				options.clear(); 
				options.push_back("mirror"); //options.push_back("mirror_2"); 
				options.push_back("periodic"); 
				//options.push_back("shifted_mirror");
				
				Value.clear();
				Value=GetValue("lowerbound_x"); if (Value.length()>0) {
					if (!In[0]->Get_string(Value,VALUE1,options,"for 'lowerbound_x' boundary condition not recognized.")) success=false;
					if (VALUE1=="mirror") BX1=1; 
					//if (VALUE1=="mirror_2") BX1=2;
					if (VALUE1=="periodic") BX1=MX;
				} else {
					VALUE1="periodic";
					BX1=MX;
				} 
				Value.clear();	
				Value=GetValue("lowerbound_y"); if (Value.length()>0) {
					if (!In[0]->Get_string(Value,VALUE2,options,"for 'lowerbound_y' boundary condition not recognized.")) success=false;
					if (VALUE2=="mirror") BY1=1; 
					//if (VALUE2=="mirror_2") BY1=2;
					if (VALUE2=="periodic") BY1=MY; 
				} else {
					VALUE2="periodic";
					BY1=MY;
				} 
				Value.clear();	
				Value=GetValue("lowerbound_z"); if (Value.length()>0) {
					if (!In[0]->Get_string(Value,VALUE3,options,"for 'lowerbound_z' boundary condition not recognized.")) success=false;
					if (VALUE3=="mirror") BZ1=1; 
					//if (VALUE3=="mirror_2") BZ1=2;
					if (VALUE3=="periodic") BZ1=MZ;  
				} else {
					VALUE3="periodic";
					BZ1=MZ;
				} 
				Value.clear();	
				Value=GetValue("upperbound_x"); if (Value.length()>0) {
					if (!In[0]->Get_string(Value,VALUE4,options,"for 'upperbound_x' boundary condition not recognized.")) success=false; 
					if (VALUE4=="mirror") BXM=MX; 
					//if (VALUE4=="mirror_2") BXM=MX-1;
					if (VALUE4=="periodic") BXM=1;
				} else {
					VALUE4="periodic";
					BXM=1;
				} 
				Value.clear();	
				Value=GetValue("upperbound_y"); if (Value.length()>0) {
					if(!In[0]->Get_string(Value,VALUE5,options,"for 'upperbound_y' boundary condition not recognized.")) success=false; 
					if (VALUE5=="mirror") BYM=MY; 
					//if (VALUE5=="mirror_2") BYM=MY-1;
					if (VALUE5=="periodic") BYM=1;
				} else {
					VALUE5="periodic";
					BYM=1;
				} 
				Value.clear();	
				Value=GetValue("upperbound_z"); if (Value.length()>0) { 
					if (!In[0]->Get_string(Value,VALUE6,options,"for 'upperbound_z' boundary condition not recognized.")) success=false; 
					if (VALUE6=="mirror") BZM=MZ; 
					//if (VALUE6=="mirror_2") BZM=MZ-1;
					if (VALUE6=="periodic") BZM=1;
				} else {
					VALUE6="periodic";
					BZM=1;
				} 	
	
				if (VALUE1=="periodic" || VALUE4=="periodic") {
					if (VALUE1 != VALUE4) {cout <<"In x-direction the boundary conditions do not match:" + VALUE1 << " and " <<  VALUE4 << endl; success=false;}
				}
				if (VALUE2=="periodic" || VALUE5=="periodic") {
					if (VALUE2 != VALUE5) {cout <<"In y-direction the boundary conditions do not match:" + VALUE2 << " and " <<  VALUE5 << endl; success=false;}
				}
				if (VALUE3=="periodic" || VALUE6=="periodic") {
					if (VALUE3 != VALUE6) {cout <<"In z-direction the boundary conditions do not match:" + VALUE3 << " and " <<  VALUE6 << endl; success=false;}
				}

				break;
			default:
				cout << "gradients out of bounds " << endl; 
				break;
		}		
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

double Lattice::GetValue(double* X,string s){ //need a check to find out if s contains 3 integers separated by ','
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

void Lattice::TimesL(double* X){
	switch(gradients) {
		case 1:
			if (geometry!="planar") Times(X,X,L,M);
			break;
		case 2: 
			if (geometry!="planar") Times(X,X,L,M);			
			break;
		case 3: 
			break;
		default:
			break;
	}
}

void Lattice::DivL(double* X){
	switch(gradients) {
		case 1:
			if (geometry!="planar") Div(X,L,M);
			break;
		case 2: 
			if (geometry!="planar") Div(X,L,M);			
			break;
		case 3: 
			break;
		default:
			break;
	}
}

double Lattice::WeightedSum(double* X){
	double sum;
	switch(gradients) {
		case 1:
			if (geometry=="planar") {
				Sum(sum,X,M); return sum;
			} else {
				Dot(sum,X,L,M); return sum;
			}
			break;
		case 2: 
			if (geometry=="planar") {
				Sum(sum,X,M); return sum;
			} else {
				Dot(sum,X,L,M); return sum;
			}			
			break;
		case 3: 
			Sum(sum,X,M); return sum;
			break;
		default:
			return 0;
			break;
	}
}

void Lattice::vtk(string filename, double* X, string id) {
if (debug) cout << "vtk in lattice " << endl;
	FILE *fp;
	int i,j,k;
	fp = fopen(filename.c_str(),"w+");
	switch(gradients) {
		case 1:
			cout << "for system with one gradient there is no VTK output available " << endl; 
			break;
		case 2: 
			fprintf(fp,"# vtk DataFile Version 7.0 \nvtk output \nDATASET STRUCTURED_POINTS \nDIMENSIONS %i %i %i\n",MX,MY,1);
			fprintf(fp,"SPACING 1 1 1 \nORIGIN 0 0 0 \nPOINT_DATA %i\n",MX*MY);
			fprintf(fp,"SCALARS %s \nLOOKUP_TABLE default \n",id.c_str());
			for (i=1; i<MX+1; i++)
			for (j=1; j<MY+1; j++)
			fprintf(fp,"%f \n",X[i*JX+j]);
			break;
		case 3:	
			fprintf(fp,"# vtk DataFile Version 7.0 \nvtk output \nDATASET STRUCTURED_POINTS \nDIMENSIONS %i %i %i\n",MX,MY,MZ);
			fprintf(fp,"SPACING 1 1 1 \nORIGIN 0 0 0 \nPOINT_DATA %i\n",MX*MY*MZ);
			fprintf(fp,"SCALARS %s \nLOOKUP_TABLE default \n",id.c_str());
			for (i=1; i<MX+1; i++)
			for (j=1; j<MY+1; j++)
			for (k=1; k<MZ+1; k++)
			fprintf(fp,"%f \n",X[i*JX+j*JY+k]);
			break;
		default:
			break;
	}
	fclose(fp);
}
void Lattice::PutProfiles(FILE* pf,vector<double*> X){
if (debug) cout <<"PutProfiles in lattice " << endl; 
	int x,y,z,i;
	int length=X.size(); 	
	switch(gradients) {
		case 1:			
			for (x=1; x<MX+1; x++){
				fprintf(pf,"%i \t",x);
				for (i=0; i<length; i++) fprintf(pf,"%f \t",X[i][x]);
				fprintf(pf,"\n");
			}
			break;
		case 2:
			for (x=1; x<MX+1; x++)
			for (y=1; y<MY+1; y++){
				fprintf(pf,"%i \t %i \t",x,y);
				for (i=0; i<length; i++) fprintf(pf,"%f \t",X[i][x*JX+y]);
				fprintf(pf,"\n");
			}
			break;
		case 3: 
			for (x=1; x<MX+1; x++)
			for (y=1; y<MY+1; y++)
			for (z=1; z<MZ+1; z++) {
				fprintf(pf,"%i \t %i \t %i \t",x,y,z);
				for (i=0; i<length; i++) fprintf(pf,"%f \t",X[i][x*JX+y*JY+z]);
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

void Lattice::push(string s, double X) {
if (debug) cout <<"push (double) in lattice " << endl; 
	doubles.push_back(s);
	doubles_value.push_back(X); 
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
	doubles.clear();
	doubles_value.clear();
	ints.clear();
	ints_value.clear(); 
	string mirror="mirror"; 
	string periodic="periodic";
	//string mirror2="mirror_2";
	string surface="surface"; 
	//string s;
	push("gradients",gradients);
	if (offset_first_layer>0) push("offset_first_layer",offset_first_layer);

	push("volume",volume); 
	push("lattice_type",lattice_type);
	push("bond_length",bond_length); 
	switch (gradients) {
		case 1:
			push("n_layers",MX);
			if (BX1==1) push("lowerbound",mirror);
			if (BXM==MX-1) push("upperbound",mirror);
			//if (BX1==2) push("lowerbound",mirror_2);
			//if (BXM==MX-2) push("upperbound",mirror_2);
			if (BX1==0) push("lowerbound",surface);
			if (BXM==MX+1) push("upperbound",surface);
			break;
		case 2:
			push("n_layers_x",MX);
			push("n_layers_y",MY);
			if (BX1==1) push("lowerbound_x",mirror);
			if (BXM==MX-1) push("upperbound_x",mirror);
			//if (BX1==2) push("lowerbound_x",mirror_2);
			//if (BXM==MX-2) push("upperbound_x",mirror_2);
			if (BX1==0) push("lowerbound_x",surface);
			if (BXM==MX+1) push("upperbound_x",surface);
			if (BY1==1) push("lowerbound_x",mirror);
			if (BYM==MY-1) push("upperbound_x",mirror);
			//if (BY1==2) push("lowerbound_x",mirror_2);
			//if (BYM==MY-2) push("upperbound_x",mirror_2);
			if (BY1==0) push("lowerbound_x",surface);
			if (BYM==MY+1) push("upperbound_x",surface);
			break;
		case 3:	
			push("n_layers_x",MX);
			push("n_layers_y",MY);
			push("n_layers_z",MZ);	
			if (BX1==1) push("lowerbound_x",mirror);
			if (BXM==MX-1) push("upperbound_x",mirror);
			//if (BX1==2) push("lowerbound_x",mirror_2);
			//if (BXM==MX-2) push("upperbound_x",mirror_2);
			if (BX1==MX) push("lowerbound_x",periodic);
			if (BXM==1) push("upperbound_x",periodic);
			if (BY1==1) push("lowerbound_y",mirror);
			if (BYM==MY-1) push("upperbound_y",mirror);
			//if (BY1==2) push("lowerbound_y",mirror_2);
			//if (BYM==MY-2) push("upperbound_y",mirror_2);
			if (BY1==MY) push("lowerbound_y",periodic);
			if (BYM==1) push("upperbound_y",periodic);	
			if (BZ1==1) push("lowerbound_z",mirror);
			if (BZM==MZ-1) push("upperbound_z",mirror);
			//if (BZ1==2) push("lowerbound_z",mirror_2);
			//if (BZM==MZ-2) push("upperbound_z",mirror_2);
			if (BZ1==MZ) push("lowerbound_z",periodic);
			if (BZM==1) push("upperbound_z",periodic);
			break;
		default:
			break;
	} 
}

int Lattice::GetValue(string prop,int &int_result,double &double_result,string &string_result){
if (debug) cout <<"GetValue (long)  in lattice " << endl; 
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


void Lattice::Side(double *X_side, double *X, int M) { //this procedure should use the lambda's according to lattice_type;
if (debug) cout <<" Side in lattice " << endl; 
	double* SIDE;
	switch(gradients) {
		case 1:
			//set_bounds(X);
			Zero(X_side,M);
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
			break;
		case 2:
			//set_bounds(X);
			Zero(X_side,M); 
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
					Add(X_side,X,M);
					AddTimes(X_side,X,lambda0,M);
					AddTimes(X_side+JX,X,lambda_1+JX,M);
					AddTimes(X_side,X+JX,lambda1,M);
					Norm(X_side,lambda*lambda,M);
					
				} else {
					SIDE=new double[M];
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
			Zero(X_side,M); //set_bounds(X);
			Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
			Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
			Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
			Norm(X_side,1.0/6.0,M);
			break;
		default:
			break;
	}
}

void Lattice::propagate(double *G, double *G1, int s_from, int s_to) { //this procedure should function on simple cubic lattice. 
if (debug) cout <<" propagate in lattice " << endl; 
	double *gs = G+M*(s_to), *gs_1 = G+M*(s_from);
	switch(gradients) {
		case 1:
			Zero(gs,M); set_bounds(gs_1);
			if (geometry=="planar") {
				Add(gs+1,gs_1,M-1); Add(gs,gs_1+1,M-1);	
				YplusisCtimesX(gs,gs_1,4.0,M);
				Norm(gs,lambda,M); Times(gs,gs,G1,M);
//cout <<" from - to " << s_from << " - " << s_to << endl; 
			} else {
				AddTimes(gs,gs_1,lambda0,M);
				AddTimes(gs+1,gs_1,lambda_1+1,M-1);
				AddTimes(gs,gs_1+1,lambda1,M-1);
				Norm(gs,lambda,M); Times(gs,gs,G1,M);
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
				Add(gs,gs_1,M);
				AddTimes(gs,gs_1,lambda0,M);
				AddTimes(gs+JX,gs_1,lambda_1+JX,M);
				AddTimes(gs,gs_1+JX,lambda1,M);
				Norm(gs,lambda*lambda,M);Times(gs,gs,G1,M);
			}
			break;
		case 3:
			Zero(gs,M); set_bounds(gs_1);
			Add(gs+JX,gs_1,M-JX); Add(gs,gs_1+JX,M-JX);
			Add(gs+JY,gs_1,M-JY); Add(gs,gs_1+JY,M-JY);
			Add(gs+1,gs_1,M-1);  Add(gs,gs_1+1, M-1);
			Norm(gs,1.0/6.0,M); Times(gs,gs,G1,M);
			break;
		default:	
			break; 			
	}
}

void Lattice::remove_bounds(double *X){ 
if (debug) cout <<" remove_bounds (double) in lattice " << endl; 
	int x,y;
	switch(gradients) {
		case 1:
			X[0]=0; X[MX+1]=0;
			break;
		case 2: 
			for (x=1; x<MX+1; x++) {
				X[x*JX+0] = 0;
				X[x*JX+MY+1]=0;
			}
			for (y=1; y<MY+1; x++) {
				X[0+y] = 0;
				X[(MX+1)*JX+y]=0;
			}
			X[0]             =0;
			X[(MX+1)*JX]     =0;
			X[0+MY+1]        =0;
			X[(MX+1)*JX+MY+1]=0;	
			break;
		case 3: 
			RemoveBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
			break;
		default:
			
			break;
	}
}

void Lattice::remove_bounds(int *X){ 
if (debug) cout <<" remove_bounds (int) in lattice " << endl; 
	int x,y;
	switch(gradients) {
		case 1:
			X[0]=0; X[MX+1]=0;
			break;
		case 2: 
			for (x=1; x<MX+1; x++) {
				X[x*JX+0] =0;
				X[x*JX+MY+1]=0;
			}
			for (y=1; y<MY+1; x++) {
				X[0+y] = 0;
				X[(MX+1)*JX+y]=0;
			}
			X[0]             =0;
			X[(MX+1)*JX]     =0;
			X[0+MY+1]        =0;
			X[(MX+1)*JX+MY+1]=0;			

			break;
		case 3:
			RemoveBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
			break;
		default:
			break;
	}
}
 
void Lattice::set_bounds(double *X){  
if (debug) cout <<"set_bounds (doubles) in lattice " << endl; 
	int x,y; 
	switch(gradients) {
		case 1:
			X[0]=X[BX1]; X[MX+1]=X[BXM];
			break;
		case 2: 
			for (x=1; x<MX+1; x++) {
				X[x*JX+0] = X[x*JX+BY1];
				X[x*JX+MY+1]=X[x*JX+BYM];
			}
			for (y=1; y<MY+1; x++) {
				X[0+y] = X[BX1*JX+y];
				X[(MX+1)*JX+y]=X[BXM*JX+y];
			}
			X[0]             =X[1]*X[JX]; if (X[0]>0) X[0]=sqrt(X[0]);
			X[(MX+1)*JX]     =X[MX*JX]*X[(MX+1)*JX+1]; if (X[(MX+1)*JX]>0) X[(MX+1)*JX]=sqrt(X[(MX+1)*JX]);
			X[0+MY+1]        =X[JX+MY+1]*X[0+MY]; if (X[MY+1]>0) X[MY+1]=sqrt(X[MY+1]);
			X[(MX+1)*JX+MY+1]=X[MX*JX+MY+1]*X[(MX+1)*JX+MY]; if (X[(MX+1)*JX+MY+1]>0) X[(MX+1)*JX+MY+1]=sqrt(X[(MX+1)*JX+MY+1]);						
			break;
		case 3:
			SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
			break;
		default:
			break;
	}
}

void Lattice::set_bounds(int *X){  
if (debug) cout <<"set_bounds (int) in lattice " << endl; 
	int x,y; 	
	switch(gradients) {
		case 1:
			X[0]=X[BX1]; X[MX+1]=X[BXM];
			break;
		case 2: 
			for (x=1; x<MX+1; x++) {
				X[x*JX+0] = X[x*JX+BY1];
				X[x*JX+MY+1]=X[x*JX+BYM];
			}
			for (y=1; y<MY+1; x++) {
				X[0+y] = X[BX1*JX+y];
				X[(MX+1)*JX+y]=X[BXM*JX+y];
			}
			X[0]             =X[1]*X[JX]; if (X[0]>0) X[0]=sqrt(X[0]);
			X[(MX+1)*JX]     =X[MX*JX]*X[(MX+1)*JX+1]; if (X[(MX+1)*JX]>0) X[(MX+1)*JX]=sqrt(X[(MX+1)*JX]);
			X[0+MY+1]        =X[JX+MY+1]*X[0+MY]; if (X[MY+1]>0) X[MY+1]=sqrt(X[MY+1]);
			X[(MX+1)*JX+MY+1]=X[MX*JX+MY+1]*X[(MX+1)*JX+MY]; if (X[(MX+1)*JX+MY+1]>0) X[(MX+1)*JX+MY+1]=sqrt(X[(MX+1)*JX+MY+1]);			
			break;
		case 3: 
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
	int a; if (range_type=="frozen_range") a=1; else a=0;
	In[0]->split(range,';',set);
	if (set.size()==2) { 
		coor.clear(); 
		block=true; In[0]->split(set[0],',',coor);		
		switch(gradients) {
			case 1:
				if (coor.size()!=1) {cout << "In mon " + 	seg_name + ", for 'pos 1', in '" + range_type + "' the coordiantes do not come as a single coordinate 'x'" << endl; success=false;}
				else {
					diggit=coor[0].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[0]=In[0]->Get_int(coor[0],0); else {
						recognize_keyword=false; 
						if (coor[0]=="firstlayer") {recognize_keyword=true; r[0] = 1;}
						if (coor[0]=="lowerbound") {recognize_keyword=true; r[0] = 0;}
						if (coor[0]=="upperbound") {recognize_keyword=true; r[0] = MX+1;}
						if (coor[0]=="lastlayer")  {recognize_keyword=true; r[0] = MX;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
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
						if (coor[0]=="lowerbound") {recognize_keyword=true; r[3] = 0;}
						if (coor[0]=="upperbound") {recognize_keyword=true; r[3] = MX+1;}
						if (coor[0]=="lastlayer")  {recognize_keyword=true; r[3] = MX;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
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
						if (coor[0]=="lowerbound") {recognize_keyword=true; r[0] = 0;}
						if (coor[0]=="upperbound") {recognize_keyword=true; r[0] = MX+1;}
						if (coor[0]=="lastlayer")  {recognize_keyword=true; r[0] = MX;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
							success=false; 
						}
					}
					if (r[0] < 1-a || r[0] > MX+a) {cout << "In mon " + seg_name + ", for 'pos 1', the x-coordinate in '" + range_type + "' is out of bounds: "<< 1-a <<" .." << MX+a << endl; success =false;}
					diggit=coor[1].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[1]=In[0]->Get_int(coor[1],0); else {
						recognize_keyword=false; 
						if (coor[1]=="firstlayer") {recognize_keyword=true; r[1] = 1;}
						if (coor[1]=="lowerbound") {recognize_keyword=true; r[1] = 0;}
						if (coor[1]=="upperbound") {recognize_keyword=true; r[1] = MY+1;}
						if (coor[1]=="lastlayer")  {recognize_keyword=true; r[1] = MY;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
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
						if (coor[0]=="lowerbound") {recognize_keyword=true; r[3] = 0;}
						if (coor[0]=="upperbound") {recognize_keyword=true; r[3] = MX+1;}
						if (coor[0]=="lastlayer")  {recognize_keyword=true; r[3] = MX;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
							success=false;
						}
					}
					if (r[3] <1-a || r[3] > MX+a) {cout << "In mon " + seg_name+ ", for 'pos 2', the x-coordinate in '" + range_type + "' is out of bounds: " << 1-a << " .. " << MX+a<< endl; success=false;} 
					diggit=coor[1].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[4]=In[0]->Get_int(coor[1],0); else {
						recognize_keyword=false;
						if (coor[1]=="firstlayer") {recognize_keyword=true; r[4] = 1;}
						if (coor[1]=="lowerbound") {recognize_keyword=true; r[4] = 0;}
						if (coor[1]=="upperbound") {recognize_keyword=true; r[4] = MY+1;}
						if (coor[1]=="lastlayer")  {recognize_keyword=true; r[4] = MY;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
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
						if (coor[0]=="lowerbound") {recognize_keyword=true; r[0] = 0;}
						if (coor[0]=="upperbound") {recognize_keyword=true; r[0] = MX+1;}
						if (coor[0]=="lastlayer")  {recognize_keyword=true; r[0] = MX;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
							success=false; 
						}
					}
					if (r[0] < 1-a || r[0] > MX+a) {cout << "In mon " + seg_name + ", for 'pos 1', the x-coordinate in '" + range_type + "' is out of bounds: "<< 1-a <<" .."<< MX+a << endl; success =false;}
					diggit=coor[1].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[1]=In[0]->Get_int(coor[1],0); else {
						recognize_keyword=false; 
						if (coor[1]=="firstlayer") {recognize_keyword=true; r[1] = 1;}
						if (coor[1]=="lowerbound") {recognize_keyword=true; r[1] = 0;}
						if (coor[1]=="upperbound") {recognize_keyword=true; r[1] = MY+1;}
						if (coor[1]=="lastlayer")  {recognize_keyword=true; r[1] = MY;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
							success=false; 
						}
					}					
					if (r[1] < 1-a || r[1] > MY+a) {cout << "In mon " + seg_name+ ", for 'pos 1', the y-coordinate in '" + range_type + "' is out of bounds: "<< 1-a <<" .." << MY+a << endl; success =false;}
					diggit=coor[2].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[2]=In[0]->Get_int(coor[2],0); else {
						recognize_keyword=false; 
						if (coor[2]=="firstlayer") {recognize_keyword=true; r[2] = 1;}
						if (coor[2]=="lowerbound") {recognize_keyword=true; r[2] = 0;}
						if (coor[2]=="upperbound") {recognize_keyword=true; r[2] = MZ+1;}
						if (coor[2]=="lastlayer")  {recognize_keyword=true; r[2] = MZ;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
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
						if (coor[0]=="lowerbound") {recognize_keyword=true; r[3] = 0;}
						if (coor[0]=="upperbound") {recognize_keyword=true; r[3] = MX+1;}
						if (coor[0]=="lastlayer")  {recognize_keyword=true; r[3] = MX;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
							success=false; 
						}
					}
					if (r[3] < 1-a || r[3] > MX+a) {cout << "In mon " + seg_name+ ", for 'pos 2', the x-coordinate in '" + range_type + "' is out of bounds; "<< 1-a <<" .."<< MX+a << endl; success =false;}
					diggit=coor[1].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[4]=In[0]->Get_int(coor[1],0); else {
						recognize_keyword=false; 
						if (coor[1]=="firstlayer") {recognize_keyword=true; r[4] = 1;}
						if (coor[1]=="lowerbound") {recognize_keyword=true; r[4] = 0;}
						if (coor[1]=="upperbound") {recognize_keyword=true; r[4] = MY+1;}
						if (coor[1]=="lastlayer")  {recognize_keyword=true; r[4] = MY;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
							success=false; 
						}
					}					
					if (r[4] < 1-a || r[4] > MY+a) {cout << "In mon " + seg_name+ ", for 'pos 2', the y-coordinate in '" + range_type + "' is out of bounds; "<< 1-a <<" .." << MY+a << endl; success =false;}
					diggit=coor[2].substr(0,1);
					if (In[0]->IsDigit(diggit)) r[5]=In[0]->Get_int(coor[2],0); else {
						recognize_keyword=false; 
						if (coor[2]=="firstlayer") {recognize_keyword=true; r[5] = 1;}
						if (coor[2]=="lowerbound") {recognize_keyword=true; r[5] = 0;}
						if (coor[2]=="upperbound") {recognize_keyword=true; r[5] = MZ+1;}
						if (coor[2]=="lastlayer")  {recognize_keyword=true; r[5] = MZ;}
						if (!recognize_keyword) {
							cout << "In mon " + seg_name + " and  range_type " + range_type + ", the first 'pos' of x-coordinate is not a number or does not contain the keywords: 'firstlayer' 'lastlayer' 'upperbound' or 'lowerbound'" << endl;
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
				if (coor[0]=="firstlayer")         {recognize_keyword=true; r[0] = r[4] = 1; }
				if (coor[0]=="lowerbound" && a==1) {recognize_keyword=true; r[0] = r[4] = 0; }
				if (coor[0]=="upperbound" && a==1) {recognize_keyword=true; r[0] = r[4] = MX+1; }
				if (coor[0]=="lastlayer")          {recognize_keyword=true; r[0] = r[4] = MX;}
				if (!recognize_keyword) {
					cout << "In mon " + seg_name + " and  range_type " + range_type + ", the input: 'firstlayer' or 'lastlayer' ";
					if (a==1) cout << "or 'upperbound' or 'lowerbound'"; cout <<" is expected. " << endl;
					success=false; 
				}
			} else {
				cout << "In mon " + seg_name + " and  range_type " + range_type + ", the info was not recognised because when 'gradients>1' the lonely keywords 'firstlayer' 'lastlayers', 'upperbound' 'lowerbound' do not work." << endl;
				success=false; 
			}
		} else {
			int length;
			int length_xyz;
			int px,py,pz;
			string s;
			int i; 
			switch(gradients) {
				case 1: 
					if (n_pos==0) {
						block=false;  
						length=coor.size(); n_pos=length; 
					} else {

						length=coor.size(); n_pos=length;
						px=0,py=0,pz=0;
						i=0;	
						while (i<length) { 
							s=coor[i].substr(1,coor[i].size()-1); 
							In[0]->split(s,',',xyz);
							length_xyz=xyz.size();
							if (length_xyz!=1) { 
								cout << "In mon " + seg_name+ " pinned_range  the expected 'single coordinate' -with brackets- structure '(x)' was not found. " << endl;  success = false;
							} else {   
								px=In[0]->Get_int(xyz[0],0);
								if (px < 0 || px > MX+1) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the x-coordinate in pinned_range out of bounds: 0.." << MX+1 << endl; success =false;}					
							}
							H_p[i]=px;
							i++;
						}
					}				
					break;
				case 2: 
					if (n_pos==0) {
						block=false;  
						length=coor.size(); n_pos=length; 
					} else {
						length=coor.size(); n_pos=length;
						px=0,py=0,pz=0;
						i=0;	
						while (i<length) { 
							s=coor[i].substr(1,coor[i].size()-1); 
							In[0]->split(s,',',xyz);
							length_xyz=xyz.size();
							if (length_xyz!=2) { 
								cout << "In mon " + seg_name+ " pinned_range  the expected 'pair of coordinate' -with brackets- structure '(x,y)' was not found. " << endl;  success = false;
							} else {   
								px=In[0]->Get_int(xyz[0],0);
								if (px < 0 || px > MX+1) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the x-coordinate in pinned_range out of bounds: 0.." << MX+1 << endl; success =false;}
								py=In[0]->Get_int(xyz[1],0);
								if (py < 0 || py > MY+1) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 0.." << MY+1 << endl; success =false;}									
							}
							H_p[i]=px*JX+py;
							i++;
						}
					}				
					break;
				case 3: 		
					if (n_pos==0) {
						block=false;  
						length=coor.size(); n_pos=length; 
					} else {
						length=coor.size(); n_pos=length;
						px=0,py=0,pz=0;
						i=0;	
						while (i<length) { 
							s=coor[i].substr(1,coor[i].size()-1); 
							In[0]->split(s,',',xyz);
							length_xyz=xyz.size();
							if (length_xyz!=3) { 
								cout << "In mon " + seg_name+ " pinned_range  the expected 'triple coordinate' -with brackets- structure '(x,y,z)' was not found. " << endl;  success = false;
							} else {   
								px=In[0]->Get_int(xyz[0],0);
								if (px < 1 || px > MX) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the x-coordinate in pinned_range out of bounds: 1.." << MX << endl; success =false;}
								py=In[0]->Get_int(xyz[1],0);
								if (py < 1 || py > MY) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 1.." << MY << endl; success =false;}								
								pz=In[0]->Get_int(xyz[2],0);
								if (pz < 1 || pz > MZ) {cout << "In mon " + seg_name+ ", for 'pos' "<< i << ", the y-coordinate in pinned_range out of bounds: 1.." << MZ << endl; success =false;}	
							}
							H_p[i]=px*JX+py*JY+pz;
							i++;
						}
					}
					break;
				default:
					break;
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
	switch(gradients) {
		case 1:
			if (!In[0]->ReadFile(sub[0].append(".").append(filename),content)) {success=false;} else {
				In[0]->split(content,'#',lines);
				length = lines.size();
				if (length == MX) { //expect to read 'mask file';
					i=0;
					if (n_pos==0) {
						while (i<length){
							if (In[0]->Get_int(lines[i],0)==1) n_pos++; 
							i++; 
						}; 
						if (n_pos==0) {cout << "Warning: Input file for locations of 'particles' does not contain any unities." << endl;}
					} else {  
						i=0; p_i=0;
						for (x=1; x<MX+1; x++) {
							if (In[0]->Get_int(lines[i],0)==1) {H_p[p_i]=x; p_i++;}
							i++; 
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
			}
			break;
		case 2: 
			if (!In[0]->ReadFile(sub[0].append(".").append(filename),content)) {success=false;} else {
				In[0]->split(content,'#',lines);
				length = lines.size();
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
			} 
			break;
		case 3: 	
			if (!In[0]->ReadFile(sub[0].append(".").append(filename),content)) {success=false;} else {
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
	int x,y,z;
	switch(gradients) {
		case 1:
			H_Zero(H_MASK,M);
			if (block) {
				for (x=r[1]; x<r[4]+1; x++)  H_MASK[x]=1;
			} else {
				for (int i=0; i<n_pos; i++) H_MASK[H_P[i]]=1; 	
			}
			break;
		case 2:
			H_Zero(H_MASK,M);
			if (block) {
				for (x=r[1]; x<r[4]+1; x++) for (y=r[2]; y<r[5]+1; y++) H_MASK[x*JX+y]=1;
			} else {
				for (int i=0; i<n_pos; i++) H_MASK[H_P[i]]=1; 	
			}
			break;
		case 3: 	
			H_Zero(H_MASK,M);
			if (block) {
				for (x=r[1]; x<r[4]+1; x++) for (y=r[2]; y<r[5]+1; y++) for (z=r[3]; z<r[6]+1; z++) H_MASK[x*JX+y*JY+z]=1;
			} else {
				for (int i=0; i<n_pos; i++) H_MASK[H_P[i]]=1; 	
			}
			break;
		default:
			break;
	}
	return success; 
}

bool Lattice::GenerateGuess(double* x, string CalculationType, string GuessType, double A_value, double B_value) {
if (debug) cout <<"GenerateGuess in lattice " << endl;
//GuessType: lamellae,Im3m,FCC,BCC,HEX,gyroid,double_gyroid,double_diamond,perforated_lamellae
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
				for (i=0; i<MX+2; i++) for (j=0; j<MY+2; j++) {
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

bool Lattice::ReadGuess(string filename,double *xx) {
cout <<"ReadGuess not yet implemented in lattice" << endl; 
	bool success=true;
	return success;
}

bool Lattice::StoreGuess(string filename,double *xx) {
cout <<"StoreGuess not yet implemented in lattice" << endl;
	bool success=true;
	return success; 
}

/* 

All below is commented out. This is stored here to recover from earlier program. 

void Sideh(double *X_side, double *X, int M) {
	Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);

	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Add(X_side+JX,X+JY,MM-JX-JY); Add(X_side+JY,X+JX,MM-JY-JX);
	Add(X_side+1,X+JY,MM-JY-1); Add(X_side+JY,X+1,MM-JY-1);
	Add(X_side+1,X+JX,MM-JX-1); Add(X_side+JX,X+1,MM-JX-1);
	Norm(X_side,1.0/12.0,M);
}



void Side(double *X_side, double *X, int M) {
	Zero(X_side,M); SetBoundaries(X,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(X_side+JX,X,M-JX); Add(X_side,X+JX,M-JX);
	Add(X_side+JY,X,M-JY); Add(X_side,X+JY,M-JY);
	Add(X_side+1,X,M-1);  Add(X_side,X+1, M-1);
	Norm(X_side,1.0/6.0,M);
}

void advanced_average(double *X_side, double *X, int M){
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

void Propagateh(double* G, double* G1, int s_from, int s_to) { //on small boxes
	int MMM=M*n_box;
	double *gs = G+MMM*s_to, *gs_1 = G+MMM*s_from, *g = G1;
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

void PROPAGATEh(double *G, double *G1, int s_from, int s_to) { //on big box
	double *gs = G+MM*(s_to), *gs_1 = G+MM*(s_from), *g = G1;
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

void Propagate(double* G, double* G1, int s_from, int s_to) { //on small boxes
	int MMM=M*n_box;
	double *gs = G+MMM*s_to, *gs_1 = G+MMM*s_from, *g = G1;
	Zero(gs,MMM);
	for (int p=0; p<n_box; p++) SetBoundaries(gs_1+M*p,jx,jy,bx1,bxm,by1,bym,bz1,bzm,Mx,My,Mz);
	Add(gs+jx,gs_1,MMM-jx); Add(gs,gs_1+jx,MMM-jx);
	Add(gs+jy,gs_1,MMM-jy); Add(gs,gs_1+jy,MMM-jy);
	Add(gs+1,gs_1,MMM-1);  Add(gs,gs_1+1, MMM-1);
	Norm(gs,1.0/6.0,MMM); Times(gs,gs,g,MMM);
}

void PROPAGATE(double *G, double *G1, int s_from, int s_to) { //on big box
	double *gs = G+MM*(s_to), *gs_1 = G+MM*(s_from), *g = G1;
	Zero(gs,MM);
	SetBoundaries(gs_1,JX,JY,BX1,BXM,BY1,BYM,BZ1,BZM,MX,MY,MZ);
	Add(gs+JX,gs_1,MM-JX); Add(gs,gs_1+JX,MM-JX);
	Add(gs+JY,gs_1,MM-JY); Add(gs,gs_1+JY,MM-JY);
	Add(gs+1,gs_1,MM-1);  Add(gs,gs_1+1, MM-1);
	Norm(gs,1.0/6.0,MM); Times(gs,gs,g,MM);
}

#ifdef CUDA
void DistributeG1(double *G1, double *g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	//int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
	//distributeg1<<<n_blocks,block_size>>>(G1,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
        dim3 blocks(ceil(Mx/8.0),ceil(My/8.0),ceil(Mz/8.0));
        dim3 blockdim(8,8,8);
        distributeg1<<<blocks,blockdim>>>(G1,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);

}
#else
void DistributeG1(double *G1, double *g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int pos_l=-M;
	int pos_x,pos_y,pos_z;
	int Bxp,Byp,Bzp;
	int ii=0,jj=0,kk=0;

	for (int p=0; p<n_box; p++) { pos_l +=M; ii=0; Bxp=H_Bx[p]; Byp=H_By[p]; Bzp=H_Bz[p];
		for (int i=1; i<Mx+1; i++) { ii+=jx; jj=0; if (Bxp+i>MX) pos_x=(Bxp+i-MX)*JX; else pos_x = (Bxp+i)*JX;
			for (int j=1; j<My+1; j++) {jj+=jy;  kk=0; if (Byp+j>MY) pos_y=(Byp+j-MY)*JY; else pos_y = (Byp+j)*JY;
				for (int k=1; k<Mz+1; k++) { kk++; if (Bzp+k>MZ) pos_z=(Bzp+k-MZ); else pos_z = (Bzp+k);
					g1[pos_l+ii+jj+kk]=G1[pos_x+pos_y+pos_z];
				}
			}
		}
	}
}
#endif

#ifdef CUDA
void CollectPhi(double* phi, double* GN, double* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
 	 dim3 blocks(ceil(Mx/8.0),ceil(My/8.0),ceil(Mz/8.0));
        dim3 blockdim(8,8,8);
        collectphi<<<blocks,blockdim>>>(phi,GN,rho,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
}
#else
void CollectPhi(double* phi, double* GN, double* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int pos_l=-M;
	int pos_x,pos_y,pos_z;
	int Bxp,Byp,Bzp;
	double Inv_H_GNp;
	int ii=0,jj=0,kk=0;

	for (int p=0; p<n_box; p++) { pos_l +=M; ii=0; Bxp=Bx[p]; Byp=By[p]; Bzp=Bz[p]; Inv_H_GNp=1.0/GN[p];
		for (int i=1; i<Mx+1; i++) {ii+=jx; jj=0;  if (Bxp+i>MX) pos_x=(Bxp+i-MX)*JX; else pos_x = (Bxp+i)*JX;
			for (int j=1; j<My+1; j++) {jj+=jy;  kk=0; if (Byp+j>MY) pos_y=(Byp+j-MY)*JY; else pos_y = (Byp+j)*JY;
				for (int k=1; k<Mz+1; k++) { kk++; if (Bzp+k>MZ) pos_z=(Bzp+k-MZ); else pos_z = (Bzp+k);
					phi[pos_x+pos_y+pos_z]+=rho[pos_l+ii+jj+kk]*Inv_H_GNp;
				}
			}
		}
	}
}

#endif

//#ifdef CUDA
//		if (charges) {
//			phi_na =(double*)AllOnDev(MM);
//			phi_cl =(double*)AllOnDev(MM);
//			psi_0 =(double*)AllOnDev(MM);
//			psi_side =(double*)AllOnDev(MM);
//			psi=x+MM;
//			q =(double*)AllOnDev(MM);
//		}
#else
//		if (charges) {
//			psi=x+MM;
//			q = new double[MM];
//			psi_0 = new double[MM];
//			psi_side = new double[MM];
//		}
//#endif
*/
