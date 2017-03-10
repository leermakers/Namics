#include "output.h"
#include "time.h"
Output::Output(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,vector<System*> Sys_,vector<Newton*> New_,vector<Engine*> Eng_,string name_,int outnr,int N_out) {
if (debug) cout <<"constructor in Output "<< endl; 
	In=In_; Lat = Lat_; Seg=Seg_; Mol=Mol_; Sys=Sys_; name=name_; n_output=N_out; output_nr=outnr;  New=New_; Eng=Eng_;
	KEYS.push_back("write_bounds");
	KEYS.push_back("append"); 
	input_error=false;
	//if (!CheckOutInput()) {input_error = true; cout << "Error found in ChcekOutInput in output module "<<endl;}
	//if (!Load()) {input_error=true;  cout <<"Error found in load output items in output module " << endl; }

}
Output::~Output() {
if (debug) cout <<"destructor in output " << endl; 
}
void Output::PutParameter(string new_param) {
if (debug) cout << "PutParameter in Output " << endl; 
	KEYS.push_back(new_param); 
}

bool Output::Load() {
if (debug) cout <<"Load in output " << endl; 
	bool success=true;
	int molnr=0;
	success= In[0]->LoadItems(name, OUT_key, OUT_name, OUT_prop);
	if (success) {	
		int length=OUT_key.size();

		for (int i=0; i<length; i++) {
			if (OUT_key[i]=="mol"){
				vector<string> sub;
				In[0]->split(OUT_prop[i],'*',sub);


				if (!(sub[0]==OUT_prop[i])){
					bool wildmon;
					if (sub[0] =="") wildmon=false; else wildmon=true; 

					int k=0; int mollength=In[0]->MolList.size();
					while (k<mollength) {
						if (In[0]->MolList[k]==OUT_name[i]) molnr=k; 
						k++;
					}
					if (wildmon) {
						int monlength=Mol[molnr]->MolMonList.size();
						for (int j=0; j<monlength; j++) {
							OUT_key.push_back(OUT_key[i]);
							OUT_name.push_back(OUT_name[i]);
							string s=sub[0];
							s=s.append(Seg[Mol[molnr]->MolMonList[j]]->name);
							s=s.append(sub[1]);
							OUT_prop.push_back(s);
						} 
					} else {
						int AlListlength=Mol[molnr]->MolAlList.size();
						for (int j=0; j<AlListlength; j++) {
							if (Mol[molnr]->Al[j]->value <1) { 
								OUT_key.push_back(OUT_key[i]);
								OUT_name.push_back(OUT_name[i]);
								string s="";
								s=s.append(Mol[molnr]->Al[j]->name);
								s=s.append(sub[1]);
								OUT_prop.push_back(s);
							}
						} 
					}
					OUT_key.erase(OUT_key.begin()+i);
					OUT_name.erase(OUT_name.begin()+i);
					OUT_prop.erase(OUT_prop.begin()+i);

				}

			}
		}

		if (name=="vtk" && OUT_key.size()>1) {
			cout << "vtk output can have only one member. Multiple entries were found namely: " << endl; 
			success=false;
			int length=OUT_key.size();
			for (int i=0; i<length; i++) cout << name << " : " << OUT_key[i] << " : " << OUT_name[i] << " : " << OUT_prop[i] << endl; 
		}

	}
	return success; 
}

bool Output::CheckInput(int start) {
if (debug) cout << "CheckInput in output " << endl; 
	bool success=true;
	success=In[0]->CheckParameters("output",name,start,KEYS,PARAMETERS,VALUES);
	if (success) {
		if (GetValue("append").size()>0) {
			if (name=="ana") append=true;
			append=In[0]->Get_bool(GetValue("append"),append); 
		} else {
			if (name=="ana") append=true;
			if (name=="kal") append=true;
			if (name=="pro") append=false;
			if (name=="vtk") append=false;
		} 
		if (GetValue("write_bounds").size()>0) {
			In[0]->Get_bool(GetValue("write_bounds"),write_bounds);
		} else write_bounds=false; 
		if (success) {if (!Load()) {cout <<"Error in Load() in output" << endl; success=false; }}
	} else cout <<"Error in CheckParameters in output" << endl; 
	return success; 
}

string Output::GetValue(string parameter) {
if (debug) cout << "GetValue in output " << endl; 
	int length = PARAMETERS.size(); 
	int i=0;
	while (i<length) {
		if (PARAMETERS[i]==parameter) {return VALUES[i];}		
		i++;
	} 
	return ""; 
}

Real* Output::GetPointer(string key, string name, string prop) {
if (debug) cout << "GetPointer in output " << endl; 
	int monlistlength=In[0]->MonList.size();
	int mollistlength=In[0]->MolList.size();
	//int aliaslistlength;
	int listlength; 
	int choice;
	int i,j; 
	if  (key=="sys") choice=1; if  (key=="mol") choice=2; if  (key=="mon") choice=3; if  (key=="engine") choice=4;

	switch(choice) {
		case 1:
			return Sys[0]->GetPointer(prop); 
			break;
		case 2:
			i=0;
			while (i<mollistlength){
				if (name==In[0]->MolList[i]) {
					listlength= Mol[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Mol[i]->strings[j]) return Mol[i]->GetPointer(Mol[i]->strings_value[j]); 
						j++;
					}
				}
				i++;
			}
			break;
		case 3:
			i=0;
			while (i<monlistlength){
				if (name==In[0]->MonList[i]) {
					listlength= Seg[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Seg[i]->strings[j]) return Seg[i]->GetPointer(Seg[i]->strings_value[j]); 
						j++;
					}
				}
				i++;
			}
			break;
		case 4:
			listlength= Eng[0]->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==Eng[0]->strings[j]) return Eng[0]->GetPointer(Eng[0]->strings_value[j]); 
				j++;
			}
			break;
		default:
			cout << "Program error: in Output, GetPointer reaches default...." << endl; 
	}
	return NULL; 
}
int Output::GetValue(string key, string name, string prop, int &int_result, Real &Real_result, string &string_result) {
if (debug) cout << "GetValue (long) in output " << endl; 
	int monlistlength=In[0]->MonList.size();
	int mollistlength=In[0]->MolList.size();
	int choice;
	int i; 
	if  (key=="sys") choice=1; if  (key=="mol") choice=2; if  (key=="mon") choice=3; if  (key=="newton") choice=4; if  (key=="lat") choice=5;if  (key=="engine") choice=6;

	switch(choice) {
		case 1:
			return Sys[0]->GetValue(prop,int_result,Real_result,string_result); 
			break;
		case 2:
			i=0;
			while (i<mollistlength){
				if (name==In[0]->MolList[i]) return Mol[i]->GetValue(prop,int_result,Real_result,string_result);
				i++;
			}
			break;
		case 3:
			i=0;
			while (i<monlistlength){
				if (name==In[0]->MonList[i]) return Seg[i]->GetValue(prop,int_result,Real_result,string_result);
				i++;
			}
			break;
		case 4:
			return New[0]->GetValue(prop,int_result,Real_result,string_result);
			break;
		case 5:
			return Lat[0]->GetValue(prop,int_result,Real_result,string_result);
			break;
		case 6:
			return Eng[0]->GetValue(prop,int_result,Real_result,string_result);
			break;
		default:
			cout << "Program error: in Output, GetValue reaches default...." << endl; 
	}
	return 0; 
}


void Output::WriteOutput(int subl) {
if (debug) cout << "WriteOutput in output " << endl; 
	int length;
	string s; 
	string filename;
	vector<string> sub;
	string infilename=In[0]->name;
	In[0]->split(infilename,'.',sub);
	string key;
	char numc[2];
        sprintf(numc,"%d",subl);
	filename=sub[0].append("_").append(numc).append(".").append(name); 
  	
	if (name=="pro") {
		vector<Real*> pointer;
		FILE *fp;
		fp=fopen(filename.c_str(),"w");
		int length=OUT_key.size();
		switch(Lat[0]->gradients) {
			case 1:
				fprintf(fp,"x \t"); 
				break;
			case 2: 
				fprintf(fp,"x \t y \t"); 
				break;
			case 3: 
				fprintf(fp,"x \t y \t z \t"); 
				break;
			default: 
				break;
		}
		//fprintf(fp,"x \t y \t z \t"); 
		for (int i=0; i<length; i++) {
			Real*  X = GetPointer(OUT_key[i],OUT_name[i],OUT_prop[i]);
			if (X!=NULL) {
				pointer.push_back(X); 
				key = OUT_key[i];
				string s=key.append(":").append(OUT_name[i]).append(":").append(OUT_prop[i]); 
				fprintf(fp,"%s \t",s.c_str());
			} else {cout << " Error for 'pro' output. It is only possible to ouput quantities known to be a 'profile'. That is why output quantity " + s + " is rejected. " << endl; }
		}
		fprintf(fp,"\n"); 
		Lat[0] -> PutProfiles(fp,pointer);

		fclose(fp);
	}

	if (name=="kal") {
		ifstream my_file(filename.c_str());
		FILE *fp;
		if (!(my_file||append)) {
			fp=fopen(filename.c_str(),"w");
			int length = OUT_key.size();
			for (int i=0; i<length; i++) {
				string s=OUT_key[i].append(":").append(OUT_name[i]).append(":").append(OUT_prop[i]); 
				fprintf(fp,"%s \t",s.c_str()); 
			}
			fprintf(fp,"\n"); append=true; 
		} else 	fp=fopen(filename.c_str(),"a"); 
		int length = OUT_key.size();
		for (int i=0; i<length; i++) {
			int int_result=0;
			int result_nr=0;
			Real Real_result=0;
			string string_result; 
			vector<string> sub;
			In[0] -> split(OUT_prop[i],'(',sub);
			result_nr= GetValue(OUT_key[i],OUT_name[i],sub[0],int_result,Real_result,string_result);
			if (result_nr==0) {fprintf(fp,"\t");}
			if (result_nr==1) {fprintf(fp,"%i\t",int_result);}
			if (result_nr==2) {fprintf(fp,"%f\t",Real_result);}
			if (result_nr==3) {
				if (sub[0]==OUT_prop[i]) {
					fprintf(fp,"%s\t",string_result.c_str());
				} else {
					Real* X=GetPointer(OUT_key[i],OUT_name[i],sub[0]);
					fprintf(fp,"%f\t",Lat[0]->GetValue(X,sub[1])); 
				}
			}
		}
		fprintf(fp,"\n"); 
		fclose(fp); 
	}

	if (name=="vtk") {
		Real*  X = GetPointer(OUT_key[0],OUT_name[0],OUT_prop[0]);	
		string s=OUT_key[0].append(" : ").append(OUT_name[0]).append(" : ").append(OUT_prop[0]);
		if (!(X==NULL)) 
		Lat[0]->vtk(filename,X,s); else {cout << "vtk file was not generated because 'profile' was not found" << endl;}
	}


	if (name=="ana") {
		time_t now;
		time(&now); 
		FILE *fp;
		if (append) fp=fopen(filename.c_str(),"a"); else fp=fopen(filename.c_str(),"w");
		fprintf(fp,"version: %s %s",version.c_str(),ctime(&now));
//System parameters		
		s="sys : " + Sys[0]->name + " :";  
		length = Sys[0]->ints.size();
		for (int i=0; i<length; i++) 
			fprintf(fp,"%s %s : %i \n",s.c_str(),Sys[0]->ints[i].c_str(),Sys[0]->ints_value[i]);
		length = Sys[0]->Reals.size();
		for (int i=0; i<length; i++) fprintf(fp,"%s %s : %e \n",s.c_str(),Sys[0]->Reals[i].c_str(),Sys[0]->Reals_value[i]);
		length = Lat[0]->bools.size(); 
		for (int i=0; i<length; i++) {
			if (Sys[0]->bools_value[i]) fprintf(fp,"%s %s : %s \n",s.c_str(),Sys[0]->bools[i].c_str(),"true");
			else fprintf(fp,"%s %s : %s \n",s.c_str(),Sys[0]->bools[i].c_str(),"false");
		}
		length = Sys[0]->strings.size();
		for (int i=0; i<length; i++) fprintf(fp,"%s %s : %s \n",s.c_str(),Sys[0]->strings[i].c_str(),Sys[0]->strings_value[i].c_str());

//Lattice parameters
		s="lat : " + Lat[0]->name + " :";  
		length = Lat[0]->ints.size();
		for (int i=0; i<length; i++) 
			fprintf(fp,"%s %s : %i \n",s.c_str(),Lat[0]->ints[i].c_str(),Lat[0]->ints_value[i]);
		length = Lat[0]->Reals.size();
		for (int i=0; i<length; i++) fprintf(fp,"%s %s : %e \n",s.c_str(),Lat[0]->Reals[i].c_str(),Lat[0]->Reals_value[i]);
		length = Lat[0]->bools.size(); 
		for (int i=0; i<length; i++) {
			if (Lat[0]->bools_value[i]) fprintf(fp,"%s %s : %s \n",s.c_str(),Lat[0]->bools[i].c_str(),"true");
			else fprintf(fp,"%s %s : %s \n",s.c_str(),Lat[0]->bools[i].c_str(),"false");
		}
		length = Lat[0]->strings.size();
		for (int i=0; i<length; i++) fprintf(fp,"%s %s : %s \n",s.c_str(),Lat[0]->strings[i].c_str(),Lat[0]->strings_value[i].c_str());

//Newton parameters
		s="newton : " + New[0]->name + " :";  
		length = New[0]->ints.size();
		for (int i=0; i<length; i++) 
			fprintf(fp,"%s %s : %i \n",s.c_str(),New[0]->ints[i].c_str(),New[0]->ints_value[i]);
		length = New[0]->Reals.size();
		for (int i=0; i<length; i++) fprintf(fp,"%s %s : %e \n",s.c_str(),New[0]->Reals[i].c_str(),New[0]->Reals_value[i]);
		length = New[0]->bools.size(); 
		for (int i=0; i<length; i++) {
			if (New[0]->bools_value[i]) fprintf(fp,"%s %s : %s \n",s.c_str(),New[0]->bools[i].c_str(),"true");
			else fprintf(fp,"%s %s : %s \n",s.c_str(),New[0]->bools[i].c_str(),"false");
		}
		length = New[0]->strings.size();
		for (int i=0; i<length; i++) fprintf(fp,"%s %s : %s \n",s.c_str(),New[0]->strings[i].c_str(),New[0]->strings_value[i].c_str());

//Engine parameters
		s="engine : " + Eng[0]->name + " :";  
		length = Eng[0]->ints.size();
		for (int i=0; i<length; i++) 
			fprintf(fp,"%s %s : %i \n",s.c_str(),Eng[0]->ints[i].c_str(),Eng[0]->ints_value[i]);
		length = Eng[0]->Reals.size();
		for (int i=0; i<length; i++) fprintf(fp,"%s %s : %e \n",s.c_str(),Eng[0]->Reals[i].c_str(),Eng[0]->Reals_value[i]);
		length = Eng[0]->bools.size(); 
		for (int i=0; i<length; i++) {
			if (Eng[0]->bools_value[i]) fprintf(fp,"%s %s : %s \n",s.c_str(),Eng[0]->bools[i].c_str(),"true");
			else fprintf(fp,"%s %s : %s \n",s.c_str(),Eng[0]->bools[i].c_str(),"false");
		}
		length = Eng[0]->strings.size();
		for (int i=0; i<length; i++) fprintf(fp,"%s %s : %s \n",s.c_str(),Eng[0]->strings[i].c_str(),Eng[0]->strings_value[i].c_str());

//segment parameters
		int length_A=In[0]->MonList.size();
		for (int j=0; j<length_A; j++) {
			s="seg : " + Seg[j]->name + " :";  
			length = Seg[j]->ints.size();
			for (int i=0; i<length; i++) 
				fprintf(fp,"%s %s : %i \n",s.c_str(),Seg[j]->ints[i].c_str(),Seg[j]->ints_value[i]);
			length = Seg[j]->Reals.size();
			for (int i=0; i<length; i++) fprintf(fp,"%s %s : %e \n",s.c_str(),Seg[j]->Reals[i].c_str(),Seg[j]->Reals_value[i]);
			length = Seg[j]->bools.size(); 
			for (int i=0; i<length; i++) {
				if (Seg[j]->bools_value[i]) fprintf(fp,"%s %s : %s \n",s.c_str(),Eng[j]->bools[i].c_str(),"true");
				else fprintf(fp,"%s %s : %s \n",s.c_str(),Eng[j]->bools[i].c_str(),"false");
			}
			length = Seg[j]->strings.size();
			for (int i=0; i<length; i++) fprintf(fp,"%s %s : %s \n",s.c_str(),Seg[j]->strings[i].c_str(),Seg[j]->strings_value[i].c_str());
		}
//molecule parameters
		length_A=In[0]->MolList.size();
		for (int j=0; j<length_A; j++) {
			s="mol : " + Mol[j]->name + " :";  
			length = Mol[j]->ints.size();
			for (int i=0; i<length; i++) 
				fprintf(fp,"%s %s : %i \n",s.c_str(),Mol[j]->ints[i].c_str(),Mol[j]->ints_value[i]);
			length = Mol[j]->Reals.size();
			for (int i=0; i<length; i++) fprintf(fp,"%s %s : %e \n",s.c_str(),Mol[j]->Reals[i].c_str(),Mol[j]->Reals_value[i]);
			length = Mol[j]->bools.size(); 
			for (int i=0; i<length; i++) {
				if (Mol[j]->bools_value[i]) fprintf(fp,"%s %s : %s \n",s.c_str(),Mol[j]->bools[i].c_str(),"true");
				else fprintf(fp,"%s %s : %s \n",s.c_str(),Mol[j]->bools[i].c_str(),"false");
			}
			length = Mol[j]->strings.size();
			for (int i=0; i<length; i++) fprintf(fp,"%s %s : %s \n",s.c_str(),Mol[j]->strings[i].c_str(),Mol[j]->strings_value[i].c_str());
		}

		fprintf(fp,"%s \n","system delimiter"); 
		fclose(fp); 
	}
}

void Output::vtk(string filename, Real *X){
if (debug) cout << "vtk in output " << endl; 
	int MX=Lat[0]->MX; 
	int MY=Lat[0]->MY;
	int MZ=Lat[0]->MZ; 
	int JX=Lat[0]->JX; 
	int JY=Lat[0]->JY;
	FILE *fp;
	fp = fopen(filename.c_str(),"w+");
	fprintf(fp, "# vtk DataFile Version 7.0 \nvtk output \nASCII \nDATASET STRUCTURED_POINTS \nDIMENSIONS %i %i %i\n",MX,MY,MZ);
	fprintf(fp, "SPACING 1 1 1 \nORIGIN 0 0 0 \nPOINT_DATA %i\n", MX*MY*MZ);
	fprintf(fp, "SCALARS Box_profile float\nLOOKUP_TABLE default \n");
	
	for(int i=1; i<MX+1; i++) for(int j=1; j<MY+1; j++) for(int k=1; k<MZ+1; k++)
		fprintf(fp,"%lf \n",X[i*JX+j*JY+k]);
	
	fclose(fp);
}


void Output::density(){
if (debug) cout << "density in output " << endl; 
	int length=In[0]->MolList.size();
	string fname;
	for (int i=0; i<length; i++) {
		fname = "Molecule_" +In[0]->MolList[i]+ "_Density.vtk";
		vtk(fname,Mol[i]->phitot);	
	}
}


void Output::printlist(){
if (debug) cout << "printlist in output " << endl; 
	
	ofstream writefile;
	writefile.open ("listdetails.cpp");
	writefile << "---------------------------------------------------------------" << endl;
	
	
	writefile << "-------------------- Input list details -----------------------" << endl;
	int length = In[0]->MonList.size();
	writefile << " " << endl;
	writefile << "MonList" << endl;
	writefile << " " << endl;
	for (int i=0; i<length; i++) writefile << "i: " <<  i << "\t " << In[0]->MonList[i] << endl;
	
	length = In[0]->MolList.size();
	writefile << " " << endl;
	writefile << "MolList" << endl;
	writefile << " " << endl;
	for (int i=0; i<length; i++) writefile << "i: " <<  i << "\t " << In[0]->MolList[i] << endl;
	writefile << "---------------------------------------------------------------" << endl;
	writefile << " " << endl;
	
	writefile << "-------------------- Molecule list details -----------------------" << endl;
	int n_mol = In[0]->MolList.size();
	for (int j=0; j<n_mol; j++){
		writefile << " " << endl;
		writefile << "Molecule: " << In[0]->MolList[j]<< endl;
		writefile << " " << endl;
		
		length = Mol[j]->MolMonList.size();
		writefile << "MolMonList" << endl;
		for (int i=0; i<length; i++) writefile << "i: " <<  i << "\t " << Mol[j]->MolMonList[i] << endl;
		
		length = Mol[j]->mon_nr.size();
		writefile << "mon_nr" << endl;
		for (int i=0; i<length; i++) writefile << "i: " <<  i << "\t " << Mol[j]->mon_nr[i] << endl;
		
		length = Mol[j]->n_mon.size();
		writefile << "n_mon" << endl;
		for (int i=0; i<length; i++) writefile << "i: " <<  i << "\t " << Mol[j]->n_mon[i] << endl;
		
		length = Mol[j]->molmon_nr.size();
		writefile << "molmon_nr" << endl;
		for (int i=0; i<length; i++) writefile << "i: " <<  i << "\t " << Mol[j]->molmon_nr[i] << endl;
	}
	writefile << "---------------------------------------------------------------" << endl;
	writefile << " " << endl;

	writefile << "-------------------- System list details -----------------------" << endl;
	length = Sys[0]->SysMonList.size();
	writefile << " " << endl;
	writefile << "SysMonList" << endl;
	writefile << " " << endl;
	for (int i=0; i<length; i++) writefile << "i: " <<  i << "\t " << Sys[0]->SysMonList[i] << endl;
	
	length = Sys[0]->SysTagList.size();
	writefile << " " << endl;
	writefile << "SysTagList" << endl;
	writefile << " " << endl;
	for (int i=0; i<length; i++) writefile << "i: " <<  i << "\t " << Sys[0]->SysTagList[i] << endl;
	writefile << "---------------------------------------------------------------" << endl;
	writefile << " " << endl;

	length = Sys[0]->FrozenList.size();
	writefile << " " << endl;
	writefile << "FrozenList" << endl;
	writefile << " " << endl;
	for (int i=0; i<length; i++) writefile << "i: " <<  i << "\t " << Sys[0]->FrozenList[i] << endl;
	writefile << "---------------------------------------------------------------" << endl;
	writefile << " " << endl;
	writefile.close();
	
}


/*
#ifdef CUDA
	TransferDataToHost(H_u,u,MM);
#endif
	vtk_output(fname+"_pot.vtk",H_u);
#ifdef CUDA
	TransferDataToHost(H_phi,phi,4*MM);
#endif

	Cp(PHI,H_phi+2*MM,MM); Add(PHI,phi+3*MM,MM); 
	vtk_output(fname+"_phi.vtk",PHI);

	ofstream varfile;
	varfile.open ("surf_profile.dat");
	for (int z=1; z<MZ; z++){	
		varfile << z <<"\t" << PHI[10*JX+40*JY+ z] << endl;
		}
	varfile.close();
	
	ofstream varfile2;
	varfile2.open ("surf_profile_part_1.dat");
	for (int z=1; z<MZ; z++){	
		varfile2 << z <<"\t" << PHI[15*JX+35*JY+ z] << endl;
		}
	varfile2.close();
	
	
	
	
	Cp(PHI,H_phi,MM);
	vtk_output(fname+"A_phi.vtk",PHI);
	
	ofstream varfile3;
	varfile3.open ("mol_profile.dat");
	for (int z=1; z<MZ; z++){	
		varfile3 << z <<"\t" << PHI[10*JX+40*JY+ z] << endl;
		}
	varfile3.close();
	
	ofstream varfile4;
	varfile4.open ("mol_profile_part_1.dat");
	for (int z=1; z<MZ; z++){	
		varfile4 << z <<"\t" << PHI[15*JX+35*JY+ z] << endl;
		}
	varfile4.close();

	
	Cp(PHI,H_phi+MM,MM);
	vtk_output(fname+"B_phi.vtk",PHI);
	filename=fname+".out";
	out_file.open(filename.c_str());
	out_file << "Free_energy : " << Free_energy << endl;
	out_file.close();
*/
