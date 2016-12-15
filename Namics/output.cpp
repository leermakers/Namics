
class Output {
public:
	Output(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,vector<System*>,string,int,int);

~Output();

	string name;
	vector<Input*> In; 
	vector<Lattice*> Lat;
	vector<Segment*> Seg; 
	vector<Molecule*> Mol;
	vector<System*> Sys;
	int n_output; 
	int output_nr; 
	std::string template_;
	std::string type;
	bool write_bounds;
	bool append; 
		

	std::vector<string> OUT_name;
	std::vector<string> OUT_property;
	std::vector<string> OUT_value;

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckOutInput(void);
	void PutParameter(string); 
	string GetValue(string);
	bool Load(string); 
	void vtk(string, double *);
	void density();
	void printlist();


};
Output::Output(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,vector<System*> Sys_,string name_,int outnr,int N_out) {
	In=In_; Lat = Lat_; Seg=Seg_; Mol=Mol_; Sys=Sys_; name=name_; n_output=N_out; output_nr=outnr; 
	KEYS.push_back("template"); 
	KEYS.push_back("type");
	KEYS.push_back("write_bounds");
	KEYS.push_back("append"); 

}
Output::~Output() {
}
void Output::PutParameter(string new_param) {
	KEYS.push_back(new_param); 
}

bool Output::Load(string template_) {
return In[0]->LoadItems(template_, OUT_name, OUT_property, OUT_value);
}

bool Output::CheckOutInput() {
	return In[0]->CheckParameters("output",name,KEYS,PARAMETERS,VALUES);
}

string Output::GetValue(string parameter) {
	int length = PARAMETERS.size(); 
	int i=0;
	while (i<length) {
		if (PARAMETERS[i]==parameter) {return VALUES[i];}		
		i++;
	} 
	return " "; 
}

void Output::vtk(string filename, double *X){
	MX=Lat[0]->MX; MY=Lat[0]->MY;MZ=Lat[0]->MZ; 
	JX=Lat[0]->JX; JY=Lat[0]->JY;
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
	int length=In[0]->MolList.size();
	string fname;
	for (int i=0; i<length; i++) {
		fname = "Molecule_" +In[0]->MolList[i]+ "_Density.vtk";
		vtk(fname,Mol[i]->phitot);	
	}
}


void Output::printlist(){
	
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
