#include "system.h"
System::System(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,string name_) {
	Seg=Seg_; Mol=Mol_; Lat=Lat_; In=In_; name=name_;
if (debug) cout << "Constructor for system " << endl;
	KEYS.push_back("calculation_type");
	KEYS.push_back("generate_guess");
	KEYS.push_back("initial_guess");
	KEYS.push_back("guess_inputfile");
	KEYS.push_back("final_guess");
	KEYS.push_back("guess_outputfile");
	KEYS.push_back("GPU");
	int length = In[0]->MonList.size();
	for (int i=0; i<length; i++) KEYS.push_back("guess-"+In[0]->MonList[i]);
	AllocateMemory();
	PrepareForCalculations();
}
System::~System() {
if (debug) cout << "Destructor for system " << endl;
	free(H_GrandPotentialDensity);
	free(H_FreeEnergyDensity);
	free(H_alpha);
	if (charged) {
		free(H_q);
		free(H_psi);
	}
#ifdef CUDA
	cudaFree(phitot);
	cudaFree(alpha);
	cudaFree(GrandPotentialDensity);
	cudaFree(FreeEnergyDensity);
	cudaFree(TEMP);
	cudaFree(KSAM);
	if (charged) {
		cudaFree(psi);
		cudaFree(eps);
		cudaFree(q);
		cudaFree(EE);
		cudeFree(psiMask);
	}
#else
	free(phitot);
	free(TEMP);
	free(KSAM);
	free(CHI);
	if (charged) {
		free(EE);
		free(psiMask);
		free(eps);
	}
#endif
}

void System::AllocateMemory() {
if (debug) cout << "AllocateMemory in system " << endl;
	int M=Lat[0]->M;
	//H_GN_A = new Real[n_box];
	//H_GN_B = new Real[n_box];

	H_GrandPotentialDensity = (Real*) malloc(M*sizeof(Real));
	H_FreeEnergyDensity = (Real*) malloc(M*sizeof(Real));
	H_alpha = (Real*) malloc(M*sizeof(Real));
	if (charged) {
		H_q=(Real*) malloc(M*sizeof(Real));
		H_psi=(Real*) malloc(M*sizeof(Real));

	}
#ifdef CUDA
	phitot = (Real*)AllOnDev(M);
	alpha=(Real*)AllOnDev(M);
	GrandPotentialDensity =(Real*)AllOnDev(M);
	FreeEnergyDensity=(Real*)AllOnDev(M);
	TEMP =(Real*)AllOnDev(M);
	KSAM =(int*)AllOnDev(M);
	if (charged) {
		psi=(Real*)AllOnDev(M);
		q=(Real*)AllOnDev(M);
		eps=(Real*)AllOnDev(M);
		EE=(Real*)AllOnDev(M);
		psiMask=(int*)AllOnDev(M);
	}
#else
	phitot = (Real*) malloc(M*sizeof(Real));
	alpha= H_alpha;
	if (charged) {
		psi=H_psi; q=H_q;
		eps=(Real*) malloc(M*sizeof(Real));
		EE=(Real*) malloc(M*sizeof(Real));
		psiMask = (int*) malloc(M*sizeof(int));
	}
	KSAM = (int*) malloc(M*sizeof(int));
	FreeEnergyDensity=H_FreeEnergyDensity;
	GrandPotentialDensity = H_GrandPotentialDensity;
	TEMP = (Real*) malloc(M*sizeof(Real));
#endif
	Zero(KSAM,M);
	if (charged) {
		Zero(psi,M);
		Zero(EE,M);
	}
	n_mol = In[0]->MolList.size();
	int i=0;
	Lat[0]->AllocateMemory();
	int n_mon=In[0]->MonList.size();
	while (i<n_mon) {Seg[i]->AllocateMemory(); i++;}
	i=0;
	while (i<n_mol) { Mol[i]->AllocateMemory(); i++;}
}

bool System::PrepareForCalculations() {
if (debug) cout <<"PrepareForCalculations in System " << endl;
	int M=Lat[0]->M;
	bool success=true;

	FrozenList.clear();
	int length = In[0]->MonList.size();
	for (int i=0; i<length; i++) {if (Seg[i]->freedom == "frozen") FrozenList.push_back(i); }
	if (FrozenList.size()+SysMonList.size()+SysTagList.size()+SysClampList.size() !=In[0]->MonList.size()) {cout <<" There are un-used monomers in system. Remove them before starting" << endl; success=false;}

	Zero(KSAM,M);

	length=FrozenList.size();
	for (int i=0; i < length ; ++i) {
		int*MASK=Seg[FrozenList[i]]->MASK;
		Add(KSAM,MASK,M);
	}

	length=SysTagList.size();
	for (int i=0; i < length ; ++i) {
		int* MASK=Seg[SysTagList[i]]->MASK;
		Add(KSAM,MASK,M);
	}

	length=SysClampList.size();
	for (int i=0; i < length ; ++i) {
		int* MASK=Seg[SysClampList[i]]->MASK;
		Add(KSAM,MASK,M);
	}

	Invert(KSAM,KSAM,M);

	volume = 0;
	for (int i = 0; i < M ; ++i) {
		volume += KSAM[i];
	}

	n_mol = In[0]->MolList.size();
	success=Lat[0]->PrepareForCalculations();
	int n_mon=In[0]->MonList.size();

	for (int i=0; i<n_mon; i++) {success=Seg[i]->PrepareForCalculations(KSAM);}
	for (int i=0; i<n_mol; i++) {success=Mol[i]->PrepareForCalculations(KSAM);}
	if (charged) {
		length=FrozenList.size();
		int i=0; Zero(psiMask,M); fixedPsi0=false;
		while (i<length) {
			if (Seg[FrozenList[i]]->fixedPsi0) {
				fixedPsi0=true;
				Add(psiMask,Seg[FrozenList[i]]->MASK,M);
			}
			i++;
		}
	}
	return success;
}

string System:: GetMonName(int mon_number_I){
if (debug) cout << "GetMonName for system " << endl;
	return Seg[mon_number_I]->name;
}

bool System::CheckInput(int start) {
if (debug) cout << "CheckInput for system " << endl;
	bool success=true;
	bool solvent_found=false; tag_segment=-1; solvent=-1;  //value -1 means no solvent defined. tag_segment=-1;
	Real phibulktot=0;
	success= In[0]->CheckParameters("sys",name,start,KEYS,PARAMETERS,VALUES);
	if (success) {
		success=CheckChi_values(In[0]->MonList.size());
		GPU=In[0]->Get_bool(GetValue("GPU"),false);
		if (GPU) {if (!cuda) cout << "You should compile the program using the CUDA=1 flag: GPU calculations are impossible; proceed with CPU computations..." << endl; GPU=false; }
		if (Lat[0]->gradients<3) {
			if (GPU) cout <<"GPU support is (for the time being) only available for three-gradient calculations " << endl;
		}
		if (cuda) {if (!GPU) cout <<" program expect that you are going to use the GPU, but the input is not in line with this (either gradients < 3, or GPU != 'true' : compile without CUDA=1 flag." << endl; success=false;}
		int i=0;
		int length = In[0]->MolList.size();
		while (i<length) {
			int j=0;
			int LENGTH=Mol[i]->MolMonList.size();
			while (j<LENGTH) {
				SysMolMonList.push_back(Mol[i]->MolMonList[j]);
				if (!In[0]->InSet(SysMonList,Mol[i]->MolMonList[j])) {
					if (Seg[Mol[i]->MolMonList[j]]->freedom!="tagged" && Seg[Mol[i]->MolMonList[j]]->freedom!="clamp") SysMonList.push_back(Mol[i]->MolMonList[j]);
				}
				if (Seg[Mol[i]->MolMonList[j]]->freedom=="tagged"){
					if (In[0]->InSet(SysTagList,Mol[i]->MolMonList[j])) {
						//cout <<"You can not use the 'tag monomer' " + GetMonName(Mol[i]->MolMonList[j]) + " in more than one molecule." << endl; success=false;
					} else	SysTagList.push_back(Mol[i]->MolMonList[j]);
				}
				if (Seg[Mol[i]->MolMonList[j]]->freedom=="clamp"){
					if (In[0]->InSet(SysClampList,Mol[i]->MolMonList[j])) {
						//cout <<"You can not use the 'clamp monomer' " + GetMonName(Mol[i]->MolMonList[j]) + " in more than one molecule." << endl; success=false;
					} else	SysClampList.push_back(Mol[i]->MolMonList[j]);
				}
				j++;
			}
			if (Mol[i]->freedom=="free") phibulktot +=Mol[i]->phibulk;
			if (Mol[i]->freedom=="solvent") {
				solvent_found=true; solvent=i;
			}
			if (Mol[i]->IsTagged()) {
				tag_segment = Mol[i]->tag_segment;
				Mol[i]->n=1.0*Seg[tag_segment]->n_pos;
				Mol[i]->theta= Mol[i]->n*Mol[i]->chainlength;
			}
			i++;
		}
		if (!solvent_found ||( phibulktot >0.99999999 && phibulktot < 1.0000000001)) {
			cout << "In system '" + name + "' the 'solvent' was not found, while the volume fractions of the bulk do not add up to unity. " << endl;
			success=false;
		}
		charged=false;
		if (IsCharged()) {
			charged=true;
			bool neutralizer_needed=false;
			bool neutralizer_found=false;
			Real phibulk=0;
			int length=In[0]->MolList.size();
			for (int i=0; i<length; i++) {
				if (Mol[i]->freedom=="neutralizer") {neutralizer_found=true; neutralizer=i;}
				if (Mol[i]->freedom=="restricted") neutralizer_needed=true;
				if (Mol[i]->freedom=="free") phibulk+=Mol[i]->Charge()*Mol[i]->phibulk;
				if (Mol[i]->freedom=="solvent") {
					if (Mol[i]->IsCharged()) {
						success=false;
						cout <<"The solvent is charged and this gives problems with neutralization: problem terminated." << endl;
					}
				}
			}
			if (phibulk!=0 || neutralizer_needed) {
				if (!neutralizer_found) {
					cout <<"Neutralizer needed because some molecules have freedom 'restricted' or charge density of 'free' molecules in the bulk is not equal to zero " << endl;
					success=false;
				}
			}
		}

		vector<string> options;
		options.push_back("micro_emulsion"); options.push_back("micro_phasesegregation");
		CalculationType="";
		if (GetValue("calculation_type").size()>0) {
			if (!In[0]->Get_string(GetValue("calculation_type"),CalculationType,options," Info about calculation_type rejected") ) {success=false;};
			if (CalculationType=="micro_emulsion") {
				if (In[0]->MolList.size()<3) {cout << "In 'calculation_type : micro-emulsion', we expect at least three types of molecules " << endl; success=false;  }
				int length=In[0]->MolList.size();
				int i=0;
				int number_of_solvents=0;
				int number_of_surfactants=0;
				while (i<length) {
					if (Mol[i]->IsTagged()) {number_of_surfactants++;} else {number_of_solvents++; }
					i++;
				}
				if (number_of_solvents <2 || number_of_surfactants <1) {
					cout << "In 'calculation_type : micro-emulsion', we expect at least two solvents and one surfactant which is tagged. " << endl;
					success=false;
				}
			}
			if (CalculationType=="micro_phasesegregation") {
				int length=In[0]->MolList.size();
				bool found=false;
				int i=0;
				while (i<length && !found) {
					if (Mol[i]->MolMonList.size()>1) found=true;
				}
				if (!found) {cout << "In 'calculation_type : micro_phasesegregation', we expect at least copolymer in the system " << endl; success=false;  }
			}
		}
		options.clear();
		options.push_back("lamellae"); options.push_back("Im3m");
		GuessType=""; MonA=0; MonB=0;
		if (GetValue("generate_guess").size()>0) {
			if (!In[0]->Get_string(GetValue("generate_guess"),GuessType,options," Info about 'generate_guess' rejected") ) {success=false;};
			if (GuessType=="lamellae" || GuessType=="Im3m") {
				int length = In[0]->MonList.size();
				int n_guess=0;
				for (int i=0; i<length; i++) {
					string s="guess-"+In[0]->MonList[i];
					if (GetValue(s).size()>0) {
						Seg[i]->guess_u=In[0]->Get_Real(GetValue(s),0);
						if (Seg[i]->guess_u <-2 || Seg[i]->guess_u >2) {
							cout << "Suggested 'guess value' for 'u' for mon '" + In[0]->MonList[i] + "' out of range: current range is -2, .., 2. Value ignored. " << endl;
							Seg[i]->guess_u = 0;
						} else {
							n_guess++;
							if (i==0 && n_guess==1)  MonA=i;
							if (i==1 && n_guess==2)  MonB=i;
						}
					}
				}
				if (n_guess < 2) {
					cout <<"generate_guess needs two valid non-zero (real) values for 'guess_xxx' quantities. Here xxx is a monomomer name; 'generate_guess : " + GuessType + "' is most likely ineffective.  " << endl;
					cout <<"The idea is that the first 'mon' quantity will be the 'oil' and the second 'mon' quantity the 'water' component. Each filling one half-space. Adjust input accordingly. " << endl;
				}
			}

		}
		initial_guess="previous_result";
		if (GetValue("initial_guess").size()>0) {
			options.clear();
			options.push_back("previous_result"); options.push_back("file");
			if (!In[0]->Get_string(GetValue("initial_guess"),initial_guess,options," Info about 'initial_guess' rejected; default: 'previous_result' used.") ) {initial_guess="previous_result";}
			if (initial_guess=="file") {
				if (GetValue("guess_inputfile").size()>0) {
					guess_inputfile=GetValue("guess_inputfile");
				} else {
					success=false; cout <<" When 'initial_guess' is set to 'file', you need to supply 'guess_inputfile', but this entry is missing. Problem terminated " << endl;
				}
			}
		}
		final_guess="next_problem";
		if (GetValue("final_guess").size()>0) {
			options.clear();
			options.push_back("next_problem"); options.push_back("file");
			if (!In[0]->Get_string(GetValue("final_guess"),final_guess,options," Info about 'final_guess' rejected; default: 'next_problem' used.") ) {final_guess="next_problem";}
			if (final_guess=="file") {
				if (GetValue("guess_outputfile").size()>0) {
					guess_outputfile=GetValue("guess_outputfile");
				} else {
					guess_outputfile="";
					cout <<"Filename not found for 'output_guess'. Default with inputfilename and extention '.outiv' is used. " << endl;
				}
			}
		}
	}
	return success;
}

bool System::PutVarInfo(string Var_type_, string Var_target_, Real Var_target_value_){
if (debug) cout << "System::PutVarInfo " << endl;
	bool success=true;
	Var_target=-1;
	if (Var_type_ !="target") success=false;
	if (Var_target_=="free_energy") Var_target=0;
	if (Var_target_=="grand_potential") {Var_target=1; }
	if (Var_target<0 || Var_target>1) {success=false; cout <<"Var target " + Var_target_ + " rejected in PutVarInfo in System " << endl; }
	Var_target_value=Var_target_value_;
	if (Var_target_value < -1e4 || Var_target_value > 1e4) success=false;
	return success;
}

Real System::GetError() {
if (debug) cout << "System::GetError " << endl;
	Real Error=0;
	switch (Var_target) {
		case 0:
			Error = FreeEnergy - Var_target_value;
			break;
		case 1:
			Error = -1.0*(GrandPotential- Var_target_value);
			break;
		default:
			cout <<"Program error in GetVarError" <<endl;
			break;
	}
	return Error;
}

bool System::IsCharged() {
if (debug) cout << "System::IsCharged " << endl;
	bool success=false;
	int length = In[0]->MolList.size();
	for (int i=0; i<length; i++) {if (Mol[i]->IsCharged()) success=true; }
	return success;
}

void System::PutParameter(string new_param) {
if (debug) cout << "PutParameter for system " << endl;
	KEYS.push_back(new_param);
}

string System::GetValue(string parameter){
if (debug) cout << "GetValue " + parameter + " for system " << endl;
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

void System::push(string s, Real X) {
if (debug) cout << "push (Real) for system " << endl;
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void System::push(string s, int X) {
if (debug) cout << "push (int) for system " << endl;
	ints.push_back(s);
	ints_value.push_back(X);
}
void System::push(string s, bool X) {
if (debug) cout << "push (bool) for system " << endl;
	bools.push_back(s);
	bools_value.push_back(X);
}
void System::push(string s, string X) {
if (debug) cout << "push (string) for system " << endl;
	strings.push_back(s);
	strings_value.push_back(X);
}
void System::PushOutput() {
if (debug) cout << "PushOutput for system " << endl;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	push("e",e);
	push("k_B",k_B);
	push("eps0",eps0);
	push("temperature",T);
	push("free_energy",FreeEnergy);
	push("grand_potential",GrandPotential);
	push("calculation_type",CalculationType);
	push("guess_type",GuessType);
	push("cuda",cuda);
	push("solvent",Mol[solvent]->name);
	int n_mon=In[0]->MonList.size();
	for(int i=0; i<n_mon; i++) for(int j=0; j<n_mon; j++) if (i!=j) Seg[i]->push("chi-" + Seg[j]->name,CHI[i*n_mon+j]);
	string s="profile;0"; push("alpha",s);
	s="profile;1"; push("GrandPotentialDensity",s);
	s="profile;2"; push("FreeEnergyDensity",s);
	if (charged) {
		s="profile;3"; push("psi",s);
		s="profile;4"; push("q",s);
		s="profile;5"; push("eps",s);
	}
#ifdef CUDA
	TransferDataToHost(H_alpha, alpha, M);
	TransferDataToHost(H_GrandPotentialDensity, GrandPotentialDensity, M);
	TransferDataToHost(H_FreeEnergyDensity, FreeEnergyDensity, M);
#endif
}

Real* System::GetPointer(string s,int &SIZE){
if (debug) cout << "GetPointer for system " << endl;
	vector<string>sub;
	SIZE=Lat[0]->M;
	In[0]->split(s,';',sub);
	if (sub[1]=="0") return H_alpha;
	if (sub[1]=="1") return H_GrandPotentialDensity;
	if (sub[1]=="2") return H_FreeEnergyDensity;
	if (sub[1]=="3") return psi;
	if (sub[1]=="4") return q;
	if (sub[1]=="5") return eps;
	return NULL;
}
int* System::GetPointerInt(string s,int &SIZE){
if (debug) cout << "GetPointerInt for system " << endl;
	vector<string>sub;
	In[0]->split(s,';',sub);
	if (sub[0]=="array") { //set SIZE and return pointer of int array
	}
	return NULL;
}

int System::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
if (debug) cout << "GetValue (long) for system " << endl;
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

bool System::CheckChi_values(int n_seg){
if (debug) cout << "CheckChi_values for system " << endl;
	bool success=true;
	CHI = (Real*) malloc(n_seg*n_seg*sizeof(Real));
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) {
		CHI[i*n_seg+k]=In[0]->Get_Real(Seg[i]->GetValue("chi-"+Seg[k]->name),123);
	}
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) if (CHI[i*n_seg+k] == 123) CHI[i*n_seg+k] = CHI[k*n_seg+i];
	for (int i=0; i<n_seg; i++) for (int k=0; k<n_seg; k++) if (CHI[i*n_seg+k] == 123) CHI[i*n_seg+k] = 0;
	for (int i=0; i<n_seg; i++) for (int k=i+1; k<n_seg; k++) if (CHI[i*n_seg+k] != CHI[k*n_seg+i]) {
		cout <<"CHI-value symmetry violated: chi("<<Seg[i]->name<<","<<Seg[k]->name<<") is not equal to chi("<<Seg[k]->name<<","<<Seg[i]->name<<")"<<endl; success=false;
	}
	for (int i=0; i<n_seg; i++) {if (CHI[i*n_seg+i]!=0) {cout <<"CHI-values for 'same'-segments (e.g. CHI(x,x) ) should be zero. " << endl; success = false;} }
 	return success;
}

void System::DoElectrostatics(Real* g, Real* x) {
	int M=Lat[0]->M;
	int n_seg=In[0]->MonList.size();
	Zero(q,M);
	Zero(eps,M);
	for (int i=0; i<n_seg; i++) {
		YplusisCtimesX(q,Seg[i]->phi,Seg[i]->valence,M);
		YplusisCtimesX(eps,Seg[i]->phi,Seg[i]->epsilon,M);
	}
	Div(q,phitot,M);
	Div(eps,phitot,M);
	Lat[0]->set_bounds(eps);
	Cp(psi,x,M); Cp(g,psi,M);
	Lat[0]->set_bounds(psi);
	if (fixedPsi0) {
		int length=FrozenList.size();
		for (int i=0; i<length; i++) {
			Seg[FrozenList[i]]->UpdateValence(g,psi,q,eps);
		}
	}
}

bool System::ComputePhis(){
if(debug) cout <<"ComputePhis in system" << endl;
	int M= Lat[0]->M;
	bool success=true;

	Real totphibulk=0;
	Real norm=0;
	Zero(phitot,M);
	int length=FrozenList.size();
	for (int i=0; i<length; i++) {
		Real *phi_frozen=Seg[FrozenList[i]]->phi;
		Add(phitot,phi_frozen,M);
	}

	for (int i=0; i<n_mol; i++) {
		success=Mol[i]->ComputePhi();
	}
	for (int i=0; i<n_mol; i++) {
		norm=0;
		int length=Mol[i]->MolMonList.size();
		if (Mol[i]->freedom=="free") {
			norm=Mol[i]->phibulk/Mol[i]->chainlength;
			totphibulk +=Mol[i]->phibulk;
			Mol[i]->n=norm*Mol[i]->GN;
		}

		if (Mol[i]->freedom=="restricted" ) {
			if (Mol[i]->GN>0) {
			norm = Mol[i]->n/Mol[i]->GN;
			if (Mol[i]->IsPinned()) Mol[i]->phibulk = 0; else Mol[i]->phibulk=Mol[i]->chainlength*norm;  totphibulk +=Mol[i]->phibulk;
			} else { norm = 0; cout <<"GN for molecule " << i << " is not larger than zero..." << endl; }
		}
		if (Mol[i]->IsTagged() || Mol[i]->IsPinned()) {
			if (Mol[i]->GN>0) norm=Mol[i]->n/Mol[i]->GN; else {norm=0; cout <<"GN for molecule " << i << " is not larger than zero..." << endl; }
		}
		if (Mol[i]->IsClamped() ) {
			norm=1;
			Mol[i]->phibulk=0;
		}
		int k=0;
		Mol[i]->norm=norm;
		while (k<length) {
			Real *phi=Mol[i]->phi+k*M;
			//Real *G1=Seg[Mol[i]->MolMonList[k]]->G1;
			Real *G1=Mol[i]->G1+k*M;
			Div(phi,G1,M); if (norm>0) Norm(phi,norm,M);
if (debug) {
Real sum; Sum(sum,phi,M); cout <<"Sumphi in mol " << i << " for mon " << Mol[i]->MolMonList[k] << ": " << sum << endl;
}
			k++;
		}
		if (Mol[i]->freedom=="range_restricted") {
			Real *phit=Mol[i]->phitot; Zero(phit,M);
			int k=0;
			while (k<length) {
				Real *phi=Mol[i]->phi+k*M;
				Add(phit,phi,M);
				k++;
			}
			OverwriteA(phit,Mol[i]->R_mask,phit,M);
			//Lat[0]->remove_bounds(phit);
			Real theta=Lat[0]->ComputeTheta(phit);
			norm=Mol[i]->theta_range/theta;
			Mol[i]->norm = norm;
			Mol[i]->phibulk=Mol[i]->chainlength*norm; totphibulk +=Mol[i]->phibulk;
			Mol[i]->n=norm*Mol[i]->GN;
			Mol[i]->theta=Mol[i]->n*Mol[i]->chainlength;
			Zero(phit,M);
			k=0;
			while (k<length) {
				Real *phi=Mol[i]->phi+k*M;
				Norm(phi,norm,M);
if (debug) {
Real sum=Lat[0]->ComputeTheta(phi); cout <<"Sumphi in mol " << i << " for mon " << Mol[i]->MolMonList[k] << ": " << sum << endl;
}
				k++;
			}
		}
	}

	if (charged) {
		Real phib=0;
		for (int i=0; i<n_mol; i++) {
			if (i!=neutralizer && i!=solvent) phib+=Mol[i]->phibulk*Mol[i]->Charge();
		}
		Mol[neutralizer]->phibulk = -phib/Mol[neutralizer]->Charge();
		totphibulk+=Mol[neutralizer]->phibulk;
		norm = Mol[neutralizer]->phibulk/Mol[neutralizer]->chainlength;
		Mol[neutralizer]->n = norm*Mol[neutralizer]->GN;
		Mol[neutralizer]->theta = Mol[neutralizer]->n*Mol[neutralizer]->chainlength;
		Mol[neutralizer]->norm=norm;
	}

	Mol[solvent]->phibulk=1.0-totphibulk;
	norm=Mol[solvent]->phibulk/Mol[solvent]->chainlength;
	Mol[solvent]->n=norm*Mol[solvent]->GN;
	Mol[solvent]->theta=Mol[solvent]->n*Mol[solvent]->chainlength;
	Mol[solvent]->norm=norm;

	int k=0;
	length = Mol[solvent]->MolMonList.size();
	while (k<length) {
		Real *phi=Mol[solvent]->phi+k*M;
		if (norm>0) Norm(phi,norm,M);
if (debug) {
Real sum; Sum(sum,phi,M); cout <<"Sumphi in mol " << solvent << "for mon " << k << ":" << sum << endl;
}
		k++;
	}
	if (charged) {
		k=0;
		length = Mol[neutralizer]->MolMonList.size();
		while (k<length) {
			Real *phi=Mol[neutralizer]->phi+k*M;
			if (Mol[neutralizer]->norm>0) Norm(phi,Mol[neutralizer]->norm,M);
if (debug) {
Real sum; Sum(sum,phi,M); cout <<"Sumphi in mol " << neutralizer << "for mon " << k << ":" << sum << endl;
}
			k++;
		}
	}

	for (int i=0; i<n_mol; i++) {

		int length=Mol[i]->MolMonList.size();
		int k=0;
		while (k<length) {
			Real* phi_mon=Seg[Mol[i]->MolMonList[k]]->phi;
			Real* mol_phitot=Mol[i]->phitot;
			Real* phi_molmon = Mol[i]->phi + k*M;
			Add(phi_mon,phi_molmon,M);
			//if (!(Seg[Mol[i]->MolMonList[k]]->freedom == "tagged"))
Add(phitot,phi_molmon,M);
			Add(mol_phitot,phi_molmon,M);
			Seg[Mol[i]->MolMonList[k]]->phibulk +=Mol[i]->fraction(Mol[i]->MolMonList[k])*Mol[i]->phibulk;
			k++;
		}
		length=SysTagList.size();
		k=0;
		while (k<length) {
			Cp(Seg[SysTagList[k]]->phi,Seg[SysTagList[k]]->MASK,M);
			k++;
		}
		length=SysClampList.size();
		k=0;
		while (k<length) {
			Cp(Seg[SysClampList[k]]->phi,Seg[SysClampList[k]]->MASK,M);
			k++;
		}
	}
	int n_seg=In[0]->MonList.size();
	for (int i=0; i<n_seg; i++) {
		if (Seg[i]->freedom !="frozen") Lat[0]->set_bounds(Seg[i]->phi);
		Lat[0]->Side(Seg[i]->phi_side,Seg[i]->phi,M);
	}
	return success;
}

bool System::CheckResults(bool e_info) {
if (debug) cout << "CheckResults for system " << endl;
	bool success=true;
	FreeEnergy=GetFreeEnergy();
	GrandPotential=GetGrandPotential();
	CreateMu();
	if (e_info) cout <<endl;
	if (e_info) {
		cout <<"free energy                 = " << FreeEnergy << endl;
		cout <<"grand potential             = " << GrandPotential << endl;
	}
	int n_mol=In[0]->MolList.size();
	Real n_times_mu=0;
	for (int i=0; i<n_mol; i++) {
		Real Mu=Mol[i]->Mu;
		Real n=Mol[i]->n;
		n_times_mu +=  n*Mu;
	}
	if (e_info) {
		cout <<"free energy     (GP + n*mu) = " << GrandPotential + n_times_mu << endl;
		cout <<"Grand potential (F - n*mu)  = " << FreeEnergy - n_times_mu  << endl;
//	}

	//for (int i=0; i<n_mol; i++) { //NEED FIX . densities are not yet computed correctly that is why it is turned off.....!!!
		//if (Mol[i]->MolAlList.size()>0) {
		//	Mol[i]->compute_phi_alias=true;
			//Mol[i]->ComputePhi();
		//}
		//ComputePhis();
	//}
	cout << endl;
	int M=Lat[0]->M;
	for (int i=0; i<n_mol; i++) {
		int n_molmon=Mol[i]->MolMonList.size();
		Real theta_tot=Mol[i]->n*Mol[i]->chainlength;
		for (int j=0; j<n_molmon; j++) {
			Real FRACTION=Mol[i]->fraction(Mol[i]->MolMonList[j]);
			if (Seg[Mol[i]->MolMonList[j]]->freedom !="clamp" ) {
				Real THETA=Lat[0]->WeightedSum(Mol[i]->phi+j*M);
				cout <<"MOL " << Mol[i]->name << " Fraction " << Seg[Mol[i]->MolMonList[j]]->name << ": " << FRACTION << "=?=" << THETA/theta_tot << " or " << THETA << " of " << theta_tot <<  endl;
			}
		}
	}
	cout << endl;
}
	return success;
}

Real System::GetFreeEnergy(void) {
if (debug) cout << "GetFreeEnergy for system " << endl;
	int M=Lat[0]->M;
	Real FreeEnergy=0;
	Real* F=FreeEnergyDensity;
	Real constant=0;
	int n_mol=In[0]->MolList.size();
	//for (int i=0; i<n_mol; i++) Lat[0]->remove_bounds(Mol[i]->phitot);
	int n_mon=In[0]->MonList.size();
	for (int i=0; i<n_mon; i++) {
		Lat[0]->remove_bounds(Seg[i]->phi_side);
	}

	Zero(F,M);
	if (charged) {
		Times(F,EE,eps,M); Norm(F,-2.0,M); //EE bevat een deling door 2 en de norm met -2 is two pre-correct de latere deling door 2.
		AddTimes(F,q,psi,M);
		Norm(F,0.5,M);
	}

	for (int i=0; i<n_mol; i++){
		if (Mol[i]->freedom =="clamped") {
			int n_box=Mol[i]->n_box;
			for (int p=0; p<n_box; p++) {
				FreeEnergy-=log(Mol[i]->gn[p]);
			}
		} else {
			Real n=Mol[i]->n;
			Real GN=Mol[i]->GN;
			int N=Mol[i]->chainlength;
			if (Mol[i]->IsTagged()) N--; //assuming there is just one tagged segment per molecule
			Real *phi=Mol[i]->phitot; //contains also the tagged segment
			constant = log(N*n/GN)/N;
			Cp(TEMP,phi,M); Norm(TEMP,constant,M); Add(F,TEMP,M);
		}
	}

	int n_sysmon=SysMonList.size();
	for (int j=0; j<n_sysmon; j++) {
		Real* phi=Seg[SysMonList[j]]->phi;
		Real* u=Seg[SysMonList[j]]->u;
		Times(TEMP,phi,u,M); Norm(TEMP,-1,M); Add(F,TEMP,M);
	}
	for (int j=0; j<n_mon; j++) for (int k=0; k<n_mon; k++) {
		Real chi;
		if (Seg[k]->freedom=="frozen") chi=CHI[j*n_mon+k]; else  chi = CHI[j*n_mon+k]/2;
		Real *phi_side = Seg[k]->phi_side;
		Real *phi = Seg[j]->phi;
		if (Seg[j]->freedom !="frozen") {Times(TEMP,phi,phi_side,M); Norm(TEMP,chi,M); Add(F,TEMP,M);}
	}

	for (int i=0; i<n_mol; i++){
		constant=0;
		int n_molmon=Mol[i]->MolMonList.size();
		for (int j=0; j<n_molmon; j++) for (int k=0; k<n_molmon; k++) {
			Real fA=Mol[i]->fraction(Mol[i]->MolMonList[j]);
			Real fB=Mol[i]->fraction(Mol[i]->MolMonList[k]);
			if (Mol[i]->IsTagged()) {
				int N=Mol[i]->chainlength;
				if (N>1) {fA=fA*N/(N-1); fB=fB*N/(N-1); } else {fA =0; fB=0;}
			}
			Real chi = CHI[Mol[i]->MolMonList[j]*n_mon+Mol[i]->MolMonList[k]]/2;
			constant -=fA*fB*chi;
		}
		Real* phi=Mol[i]->phitot;
		Cp(TEMP,phi,M); Norm(TEMP,constant,M); Add(F,TEMP,M);
	}
	Lat[0]->remove_bounds(F); Times(F,F,KSAM,M);
	return FreeEnergy+Lat[0]->WeightedSum(F);
}

Real System::GetGrandPotential(void) {
if (debug) cout << "GetGrandPotential for system " << endl;
	int M=Lat[0]->M;
	Real* GP =GrandPotentialDensity;
	int n_mol=In[0]->MolList.size();
	int n_mon=In[0]->MonList.size();
	Zero(GP,M);
	if (charged) {Times(GP,q,psi,M); Norm(GP,0.5,M);}
	for (int i=0; i<n_mol; i++){
		Real *phi=Mol[i]->phitot;
		Real phibulk = Mol[i]->phibulk;
		int N=Mol[i]->chainlength;
		if (Mol[i]->IsTagged()) N--; //One segment of the tagged molecule is tagged and then removed from GP through KSAM
		if (Mol[i]->IsClamped()) N=N-2;
		Cp(TEMP,phi,M); YisAplusC(TEMP,TEMP,-phibulk,M); Norm(TEMP,1.0/N,M); //GP has wrong sign. will be corrected at end of this routine;
		Add(GP,TEMP,M);
	}
	Add(GP,alpha,M);
	int n_sysmon=SysMonList.size();
	for (int j=0; j<n_sysmon; j++)for (int k=0; k<n_sysmon; k++){
		if (!(Seg[SysMonList[j]]->freedom=="tagged" || Seg[SysMonList[j]]->freedom=="clamp"  )) { //not sure about this line...
		Real phibulkA=Seg[SysMonList[j]]->phibulk;
		Real phibulkB=Seg[SysMonList[k]]->phibulk;
		Real chi = CHI[SysMonList[j]*n_mon+SysMonList[k]]/2;
		Real *phi=Seg[SysMonList[j]]->phi;
		Real *phi_side=Seg[SysMonList[k]]->phi_side;
		Times(TEMP,phi,phi_side,M); YisAplusC(TEMP,TEMP,-phibulkA*phibulkB,M); Norm(TEMP,chi,M); Add(GP,TEMP,M);
	}
	}
	Norm(GP,-1.0,M); //correct the sign.
	Lat[0]->remove_bounds(GP); Times(GP,GP,KSAM,M);
	return  Lat[0]->WeightedSum(GP);

}

bool System::CreateMu() {
if (debug) cout << "CreateMu for system " << endl;
	bool success=true;
	Real constant;
	Real n;
	Real GN;
	int n_mol=In[0]->MolList.size();
	int n_mon=In[0]->MonList.size();
	for (int i=0; i<n_mol; i++) {
		Real Mu=0;
		Real NA=Mol[i]->chainlength;
		if (Mol[i]->IsTagged()) NA=NA-1;
		if (Mol[i]->IsClamped()) NA=NA-2;
		if (Mol[i]->IsClamped()) {
			n=Mol[i]->n_box;
			int n_box = Mol[i]->n_box;
			GN=0;
			for (int p=0; p<n_box; p++) {
				GN+=log(Mol[i]->gn[p]);
			}
			Mu= -GN/n +1;
		} else {
			n=Mol[i]->n;
			GN=Mol[i]->GN;
			Mu=log(NA*n/GN)+1;
		}
		constant=0;
		for (int k=0; k<n_mol; k++) {
			Real NB = Mol[k]->chainlength;
			if (Mol[k]->IsTagged()) NB=NB-1;
			if (Mol[k]->IsClamped()) NB=NB-2;
			Real phibulkB=Mol[k]->phibulk;
			constant +=phibulkB/NB;
		}
		Mu = Mu - NA*constant;

		for (int j=0; j<n_mon; j++) for (int k=0; k<n_mon; k++) {
			Real chi= CHI[j*n_mon+k]/2;
			Real phibulkA=Seg[j]->phibulk;
			Real phibulkB=Seg[k]->phibulk;
			Real Fa=Mol[i]->fraction(j);
			Real Fb=Mol[i]->fraction(k);
			if (Mol[i]->IsTagged()) {Fa=Fa*(NA+1)/(NA); Fb=Fb*(NA+1)/NA;}
			if (Mol[i]->IsClamped()) {Fa=Fa*(NA+2)/(NA); Fb=Fb*(NA+2)/NA;}
			Mu = Mu-NA*chi*(phibulkA-Fa)*(phibulkB-Fb);
		}

		Mol[i]->Mu=Mu;
//cout <<"mol" << i << " = " << Mu << endl;
	}
	return success;
}

//TODO //make sure that the tagged positions are not in frozen range.
