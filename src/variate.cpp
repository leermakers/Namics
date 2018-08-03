#include "variate.h"

Variate::Variate(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, vector<Molecule*> Mol_,vector<System*>Sys_, string name_) {
	In=In_; name=name_;   Lat=Lat_; Seg=Seg_, Mol=Mol_; Sys=Sys_;
	KEYS.push_back("scan");
	KEYS.push_back("search");
	KEYS.push_back("step");
	KEYS.push_back("steps");
	KEYS.push_back("end_value");
	KEYS.push_back("scale");
	KEYS.push_back("grand_potential");
	KEYS.push_back("free_energy");
	KEYS.push_back("mu");
	KEYS.push_back("theta");
	KEYS.push_back("n");
	KEYS.push_back("phibulk");
}
Variate::~Variate() {
	DeAllocateMemory();
}
void Variate::DeAllocateMemory(){
if (debug) cout <<"Destructor for variate " + name << endl;

#ifdef CUDA

#else

#endif
}

void Variate::AllocateMemory() {
#ifdef CUDA

#else

#endif

}

void Variate::PrepareForCalculations() {
if (debug) cout <<"PrepareForCalculations in Variate " + name << endl;


}


void Variate::PutParameter(string new_param) {
if (debug) cout <<"PutParameter in Variate " + name << endl;
	KEYS.push_back(new_param);
}

bool Variate::CheckInput(int start) {
if (debug) cout <<"CheckInput in Variate " + name << endl;
	bool success=true;
	num_of_cals=0;
	success= In[0]->CheckParameters("var",name,start,KEYS,PARAMETERS,VALUES);
	if (success && name != "noname") {
		vector<string>sub;
		In[0]->split(name,'-',sub);
		if (sub.size()!=2) {
			cout <<"in 'var' the second argument should be composed of two parts separated by a '-'. The first element is an 'item'={sys, mol, mon, lat}. The second element is a valid corresponding 'name'. For example: mol-A . " <<endl;
			return false;
		}
		int pos=0;
		Real R_target= -123.0;
		bool dubbel=false;
		int choice=4;
		scanning =-1;
		targeting =-1;
		searching=-1;
		eq_to_solvating=-1;
		eq_to_mu =-1;
		ets_nr=-1;
		scan_nr=0;
		search_nr=0;
		if (sub[0]=="sys") choice=0;
		if (sub[0]=="lat") choice=1;
		if (sub[0]=="mol") choice=2;
		if (sub[0]=="mon") choice=3;
		switch(choice) {
			case 0:
				if (!In[0]->InSet(In[0]->SysList,pos,sub[1])) {
					cout <<"In 'var' sys name " + sub[1] + " not found." << endl; success=false;
				} else {
					if (GetValue("free_energy").size()>0) {
						R_target=In[0]->Get_Real(GetValue("free_energy"),0); targeting=0; target_nr=0;
						if (!Sys[0]->PutVarInfo("target","free_energy",R_target)) {
							success=false;
							cout <<"In var:" + name + ":free_energy, Target value rejected" << endl;
						}
					}
					if (GetValue("grand_potential").size()>0) {
						if (R_target ==-123.0) {
							R_target=In[0]->Get_Real(GetValue("grand_potential"),0); targeting=0; target_nr=0;
							if (!Sys[0]->PutVarInfo("target","grand_potential",R_target)) {
								success=false;
								cout <<"In var:" + name + ":grand_potential, target value rejected." << endl;
							}
						} else {dubbel=true;}

					}
					if (R_target==-123.0 || dubbel) {
						success=false;
						cout <<"In var:" + name + " we expect either 'free_energy' or 'grand_potential' as target a function. " << endl;
					}
				}
				break;
			case 1:
				if (!In[0]->InSet(In[0]->LatList,pos,sub[1])) {
					cout <<"In 'var' lat name " + sub[1] + " not found" << endl; success=false;
				} else {
					if (GetValue("scan").size()==0) {
						success=false; cout <<"in var:" + name + ", the parameter 'scan' is expected. " << endl;
						if (GetValue("search").size()>0) cout <<"In var:" + name + ", the parameter 'search' is not allowed. " <<endl;
					} else {
						if (!Lat[0]->PutVarInfo("scan",GetValue("scan"),0)) {
							success=false; cout <<"In var:" + name + "scan, the value for 'target' is rejected" << endl;
						}
						scanning=1; scan_nr=0;
					}
				}
				break;
			case 2:
				if (!In[0]->InSet(In[0]->MolList,pos,sub[1])) {
					cout <<"In 'var' mol name " + sub[1] + " not found" << endl; success=false;
				} else {
					if (GetValue("scan").size()==0 && GetValue("search").size() ==0) {
						if (GetValue("mu").size()>0) {
							string s = GetValue("mu");
							vector<string>sub_;
							In[0]->split(s,'-',sub_);
							if (sub_[0]==s) {
								R_target=In[0]->Get_Real(GetValue("mu"),0); targeting=2; target_nr=pos;
								if (!Mol[pos]->PutVarInfo("target","mu",R_target)) {
									success=false; cout <<"In var:"+name+":mu, Target value rejected" << endl;
								}
							} else {
//cout <<"punt 1 "<< endl;
								target_nr = pos; R_target=12345.0;
								if (sub_[0]!="mol") {success=false; cout <<"in var:"+name+":mu, target should be specified as mol-'name'. " << endl; }
								if (!In[0]->InSet(In[0]->MolList,eq_to_mu,sub_[1])) { success=false; cout <<"In 'var:"+name+":mu: mol-'name', 'name' is not valid mol-name. " << endl; }
								if (success) Mol[pos]->PutVarInfo("target","mu",R_target);
							}
//cout <<"punt 2 " << endl;
						}
						if (GetValue("theta").size()>0) {
							if (R_target ==-123.0) {
								R_target=In[0]->Get_Real(GetValue("theta"),0); targeting=2; target_nr=pos;
								if (!Mol[0]->PutVarInfo("target","theta",R_target)) {
									success=false;
									cout <<"In var:"+name+":theta, Target value rejected" << endl;
								}
							} else {dubbel = true;}
						}
						if (GetValue("n").size()>0) {
							if (R_target ==-123.0) {
								R_target=In[0]->Get_Real(GetValue("n"),0); targeting=2; target_nr=pos;
								if (!Mol[0]->PutVarInfo("target","n",R_target)) {
									success=false;
									cout <<"In var:"+name+":n, Target value rejected" << endl;
								}
							} else {dubbel = true;}
						}
						if (GetValue("theta").size()>0) {
							if (R_target ==-123.0) {
								R_target=In[0]->Get_Real(GetValue("phibulk"),0); targeting=2; target_nr=pos;
								if (!Mol[0]->PutVarInfo("target","phibulk",R_target)) {
									success=false;
									cout <<"In var:"+name+":phibulk, Target value rejected" << endl;
								}
							} else {dubbel = true;}
						}
						if (R_target == -123.0) {
							success=false;
							cout <<"In var:" + name + ", we need either 'scan' or 'search' as third parameter." << endl;
							cout <<"or targets selected from {mu, theta, n, phibulk} "  << endl;
						}
					} else {
						if (GetValue("scan").size()>0) {
							if (!Mol[pos]->PutVarInfo("scan",GetValue("scan"),0)) {
								success=false; cout <<"In var:" + name + ":scan, the target is rejected " << endl;
							} else {scanning = 2; scan_nr=pos;}
						}
						if (GetValue("search").size()>0) {
							if (!Mol[pos]->PutVarInfo("search",GetValue("search"),0)) {
								success=false; cout <<"In var:" + name + ":search, the target is rejected " << endl;
							} else {
								if (Mol[pos]->Var_search_value==3) {
									eq_to_solvating=2; ets_nr=pos;
								} else {
									searching=2; search_nr=pos;
								}
							}

						}
					}
				}
				break;
			case 3:
				if (!In[0]->InSet(In[0]->MonList,pos,sub[1])) {
					cout <<"in 'var' mon name " + sub[1] + " not found" << endl; success=false;
				} else {
					if (GetValue("scan").size()==0) {
						success=false;
						cout <<"In var:" + name + ", we need a 'scan' parameter." << endl;
					} else {
						if (GetValue("scan").size()>0) {
							scanning=3; scan_nr=pos;
							if (!Seg[pos]->PutVarInfo("scan",GetValue("scan"),0)) {
								success=false; cout <<"In var:" + name + ":scan, the target is rejected " << endl;
							}
						}
					}
				}
				break;
			default:
				success=false;
				cout <<"In var: keyword info from 'name' " + name + " is not recognised. Choose from 'sys, lat, mol, mon' " <<endl;
				break;
		}
		if (scanning>-1) {
			scale="";
			if (GetValue("end_value").size() == 0) {
				success=false; cout <<"In var: the 'scan' item is set and therefore the system expect an 'end_value'"<<endl;
			}  else {end_value=In[0]->Get_Real(GetValue("end_value"),0);}
			if (GetValue("scale").size() == 0) {scale="linear";} else {scale=GetValue("scale"); }
			if (scale =="exponential") {
				if (GetValue("steps").size() ==0) {
					steps = 1; cout <<"In var: the property 'steps' is set to the default value of '1'"<< endl;
				} else {
					steps = In[0]->Get_int(GetValue("steps"),1);
					if (steps < 1) {
						success = false;
						cout <<"In var: the property 'steps' should be a positive integer, indicating the number of 'steps' per 'decade' " << endl;
					}
				}
			} else {
				if (scale=="linear"){
					if (GetValue("step").size() == 0) {success=false; cout <<"In var: while issueing a 'scan' you need to supply the value for 'step'. "<< endl; }
					else {
						step=In[0]->Get_Real(GetValue("step"),0);
						if (step==0) {success=false; cout <<"In var: while issuing a 'scan' the value for 'step' is not recognised or equal to zero. " << endl;}
					}
					if (GetValue("steps").size()!=0) {
						cout <<"in var: while issuing a 'scan' the value for 'steps' is ignored. You are advised to use 'step' " << endl;
					}
				} else {
					success=false; cout <<"In var: property 'scale' should be either 'exponential' or 'linear' and not " + scale << endl;
				}
			}
			if (success) {
				switch(scanning) {
					case 0:
						cout <<"programming error in variate" <<endl;
						break;
					case 1:
						num_of_cals=Lat[0]->PutVarScan(step,end_value);
						break;
					case 2:
						num_of_cals=Mol[scan_nr]->PutVarScan(step,end_value,steps,scale);
						break;
					case 3:
						num_of_cals=Seg[scan_nr]->PutVarScan(step,end_value,steps,scale);
						break;
					default:
						cout <<"progamming error in variate " << endl;
						break;
				}
				success=num_of_cals>0;
				cout <<"num of calculations " << num_of_cals << endl;

			}
		} else {
			if (GetValue("step").size()>0) { cout <<"In var: the 'scan' property is not set and therefore the 'step'-variable is ignored. " << endl; }
			if (GetValue("steps").size()>0) { cout <<"In var: the 'scan' property is not set and therefore the 'steps'-variable is ignored. " << endl; }
			if (GetValue("end_value").size()>0) { cout <<"In var: the 'scan' property is not set and therefore the 'end_value'-variable is ignored. " << endl; }
			if (GetValue("scale").size()>0) { cout <<"In var: the 'scan' property is not set and therefore the 'scale'-variable is ignored. " << endl; }
		}
		if (scanning>0 && searching>0) {
			if (scanning==2 && searching==2) {
				if (scan_nr==search_nr) {
					if (Mol[scan_nr]->Var_scan_value==Mol[search_nr]->Var_search_value) {
						cout <<"In var: The search quantity can not be equal to the scan quantity " << endl;
						success=false;
					}
				}
			}
		}

	}
	return success;
}

bool Variate::PutVarScan(int cal_nr) {
	bool success=true;
	switch(scanning) {
		case 0:
			success=false;
			cout <<"programming error in PutVarScan" << endl;
			break;
		case 1:
			Lat[scan_nr]->UpdateVarInfo(cal_nr);
			break;
		case 2:
			Mol[scan_nr]->UpdateVarInfo(cal_nr);
			break;
		case 3:
			Seg[scan_nr]->UpdateVarInfo(cal_nr);
			if (Seg[scan_nr]->chi_var_seg>0) {
				int n_seg=In[0]->MonList.size();
				Sys[0]->CHI[scan_nr*n_seg+Seg[scan_nr]->chi_var_seg] = Seg[scan_nr]->chi_value;
				Sys[0]->CHI[scan_nr+n_seg*Seg[scan_nr]->chi_var_seg] =Sys[0]->CHI[scan_nr*n_seg+Seg[scan_nr]->chi_var_seg];
			} else Seg[scan_nr]->ResetInitValue();
			break;
		default:
			break;
	}
	return success;
}

void Variate::PutValue(Real X) {
	switch(searching) {
		case 0:
			cout <<"programming error in PutValue" << endl;
			break;
		case 1:
			Seg[search_nr]->PutValue(X);
			break;
		case 2:
			Mol[search_nr]->PutValue(X);
			break;
		case 3:
			cout <<"programming error in PutValue" << endl;
			break;
		default:
			break;
	}
	if (ets_nr>-1) {Mol[ets_nr]->theta=X; Mol[ets_nr]->n=X/Mol[ets_nr]->chainlength;}
	if (eq_to_mu>-1) {Mol[eq_to_mu]->theta=X; Mol[eq_to_mu]->n=X/Mol[eq_to_mu]->chainlength;}
}
Real Variate::GetValue(void) {
	Real X=0;
	switch(searching) {
		case 0:
			cout <<"programming error in GetValue" << endl;
			break;
		case 1:
			X=Seg[search_nr]->GetValue();
			break;
		case 2:
			if (eq_to_mu>-1) X=Mol[eq_to_mu]->GetValue(); else
			X=Mol[search_nr]->GetValue();
			break;
		case 3:
			cout <<"programming error in GetValue" << endl;
			break;
		default:
			break;
	}
	if (ets_nr>-1) X=Mol[ets_nr]->theta;
	if (eq_to_mu>-1) X=Mol[eq_to_mu]->theta;
	return X;
}

Real Variate::GetError(void) {
	Real X=0;
	switch(targeting) {

		case 0:
			X=Sys[target_nr]->GetError();
			break;
		case 1:
			cout <<"programming error in GetError" << endl;
			break;
		case 2:
			if (target_nr>-1) X=Mol[target_nr]->GetError();
			break;
		case 3:
			cout <<"programming error in GetError" << endl;
			break;
		default:
			break;
	}
	if (ets_nr>-1) X=-(Mol[ets_nr]->theta/Sys[0]->Mol[Sys[0]->solvent]->theta-1.0);
	if (eq_to_mu>-1) X=(Mol[eq_to_mu]->Mu-Mol[target_nr]->Mu)*-1.0;
	return X;
}

bool Variate::ResetScanValue(void) {
	if (debug) cout <<"ResetScanValue in Variate" << endl;
	bool success=true;
	switch(scanning) {
		case 0:
			success=false;
			cout <<"programming error in ResetVarScan" << endl;
			break;
		case 1:
			Lat[scan_nr]->ResetInitValue();
			break;
		case 2:
			Mol[scan_nr]->ResetInitValue();
			break;
		case 3:
			if (Seg[scan_nr]->chi_var_seg>0) {
				int n_seg=In[0]->MonList.size();
				Sys[0]->CHI[scan_nr*n_seg+Seg[scan_nr]->chi_var_seg] = Seg[scan_nr]->Var_start_value;
				Sys[0]->CHI[scan_nr+n_seg*Seg[scan_nr]->chi_var_seg] =Sys[0]->CHI[scan_nr*n_seg+Seg[scan_nr]->chi_var_seg];
			} else Seg[scan_nr]->ResetInitValue();
			break;
		default:
			cout <<"programming error in ResetVarScan" << endl;
			break;
	}
	return success;
}

string Variate::GetValue(string parameter){
if (debug) cout <<"GetValue in Variate " + name << " " <<parameter << endl;
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

void Variate::push(string s, Real X) {
if (debug) cout <<"push (Real) in Variate " + name << endl;
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void Variate::push(string s, int X) {
if (debug) cout <<"push (int) in Variate " + name << endl;
	ints.push_back(s);
	ints_value.push_back(X);
}
void Variate::push(string s, bool X) {
if (debug) cout <<"push (boool) in Variate " + name << endl;
	bools.push_back(s);
	bools_value.push_back(X);
}
void Variate::push(string s, string X) {
if (debug) cout <<"push (string) in Variate " + name << endl;
	strings.push_back(s);
	strings_value.push_back(X);
}
void Variate::PushOutput() {
if (debug) cout <<"PushOutput in Variate " + name << endl;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
#ifdef CUDA
	TransferDataToHost(H_phi,phi,M);
#endif
}

Real* Variate::GetPointer(string s) {
if (debug) cout <<"GetPointer in Variate " + name << endl;
	return NULL;
}


int Variate::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
if (debug) cout <<"GetValue (long)  in Variate " + name << endl;
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
