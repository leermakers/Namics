#include "cleng.h"
#include "output.h"

Cleng::Cleng(vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Newton*> New_, vector<Engine*> Eng_,  string name_)
    : name{name_},
      In{In_},
      Lat{Lat_},
      Mol{Mol_},
      Seg{Seg_},
      Sys{Sys_},
      New{New_},
      Eng{Eng_}

{
  KEYS.push_back("MCS");
  KEYS.push_back("save_interval");
  KEYS.push_back("save_filename");
  KEYS.push_back("seed");
  if (debug)
    cout << "Cleng initialized" << endl;


}

Cleng::~Cleng() {
}

bool Cleng::CheckInput(int start) {
  if (debug)
    cout << "CheckInput in Cleng" << endl;
  bool success = true;


  /* When adding new items here, please add them to the function prepareOutputFile()*/
  success = In[0]->CheckParameters("cleng", name, start, KEYS, PARAMETERS, VALUES);
  if (success) {
    vector<string> options;
    if (GetValue("MCS").size() > 0) {
      success = In[0]->Get_int(GetValue("MCS"), MCS, 1, 10000, "The number of timesteps should be between 1 and 10000");
    }
    if (debug)
      cout << "MCS is " << MCS << endl;

    if (GetValue("save_interval").size() > 0) {
      success = In[0]->Get_int(GetValue("save_interval"), save_interval,1,100,"The save interval nr should be between 1 and 100. ");
    }
    if (debug) cout << "Save_interval " << save_interval << endl;
    if (Sys[0]->SysClampList.size() <1) {cout <<"Cleng needs to have clamped molecules in the system" << endl; success=false;}
	else {clamp_seg=Sys[0]->SysClampList[0]; if (Sys[0]->SysClampList.size()>1) {success=false; cout <<"Currently the clamping is limited to one molecule per system. " << endl; }}
//cout <<"clamp_seg : " << clamp_seg << endl; 
	if (success) {
		n_boxes = Seg[clamp_seg]->n_box;
		sub_box_size=Seg[clamp_seg]->mx;
//cout <<"subboxsize " << sub_box_size << endl; 
		//success=CP(to_cleng);
		//success=CP(to_segment);
		New[0]->Solve(true);  
	}
	
  }

   int n_out = In[0]->OutputList.size();
    if (n_out == 0)
      cout << "Warning: no output defined!" << endl;

    // Create output class instance and check inputs (reference above)
    for (int i = 0; i < n_out; i++) {

      Out.push_back(new Output(In, Lat, Seg, Mol, Sys, New, Eng, In[0]->OutputList[i], i, n_out));
      if (!Out[i]->CheckInput(start)) {
        cout << "input_error in output " << endl;
        success=false;
      }
    }

  return success;
}


bool Cleng::CP(transfer tofrom) {
	int MX=Lat[0]->MX;
	int MY=Lat[0]->MY;
	int MZ=Lat[0]->MZ;
	int JX=Lat[0]->JX;
	int JY=Lat[0]->JY;
	int M = Lat[0]->M;
	bool success=true;
	int i,j;
	int length; 
	bool found;
	switch(tofrom) {
		case to_cleng:
			for (i=0; i<n_boxes; i++) {
				PX=Seg[clamp_seg]->px1[i]; if (PX>MX) {PX-=MX; Sx.push_back(1);} else {Sx.push_back(0);}
				PY=Seg[clamp_seg]->py1[i]; if (PY>MY) {PY-=MY; Sy.push_back(1);} else {Sy.push_back(0);}
				PZ=Seg[clamp_seg]->pz1[i]; if (PZ>MZ) {PZ-=MZ; Sz.push_back(1);} else {Sz.push_back(0);}
				length=X.size();
				found=false; j=0;
				while (!found && j<length) {
					if (X[j]==PX&&Y[j]==PY && Z[j]==PZ) {found = true; P.push_back(j);} else j++;
				}
				if (found==false) {
					X.push_back(PX);Y.push_back(PY);Z.push_back(PZ); P.push_back(X.size()-1);
				}
//cout << "("<<PX<<","<<PY<<","<<PZ<<")"; 
				PX=Seg[clamp_seg]->px2[i]; if (PX>MX) {PX-=MX; Sx.push_back(1);} else {Sx.push_back(0);}
				PY=Seg[clamp_seg]->py2[i]; if (PY>MY) {PY-=MY; Sy.push_back(1);} else {Sy.push_back(0);}
				PZ=Seg[clamp_seg]->pz2[i]; if (PZ>MZ) {PZ-=MZ; Sz.push_back(1);} else {Sz.push_back(0);}
				length=X.size();
				found=false; j=0;
				while (!found&& j<length) {
					if (X[j]==PX&&Y[j]==PY && Z[j]==PZ) {found = true; P.push_back(j);} else j++;
				}
				if (found==false) {
					X.push_back(PX);Y.push_back(PY);Z.push_back(PZ); P.push_back(X.size()-1);
				}
				
//cout << "("<<PX<<","<<PY<<","<<PZ<<")"<< endl; 

			}
//cout <<"number of unique pos " << X.size() << endl;
//cout <<"size of P " << P.size()<< endl; 
 
		break;
		case to_segment:
			Zero(Seg[clamp_seg]->H_MASK,M); 
			for (i=0; i<n_boxes; i++) {
				Seg[clamp_seg]->px1[i]=X[P[2*i]]+Sx[2*i]*MX;
				Seg[clamp_seg]->px2[i]=X[P[2*i+1]]+Sx[2*i+1]*MX;
				Seg[clamp_seg]->py1[i]=Y[P[2*i]]+Sy[2*i]*MY;
				Seg[clamp_seg]->py2[i]=Y[P[2*i+1]]+Sy[2*i+1]*MY;
				Seg[clamp_seg]->pz1[i]=Z[P[2*i]]+Sz[2*i]*MZ;
				Seg[clamp_seg]->pz2[i]=Z[P[2*i+1]]+Sz[2*i+1]*MZ;
				Seg[clamp_seg]->bx[i]=(Seg[clamp_seg]->px2[i]+Seg[clamp_seg]->px1[i]-sub_box_size)/2;
				Seg[clamp_seg]->by[i]=(Seg[clamp_seg]->py2[i]+Seg[clamp_seg]->py1[i]-sub_box_size)/2;
				Seg[clamp_seg]->bz[i]=(Seg[clamp_seg]->pz2[i]+Seg[clamp_seg]->pz1[i]-sub_box_size)/2;
//cout << "("<<Seg[clamp_seg]->bx[i]<<","<<Seg[clamp_seg]->by[i]<<","<<Seg[clamp_seg]->bz[i]<<")";
				if (Seg[clamp_seg]->bx[i]<1) {Seg[clamp_seg]->bx[i] +=MX; Seg[clamp_seg]->px1[i] +=MX; Seg[clamp_seg]->px2[i] +=MX;} 
				if (Seg[clamp_seg]->by[i]<1) {Seg[clamp_seg]->by[i] +=MY; Seg[clamp_seg]->py1[i] +=MY; Seg[clamp_seg]->py2[i] +=MY;} 
				if (Seg[clamp_seg]->bz[i]<1) {Seg[clamp_seg]->bz[i] +=MZ; Seg[clamp_seg]->pz1[i] +=MZ; Seg[clamp_seg]->pz2[i] +=MZ;} 
				Seg[clamp_seg]->H_MASK[((Seg[clamp_seg]->px1[i]-1)%MX+1)*JX + ((Seg[clamp_seg]->py1[i]-1)%MY+1)*JY + (Seg[clamp_seg]->pz1[i]-1)%MZ+1]=1;
				Seg[clamp_seg]->H_MASK[((Seg[clamp_seg]->px2[i]-1)%MX+1)*JX + ((Seg[clamp_seg]->py2[i]-1)%MY+1)*JY + (Seg[clamp_seg]->pz2[i]-1)%MZ+1]=1;
//cout << "("<<Seg[clamp_seg]->bx[i]<<","<<Seg[clamp_seg]->by[i]<<","<<Seg[clamp_seg]->bz[i]<<")";
//cout << "("<<Seg[clamp_seg]->px1[i]<<","<<Seg[clamp_seg]->py1[i]<<","<<Seg[clamp_seg]->pz1[i]<<")("<<Seg[clamp_seg]->px2[i]<<","<<Seg[clamp_seg]->py2[i]<<","<<Seg[clamp_seg]->pz2[i]<<")"<<endl;
			}
		break;
		default:
			success=false; 
			cout <<"error in tranfer" << endl; 
		break;
	}
	return success; 
}


bool Cleng::MonteCarlo() {
  if (debug)
    cout << "Monte Carlo in Cleng" << endl;

/*
     Lat[0]->PushOutput();
      New[0]->PushOutput();
      Eng[0]->PushOutput();
      int length = In[0]->MonList.size();
      for (int i = 0; i < length; i++)
        Seg[i]->PushOutput();
      length = In[0]->MolList.size();
      for (int i = 0; i < length; i++) {
        int length_al = Mol[i]->MolAlList.size();
        for (int k = 0; k < length_al; k++) {
          Mol[i]->Al[k]->PushOutput();
        }
	PushOutput();
        Mol[i]->PushOutput();
      }
      // length = In[0]->AliasList.size();
      // for (int i=0; i<length; i++) Al[i]->PushOutput();
      Sys[0]->PushOutput(); // needs to be after pushing output for seg.
*/
	  
	return true;
}


void Cleng::prepareOutputFile() {
  /* Open filestream and set filename to "mesodyn-datetime.csv" */

 }


void Cleng::PutParameter(string new_param) {
  KEYS.push_back(new_param);
}
string Cleng::GetValue(string parameter) {
  int i = 0;
  int length = PARAMETERS.size();
  while (i < length) {
    if (parameter == PARAMETERS[i]) {
      return VALUES[i];
    }
    i++;
  }
  return "";
}
void Cleng::push(string s, Real X) {
  Reals.push_back(s);
  Reals_value.push_back(X);
}
void Cleng::push(string s, int X) {
  ints.push_back(s);
  ints_value.push_back(X);
}
void Cleng::push(string s, bool X) {
  bools.push_back(s);
  bools_value.push_back(X);
}
void Cleng::push(string s, string X) {
  strings.push_back(s);
  strings_value.push_back(X);
}
void Cleng::PushOutput() {
  strings.clear();
  strings_value.clear();
  bools.clear();
  bools_value.clear();
  Reals.clear();
  Reals_value.clear();
  ints.clear();
  ints_value.clear();
}

int Cleng::GetValue(string prop, int& int_result, Real& Real_result, string& string_result) {
  int i = 0;
  int length = ints.size();
  while (i < length) {
    if (prop == ints[i]) {
      int_result = ints_value[i];
      return 1;
    }
    i++;
  }
  i = 0;
  length = Reals.size();
  while (i < length) {
    if (prop == Reals[i]) {
      Real_result = Reals_value[i];
      return 2;
    }
    i++;
  }
  i = 0;
  length = bools.size();
  while (i < length) {
    if (prop == bools[i]) {
      if (bools_value[i])
        string_result = "true";
      else
        string_result = "false";
      return 3;
    }
    i++;
  }
  i = 0;
  length = strings.size();
  while (i < length) {
    if (prop == strings[i]) {
      string_result = strings_value[i];
      return 3;
    }
    i++;
  }
  return 0;
}
