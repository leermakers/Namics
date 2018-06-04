#include "teng.h"
#include "output.h"

Teng::Teng(vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Newton*> New_,  string name_)
    : name{name_},
      In{In_},
      Lat{Lat_},
      Mol{Mol_},
      Seg{Seg_},
      Sys{Sys_},
      New{New_}

{
	if (debug) cout << "Teng initialized" << endl;
 	 KEYS.push_back("MCS");
  	KEYS.push_back("save_interval");
  	KEYS.push_back("save_filename");
  	KEYS.push_back("seed");
}

Teng::~Teng() {
}

bool Teng::CheckInput(int start) {
  if (debug) cout << "CheckInput in Teng" << endl;
  bool success = true;

  success = In[0]->CheckParameters("cleng", name, start, KEYS, PARAMETERS, VALUES);
  if (success) {
    vector<string> options;
    if (GetValue("MCS").size() > 0) {
      success = In[0]->Get_int(GetValue("MCS"), MCS, 1, 10000, "The number of timesteps should be between 1 and 10000");
    }
    if (debug)
      cout << "MCS is " << MCS << endl;

    if (GetValue("save_interval").size() > 0) {
      success = In[0]->Get_int(GetValue("save_interval"), save_interval,1,MCS,"The save interval nr should be between 1 and 100");
    }
    if (debug) cout << "Save_interval " << save_interval << endl;
    if (Sys[0]->SysClampList.size() <1) {cout <<"Teng needs to have clamped molecules in the system" << endl; success=false;}
	else {clamp_seg=Sys[0]->SysClampList[0]; if (Sys[0]->SysClampList.size()>1) {success=false; cout <<"Currently the clamping is limited to one molecule per system. " << endl; }}
    if (success) {
	n_boxes = Seg[clamp_seg]->n_box;
	sub_box_size=Seg[clamp_seg]->mx;
    }
    clp_mol=-1;
    int length = In[0]->MolList.size();
    for (int i=0; i<length; i++) if (Mol[i]->freedom =="clamped") clp_mol=i; 
  }
  if (success) {
    n_out = In[0]->OutputList.size();
    if (n_out == 0) cout << "Warning: no output defined!" << endl;
    for (int i =  0; i < n_out; i++) {
      Out.push_back(new Output(In, Lat, Seg, Mol, Sys, New, In[0]->OutputList[i], i, n_out));
       if (!Out[i]->CheckInput(start)) {
        cout << "input_error in output " << endl;
        success=false;
      }
     }
     MonteCarlo();
   }
  return success;
}


bool Teng::CP(transfer tofrom) {
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
			}

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
				if (Seg[clamp_seg]->bx[i]<1) {Seg[clamp_seg]->bx[i] +=MX; Seg[clamp_seg]->px1[i] +=MX; Seg[clamp_seg]->px2[i] +=MX;} 
				if (Seg[clamp_seg]->by[i]<1) {Seg[clamp_seg]->by[i] +=MY; Seg[clamp_seg]->py1[i] +=MY; Seg[clamp_seg]->py2[i] +=MY;} 
				if (Seg[clamp_seg]->bz[i]<1) {Seg[clamp_seg]->bz[i] +=MZ; Seg[clamp_seg]->pz1[i] +=MZ; Seg[clamp_seg]->pz2[i] +=MZ;} 
				Seg[clamp_seg]->H_MASK[((Seg[clamp_seg]->px1[i]-1)%MX+1)*JX + ((Seg[clamp_seg]->py1[i]-1)%MY+1)*JY + (Seg[clamp_seg]->pz1[i]-1)%MZ+1]=1;
				Seg[clamp_seg]->H_MASK[((Seg[clamp_seg]->px2[i]-1)%MX+1)*JX + ((Seg[clamp_seg]->py2[i]-1)%MY+1)*JY + (Seg[clamp_seg]->pz2[i]-1)%MZ+1]=1;
			}
		break;
		default:
			success=false; 
			cout <<"error in tranfer" << endl; 
		break;
	}
	return success; 
}

void Teng::WriteOutput(int subloop){
      PushOutput();
      Sys[0]->PushOutput(); // needs to be after pushing output for seg.
      Lat[0]->PushOutput();
      New[0]->PushOutput();
      int length = In[0]->MonList.size();
      for (int i = 0; i < length; i++)
        Seg[i]->PushOutput();
      length = In[0]->MolList.size();
      for (int i = 0; i < length; i++) {
        int length_al = Mol[i]->MolAlList.size();
        for (int k = 0; k < length_al; k++) {
          Mol[i]->Al[k]->PushOutput();
        } 
        Mol[i]->PushOutput();
       }
      // length = In[0]->AliasList.size();
      // for (int i=0; i<length; i++) Al[i]->PushOutput();

      for (int i = 0; i < n_out; i++) {
        Out[i]->WriteOutput(subloop);
	}
}

bool Teng::MonteCarlo() {
  if (debug) cout << "Monte Carlo in Teng" << endl;
	bool success; 
	t=0;
 	success=CP(to_cleng);
	New[0]->Solve(true);
	WriteOutput(t);
	n_p = X.size();
	X[0]=X[0]++; Y[0]=Y[0]++; 
	success=CP(to_segment);
	t++;
	New[0]->Solve(true);
	WriteOutput(t);

	return success;
}

void Teng::PutParameter(string new_param) {
  KEYS.push_back(new_param);
}
string Teng::GetValue(string parameter) {
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

void Teng::PushOutput() {
	int* point;
	for (int i = 0; i < n_out; i++) {
		Out[i]->PointerVectorInt.clear();
		Out[i]->PointerVectorReal.clear();
		Out[i]->SizeVectorInt.clear();
		Out[i]->SizeVectorReal.clear();
  		Out[i]->strings.clear();
  		Out[i]->strings_value.clear();
 		Out[i]->bools.clear();
 	 	Out[i]->bools_value.clear();
 		Out[i]->Reals.clear();
 	 	Out[i]->Reals_value.clear();
  		Out[i]->ints.clear();
  		Out[i]->ints_value.clear();
		if (Out[i]->name=="ana" || Out[i]->name=="kal") Out[i]->push("t",t);
		if (Out[i]->name=="ana" || Out[i]->name=="kal") Out[i]->push("MCS",MCS);
		if (Out[i]->name=="ana" || Out[i]->name=="vec") { //example for putting an array of Reals of arbitrary length to output
			string s="vector;0"; //the keyword 'vector' is used for Reals; the value 0 is the first vector, use 1 for the next etc, 
			Out[i]->push("gn",s); //"gn" is the name that will appear in output file
			Out[i]->PointerVectorReal.push_back(Mol[clp_mol]->gn); //this is the pointer to the start of the 'vector' that is reported to output.
			Out[i]->SizeVectorReal.push_back(sizeof(Mol[clp_mol]->gn)); //this is the size of the 'vector' that is reported to output
		}
		if (Out[i]->name=="ana" || Out[i]->name=="pos") { //example for putting 3 array's of Integers of arbitrary length to output
			string s="array;0";
			Out[i]->push("X",s);
			point=X.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(X.size());
			s="array;1";
			Out[i]->push("Y",s);
			point = Y.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(Y.size());
			s="array;2";
			Out[i]->push("Z",s);
			point=Z.data();
			Out[i]->PointerVectorInt.push_back(point);
			Out[i]->SizeVectorInt.push_back(Z.size());
		}
	
	}
}

