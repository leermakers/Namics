#include "bate.h"
#include "tools.h"

Balancedtensionless::Balancedtensionless(vector<Input *> In_, vector<Output *> Out_, vector<Lattice *> Lat_, vector<Segment *> Seg_, vector<State *> Sta_, vector<Reaction *> Rea_, vector<Molecule *> Mol_, vector<System *> Sys_,vector<Solve_scf*> New_,vector<Variate*> Var_,string name_):
	name{name_}, In{In_}, Out{Out_}, Lat{Lat_}, Mol{Mol_}, Seg{Seg_}, Sta{Sta_}, Rea{Rea_}, Sys{Sys_}, New{New_}, Var{Var_}
{
	if (debug)cout << "Constructor for Bate " << endl;
	KEYS.push_back("j_tolerance");
	KEYS.push_back("j_info");
	KEYS.push_back("monT");
	KEYS.push_back("monH");
	KEYS.push_back("monW");
	KEYS.push_back("var_tuning");
	KEYS.push_back("water");
  	KEYS.push_back("surfactant");
  	KEYS.push_back("chiTH=chiTW");
  	KEYS.push_back("kJ0m");
	//KEYS.push_back("control_parameter");
	KEYS.push_back("previous_guess");
	KEYS.push_back("caution_factor");

}
Balancedtensionless::~Balancedtensionless(){}

bool Balancedtensionless::CheckInput(int start_){
	if (debug) cout << "CheckInput for system " << endl;
	start=start_;
	previous_guess=true;
	bool success = true;
	string control;
	surfactant=-1;
	water=-1;
	string molname;
	control="wrong";
	HequalW=false;
	//int M=Lat[0]->M;
	int length = In[0]->MolList.size();
	success = In[0]->CheckParameters("bate", name, start, KEYS, PARAMETERS, VALUES);
	if (success) {

		caution_factor=1.0;
		if (GetValue("caution_factor").size()>0){
			caution_factor=In[0]->Get_Real(GetValue("caution_factor"),caution_factor);
			//if (caution_factor < 0) caution_factor=1;
		}
		J0m=0;
		if (GetValue("kJ0m").size()>0){
			J0m=In[0]->Get_Real(GetValue("kJ0m"),J0m);
			if (abs(J0m)>1) {
				cout <<"Kappa (bending modulus of monolayer) is positive and order unity; Spontaneous curvature monolayer should be a small number. Abs value of kJ0m is currently expected less than unity" << endl;
				success=false;
			}
		}

		previous_guess=true;
		if (GetValue("previous_guess").size()>0) {
			previous_guess=In[0]->Get_bool(GetValue("previous_guess"),true);
		}


		if (GetValue("water").size()>0) {
			molname=GetValue("water");
			for (int i=0; i<length; i++) {
				if (In[0]->MolList[i]==molname) water=i;
			}
		}
		if (water<0) {
			success=false; cout <<"In Bates you need to specify the water: use : Bate : 'name' : water : 'molname' " << endl;
		}

		if (GetValue("surfactant").size()>0) {
			molname=GetValue("surfactant");
			for (int i=0; i<length; i++) {
				if (In[0]->MolList[i]==molname) surfactant=i;
			}
		}
		if (surfactant<0) {
			success=false; cout <<"In Bates you need to specify the surfactant: use : Bate : 'name' : surfactant : 'valid molname' " << endl;
		}

		monT=-1;
		if (GetValue("monT").size()>0) {
			string MONT=GetValue("monT");
			int length =In[0]->MonList.size();
			for (int i=0; i<length; i++) {
				if (In[0]->MonList[i]==MONT) monT=i;
			}
		}
		if (monT<0) {
			success=false; cout <<"In Bates you need to specify the monT.You need a statement like  Bate : 'name' : monT : 'valid monname' " << endl;
		}

		monH=-1;
		if (GetValue("monH").size()>0) {
			string MONH=GetValue("monH");
			int length =In[0]->MonList.size();
			for (int i=0; i<length; i++) {
				if (In[0]->MonList[i]==MONH) monH=i;
			}
		}
		if (monH<0) {
			success=false; cout <<"In Bates you need to specify the monH. You need a statement like Bate : 'name' : monH : 'valid monname' " << endl;
		}

		monW=-1;
		if (GetValue("monW").size()>0) {
			string MONW=GetValue("monW");
			int length =In[0]->MonList.size();
			for (int i=0; i<length; i++) {
				if (In[0]->MonList[i]==MONW) monW=i;
			}
		}
		if (monW<0) {
			success=false; cout <<"In Bates you need to specify the monW. You need a statement like Bate : 'name' : monW : 'valid monname' " << endl;
		}

		if (Seg[monT]->chi[monW]==Seg[monT]->chi[monH]){
			if (GetValue("chiTH=chiTW").size()>0) {
				HequalW=In[0]->Get_bool(GetValue("chiTH=chiTW"),false);
			}
		} else {
			if (GetValue("chiTH=chiTW").size()>0) {
				HequalW=In[0]->Get_bool(GetValue("chiTH=chiTW"),false);
				if (HequalW) {
					success=false; cout <<" chiTH=chiTW is requested but the initial values for chiTH not equal to chiTW. Do not know what to do. Equalize chi's first " << endl;
				}
			}
		}


/*
		if (GetValue("control_parameter").size()>0) {
			vector<string> options;
			options.push_back("chi_monT_monX");
			if (!In[0]->Get_string(GetValue("control_parameter"),control,options,"In Bate, the control_parameter has not found proper option. "))
					success=false;
			else {
				monT=-1; monW=-1;
				if (control=="chi_monT_monX") {
					if (GetValue("monT").size()>0 && GetValue("monW").size()>0){
						string MONT=GetValue("monT");
						string MONW=GetValue("monW");
						int length =In[0]->MonList.size();
						monT=-1; monW=-1;
						for (int i=0; i<length; i++) {
							if (Seg[i]->name == MONT) monT=i;
							if (Seg[i]->name == MONW) monW=i;
						}
					}
					if (monT>-1 && monW>-1) {
						if (monT==monW){
							success=false; cout <<"monT and monW can not be the same segment type in Bate. control parameter is not identified correctly"<<endl;
						}
					} else {
						success=false; cout <<"monT and/or monW not found. Unable to identify control_parameter properly in Bate" << endl;
					}
				}
			}
		} else {
			if (control == "wrong") {
				success=false; cout <<"In Bate control_parameter is unknown; not certain what to do now. .." <<endl;
			}
		}
*/


		j_tolerance = 1e-6;
		if (GetValue("j_tolerance").size()>0) {
			j_tolerance=In[0]->Get_Real(GetValue("j_tolerance"),j_tolerance);
			if (j_tolerance>1e-2 || j_tolerance <0) {
				success=false ; cout <<"j_tolerance not in proper range in Bate "<< endl;
			}
		}
		j_info=3;

		if (GetValue("j_info").size()>0) {
			j_info=In[0]->Get_int(GetValue("j_info"),j_info);
			if (j_info <0) {
				cout <<"j_info should be positive integer. default j_info = 3 is used" << endl;
				j_info=3;
			}
		}

		chi_start=-123;
		chi_step=-123;
		chi_end=-123;

		if (GetValue("var_tuning").size()>0) {
			string variate=GetValue("var_tuning");
			vector<string>sub;
			In[0]->split(variate,';',sub);
			if (sub.size()==3) {
				chi_start=In[0]->Get_Real(sub[0],-123);
				chi_step=In[0]->Get_Real(sub[1],-123);
				chi_end=In[0]->Get_Real(sub[2],-123);
				n_steps=0;
				if (chi_start !=-123 && chi_end!=-123) {
					if (chi_start > chi_end) {
						if (chi_step > 0 ) {success = false; cout << "expected chi_step <0 because chi_start > chi_end "<< endl; }
						n_steps=round(( chi_end-chi_start)/chi_step);
						//else cout<< "Number of chi-steps is " << (chi_end-chi_start)/chi_step << endl;
					}
					if (chi_end > chi_start) {
						if (chi_step < 0 ) {success = false; cout << "expected chi_step >0 because chi_end > chi_start "<< endl; }
						n_steps =round((chi_end-chi_start)/chi_step);
						//else cout<< "Number of chi-steps is " << (chi_start-chi_end)/chi_step << endl;
					}
					if (n_steps <0) {
						cout << "combination chi-start, chi-step, chi-end gives negavite number of steps....try again.." << endl;
						success =false;
					} else n_steps++;
				} else {
					success=false ; cout <<" expected chi_start > 0 and chi_end > 0" << endl;
				}
			} else {success = false; cout <<" expected three real value  separated by ';', in var_tuning " << endl; }
		}
	}

	ControlParameter=chi_C_D;
	return success;
}

string Balancedtensionless::GetValue(string parameter)
{
	if (debug) cout << "GetValue " + parameter + " for Bate " << endl;
	int length = PARAMETERS.size();
	for (int i = 0; i < length; ++i)
	{
		if (parameter == PARAMETERS[i])
		{
			//cout << "return " << VALUES[i] << endl;
			return VALUES[i];
		}
	}
	return "";
}


Real Balancedtensionless:: FixedPoint(Real Xs) { //theta surfactant
	if (debug) cout << "Fixed point in Bates " << endl;
	//int M=Lat[0]->M;
	//int fjc=Lat[0]->fjc;
	Sys[0]->MakeItsLists();
	New[0]->AllocateMemory();
	New[0]->Guess(X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);

	Mol[surfactant]->PutTheta(Xs);


	if (search_nr < 0 && ets_nr < 0 && etm_nr < 0) {
		if (debug) cout << "Fixed point to solve " << endl;
		New[0]->Solve(true);
	} else {
		if (debug) cout << "Solve to superiteration " << endl;
		//cout <<"Seg[monH]->chi[monW]" << Seg[monH]->chi[monW] << endl;
		//cout <<"Seg[monT]->chi[monW]" << Seg[monT]->chi[monW] << endl;
		New[0]->SuperIterate(search_nr, target_nr, ets_nr, etm_nr, bm_nr);
	}
	return Mol[surfactant]->theta;
}

bool Balancedtensionless::PutChi(Real CHI) {
	bool success=true;
	if (ControlParameter==chi_C_D) {
		Seg[monH]->chi[monW]=CHI;
		Seg[monW]->chi[monH]=CHI;
	}
	return success;
}


Real Balancedtensionless::zero_J0(Real guessXs, Real GuessChi) {
	if (debug) cout << "In Bate: zero_j0 " << endl;
    Real XS=FixedPoint(guessXs);
    Real cx=-0.01*caution_factor;
    Real xa=0;
    if (ControlParameter==chi_C_D) {
		PutChi(GuessChi);
		xa=GuessChi;
	}

    int j_calls=0;
    Real fxa=Sys[0]->GetSpontaneousCurvature(0.0)-J0m;
    if (abs(fxa) < j_tolerance){
        return XS;
	}
    Real xb=0;

	switch (ControlParameter) {
		case chi_C_D:
			if (fxa<0) {
				xb=xa+cx;
			} else {
				xb=xa-cx;
			}
			PutChi(xb);
			XS=FixedPoint(XS);
			break;
		default:
			cout<<"program error; report this message. " <<endl;
			break;
	}

    Real fxb=Sys[0]->GetSpontaneousCurvature(0.0)-J0m;
    while (fxa*fxb>0){
		j_calls++;
        xa=xb;
        fxa=fxb;
        switch(ControlParameter){
			case chi_C_D:
				//cout <<"not implemented yet" << endl;
				if (fxa<0) xb=xa+cx; else xb=xa-cx;
				PutChi(xb);
				XS=FixedPoint(XS-0.01);
				break;
			default:
				break;
		}

        fxb=Sys[0]->GetSpontaneousCurvature(0.0)-J0m;
        switch(ControlParameter) {

			case chi_C_D:
				if (j_calls%j_info==0) cout << "j_it = " << j_calls << " chi_"<<Seg[monH]->name<<"_"<<Seg[monW]->name << " = "  << xb << " kJ0 = " << fxb << endl;
				break;
			default:
				break;
		}
	}

    Real xc=(xa*fxb-xb*fxa)/(fxb-fxa);
    switch (ControlParameter) {
		case chi_C_D:
			PutChi(xc);
			XS=FixedPoint(XS-0.01);
			break;
		default:
		    XS=FixedPoint(XS-0.01);
			break;
	}
    Real fxc=Sys[0]->GetSpontaneousCurvature(0.0)-J0m;
    j_calls=0;
    while (abs(fxc)>j_tolerance) {
		j_calls++;
        if (fxa*fxc<0){
            xb=xc; fxb=fxc;
        } else {
            xa=xc; fxa=fxc;
		}
        xc=(xa*fxb-xb*fxa)/(fxb-fxa);
        switch (ControlParameter) {
			case chi_C_D:
				PutChi(xc);
				XS=FixedPoint(XS-0.01);
				break;
			default:
				XS=FixedPoint(XS-0.01);
				break;
		}

        fxc=Sys[0]->GetSpontaneousCurvature(0.0)-J0m;
        switch(ControlParameter) {
			case chi_C_D:
				if (j_calls%j_info==0 && j_calls>3*j_info) cout << "j_it = " << j_calls << " chi_"<<Seg[monH]->name<<"_"<<Seg[monW]->name << " = "  << xc << " kJ0 = " << fxc << endl;
				break;
			default :
				break;
		}

        if (j_calls%100==0 || j_calls%101==0) {
			cout <<"Restart J0 iteration" << endl;
			switch (ControlParameter) {
				case chi_C_D:
					PutChi(xc);
					return zero_J0(XS-0.01,xc);
					break;
				default:
					return zero_J0(XS-0.01,xc);
					break;
			}
		}
	}
	switch(ControlParameter) {
		case chi_C_D:
			return xc;
			break;
		default:
		    return xc;
			break;
	}
}


bool Balancedtensionless:: WriteResults() {
	if (debug) cout << "In Bate writeResults " << endl;
	int n_out = In[0]->OutputList.size();
	for (int ii = 0; ii < n_out; ii++) {
		Out.push_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[ii], ii, n_out));
		if (!Out[ii]->CheckInput(start)){
			cout << "input_error in output " << endl;
			return 0;
		} else {
			if (Out[ii]->name=="kal") { //this is to make sure that append is set to true when 'kal'-file is not initiated for the first time.
				if (kal_append) Out[ii]->append=true;
				else {kal_append=true;}
			}
		}
	}
	New[0]->PushOutput();


	for (int ii = 0; ii < n_out; ii++){
		Out[ii]->WriteOutput(subloop);
	}
	if (Sys[0]->final_guess == "file"){ //if iv have changed: see namics variant.
		Lat[0]->StoreGuess(Sys[0]->guess_outputfile, New[0]->xx , METHOD, MONLIST, STATELIST, CHARGED, start);
	}
	subloop++;
	return true;
}


bool Balancedtensionless::Doit(Real* X_,string METHOD_,vector<string> MONLIST_,vector<string> STATELIST_,bool CHARGED_,int MX_,int MY_,int MZ_,int fjc_old_,int search_nr_,int ets_nr_,int etm_nr_,int target_nr_,int bm_nr_,int subloop_, bool& kal_append_) {
	X=X_; METHOD=METHOD_; MONLIST=MONLIST_; STATELIST=STATELIST_; CHARGED=CHARGED_; MX=MX_; MY=MY_; MZ=MZ_; fjc_old=fjc_old_; search_nr=search_nr_; ets_nr=ets_nr_;  etm_nr=etm_nr_; target_nr=target_nr_; bm_nr=bm_nr_; subloop=subloop_; kal_append=kal_append_;
	if (debug) cout <<"In Bate: I'll do it " << endl;
	bool success=true;
	Real CHI;

	if (n_steps>0) {
		for (int i=0; i< n_steps; i++) {
			CHI=chi_start+i*chi_step;
			if (HequalW) {
				Seg[monT]->chi[monW]=CHI;
				Seg[monW]->chi[monT]=CHI;
				Seg[monT]->chi[monH]=CHI;
				Seg[monH]->chi[monT]=CHI;
			} else {
				Seg[monT]->chi[monW]=CHI;
				Seg[monW]->chi[monT]=CHI;
			}
			cout << "BaTe-problem " << i+1 << " of " << n_steps << " chi = " << CHI << endl;

			switch (ControlParameter) {
				case chi_C_D:
					if (i==0){
						ini_surf=Mol[surfactant]->theta;
						ini_chi=Seg[monH]->chi[monW];
					} else {
						if (!previous_guess) {
							Mol[surfactant]->theta=ini_surf;
							Seg[monH]->chi[monW]=ini_chi;
							Seg[monW]->chi[monH]=ini_chi;
						}
					}
					zero_J0(Mol[surfactant]->theta,Seg[monH]->chi[monW]);
					break;
				default:
					break;
			}

			WriteResults();
		}
	} else {
		switch (ControlParameter) {
			case chi_C_D:
				zero_J0(Mol[surfactant]->theta,Seg[monH]->chi[monW]);
				break;
			default:
				break;
		}

		WriteResults();
	}
	kal_append_=true;
	return success;
}

