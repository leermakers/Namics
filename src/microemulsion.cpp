#include "microemulsion.h"
#include "tools.h"

Microemulsion::Microemulsion(vector<Input *> In_, vector<Output *> Out_, vector<Lattice *> Lat_, vector<Segment *> Seg_, vector<State *> Sta_, vector<Reaction *> Rea_, vector<Molecule *> Mol_, vector<System *> Sys_,vector<Solve_scf*> New_,vector<Variate*> Var_,string name_):
	name{name_}, In{In_}, Out{Out_}, Lat{Lat_}, Mol{Mol_}, Seg{Seg_}, Sta{Sta_}, Rea{Rea_}, Sys{Sys_}, New{New_}, Var{Var_}
{
	if (debug)cout << "Constructor for microemulsion " << endl;
  	KEYS.push_back("surfactant");
	KEYS.push_back("co_solvent");
	KEYS.push_back("g_tolerance");
	KEYS.push_back("j_tolerance");
	KEYS.push_back("g_info");
	KEYS.push_back("j_info");
	KEYS.push_back("monA");
	KEYS.push_back("monB");
	KEYS.push_back("var_chi_monA_monB");
}
Microemulsion::~Microemulsion(){}

bool Microemulsion::CheckInput(int start_){
	if (debug) cout << "CheckInput for system " << endl;
	start=start_;
	bool success = true;
	surfactant=-1;
	co_solvent=-1;
	string molname;
	//int M=Lat[0]->M;
	int length = In[0]->MolList.size();
	success = In[0]->CheckParameters("micro", name, start, KEYS, PARAMETERS, VALUES);
	if (success) {
		if (GetValue("surfactant").size()>0) {
			molname=GetValue("surfactant");
			for (int i=0; i<length; i++) {
				if (In[0]->MolList[i]==molname) surfactant=i;
			}
			if (surfactant<0) {
				success=false; cout <<"In microemulsions you need to specify the surfactant: use : micro : 'name' : surfactant : 'molname' " << endl;
			}
		}
		if (GetValue("co_solvent").size()>0) {
			molname=GetValue("co_solvent");
			for (int i=0; i<length; i++) {
				if (In[0]->MolList[i]==molname) co_solvent=i;
			}
			if (co_solvent<0) {
				success=false; cout <<"In microemulsions you need to specify the co_solvent: use : micro : 'name' : co_solvent : 'molname' " << endl;
			}
		}
		g_tolerance = 1e-8;
		if (GetValue("g_tolerance").size()>0) {
			g_tolerance=In[0]->Get_Real(GetValue("g_tolerance"),g_tolerance);
			if (g_tolerance>1e-5 || g_tolerance <0) {
				success=false ; cout <<"g_tolerance not in proper range in microemulsion "<< endl;
			}
		}
		j_tolerance = 1e-6;
		if (GetValue("j_tolerance").size()>0) {
			j_tolerance=In[0]->Get_Real(GetValue("j_tolerance"),j_tolerance);
			if (j_tolerance>1e-2 || j_tolerance <0) {
				success=false ; cout <<"j_tolerance not in proper range in microemulsion "<< endl;
			}
		}
		g_info=5; j_info=3;
		if (GetValue("g_info").size()>0) {
			g_info=In[0]->Get_int(GetValue("j_tolerance"),g_info);
			if (g_info <0) {
				cout <<"g_info should be positive integer. default g_info = 5 is used" << endl;
				g_info=5;
			}
		}
		if (GetValue("j_info").size()>0) {
			j_info=In[0]->Get_int(GetValue("j_tolerance"),j_info);
			if (j_info <0) {
				cout <<"j_info should be positive integer. default j_info = 3 is used" << endl;
				j_info=3;
			}
		}


		length =In[0]->MonList.size();
		monA=-1; monB=-1;
		chi_start=-123;
		chi_step=-123;
		chi_end=-123;
		if (GetValue("monA").size()>0 && GetValue("monB").size()>0){
			string MONA=GetValue("monA");
			string MONB=GetValue("monB");
			for (int i=0; i<length; i++) {
				if (Seg[i]->name == MONA) monA=i;
				if (Seg[i]->name == MONB) monB=i;
			}
			if (monA!=monB && monA>-1 && monB>-1) {
				if (GetValue("var_chi_monA_monB").size()>0) {
					string variate=GetValue("var_chi_monA_monB");
					vector<string>sub;
					In[0]->split(variate,';',sub);
					if (sub.size()==3) {
						chi_start=In[0]->Get_Real(sub[0],-123);
						chi_step=In[0]->Get_Real(sub[1],-123);
						chi_end=In[0]->Get_Real(sub[2],-123);
						n_steps=0;

						if (chi_start !=-123 && chi_end!=-123) {

							if (chi_start > chi_end) {
								if (chi_step > 0 ) {success = false; cout << "expected chi_step <0 becaue ghi_start > chi_end "<< endl; }
								n_steps=round(( chi_end-chi_start)/chi_step);
								//else cout<< "Number of chi-steps is " << (chi_end-chi_start)/chi_step << endl;
							}
							if (chi_end > chi_start) {
								if (chi_step < 0 ) {success = false; cout << "expected chi_step >0 becaue ghi_end > chi_start "<< endl; }
								n_steps =round((chi_start-chi_end)/chi_step);
								//else cout<< "Number of chi-steps is " << (chi_start-chi_end)/chi_step << endl;
							}
							if (n_steps <0) {
								cout << "combination chi-start, chi-step, chi-end gives negivate number of steps....try again.." << endl;
								success =false;
							} else n_steps++;
						} else {
							success=false ; cout <<" expected chi_start > 0 and chi_end > 0" << endl;
						}
					} else {success = false; cout <<" expected three real value  separated by ';', in var_chi_monA_monB " << endl; }
				} else {
					success=false; cout<<" var_chi_monA_monB input is missing, while monA and monB were specified; expected value: chi_start;chi_step;chi_end " << endl;
				}
			}
		}
	}
	return success;
}

string Microemulsion::GetValue(string parameter)
{
	if (debug) cout << "GetValue " + parameter + " for microemulsion " << endl;
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

bool Microemulsion:: FixedPoint(Real Ts, Real Tc) {
	if (debug) cout << "Fixed point in microemulsions " << endl;
	Sys[0]->MakeItsLists();
	New[0]->AllocateMemory();
	New[0]->Guess(X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);
	Mol[surfactant]->PutTheta(Ts); Mol[co_solvent]-> PutTheta(Tc);
	if (search_nr < 0 && ets_nr < 0 && etm_nr < 0) {
		if (debug) cout << "Fixed point to solve " << endl;
		New[0]->Solve(true);
	} else {
		if (debug) cout << "Solve to superiteration " << endl;
		New[0]->SuperIterate(search_nr, target_nr, ets_nr, etm_nr, bm_nr);
	}
	return true;
}



bool Microemulsion::Doit(Real* X_,string METHOD_,vector<string> MONLIST_,vector<string> STATELIST_,bool CHARGED_,int MX_,int MY_,int MZ_,int fjc_old_,int search_nr_,int ets_nr_,int etm_nr_,int target_nr_,int bm_nr_,int subloop_, bool kal_append_) {
	X=X_; METHOD=METHOD_; MONLIST=MONLIST_; STATELIST=STATELIST_; CHARGED=CHARGED_; MX=MX_; MY=MY_; MZ=MZ_; fjc_old=fjc_old_; search_nr=search_nr_; ets_nr=ets_nr_;  etm_nr=etm_nr_; target_nr=target_nr_; bm_nr=bm_nr_; subloop=subloop_; kal_append=kal_append_;
	if (debug) cout <<"In microemulsion: I'll do it " << endl;
	Real CHI;
	Real result=10;
	if (n_steps>0) {
		for (int i=0; i< n_steps; i++) {
			CHI=chi_start+i*chi_step;
			Seg[monA]->chi[monB]=CHI;
			Seg[monB]->chi[monA]=CHI;
			cout << "Micro-problem " << i+1 << " of " << n_steps << " chi = " << CHI << endl;
			if (result>0) {
				result=zero_J0(Mol[surfactant]->theta,Mol[co_solvent]->theta);
				WriteResults();
			}
		}
	} else {
		zero_J0(Mol[surfactant]->theta,Mol[co_solvent]->theta);
		WriteResults();
	}

	return true;
}

Real Microemulsion::get_gamma(Real Ts, Real Tc){
	if (debug) cout <<"In microemulsion: get_gamma" << endl;
    FixedPoint(Ts,Tc);
    Real gamma=Sys[0]->GetGrandPotential();
    return gamma;
}

Real Microemulsion::zero_gamma(Real guess, Real Tc) {
    if (debug) cout <<"In Microemulsion: zero_gamma. Gusess = " << guess << " Tc = " << Tc << endl;
	int sign=0;
	int g_calls=0;

    Real dx=0.001;
    Real x1=guess;
    Real x2=guess+dx;
    Real x3=guess+2*dx;
    Real fx1=get_gamma(x1,Tc); g_calls++;
    Real fx2=get_gamma(x2,Tc); g_calls++;
    Real fx3=get_gamma(x3,Tc); g_calls++;
    Real gradient=(fx3-fx1)/(2*dx);
    //if (gradient>0) {
    //    cout << "Error in zero_gamma. Gradient is positive" << endl;
	//}
    Real step=-1.0*fx2/gradient;
    //Real hessian=(fx1-2*fx2+fx3)/(dx*dx);
    if (step>guess){
        cout << "Warning in zero_gamma? Step is very large... proceeding in careful mode" << endl;
        step=guess/2;
	}
    Real x4=x2+step*9.5/10.0;
    Real fx4=get_gamma(x4,Tc); g_calls++;
    while (fx4>0) {
        //cout << "searching gamma < 0 in zero_gamma"<< endl;
        x2=x4;
        if (fx2<fx4){
            cout << "walking the wrong way... possible no zero for gamma?" << endl;
            fx4=-1;
        } else {
            fx2=fx4;
            x4=x4+step/10.0;
            fx4=get_gamma(x4,Tc); g_calls++;
		}
	}
    Real xa=x2;
    Real fxa=fx2;
    Real xb=x4;
    Real fxb=fx4;
    //cout << "here we go...for gamma..." << endl;
    Real xc=(xa*fxb-xb*fxa)/(fxb-fxa);
    Real fxc=get_gamma(xc,Tc); g_calls=0;
    while (abs(fxc)>g_tolerance) {
        if (fxa*fxc<0){
            xb=xc;
        } else {
            xa=xc;
		}
        xc=(xa*fxb-xb*fxa)/(fxb-fxa);
        fxc=get_gamma(xc,Tc);g_calls++;
        if (fxc>0) {
            sign++;
            if (sign>5) {
                sign=0;
                xc=(2*xb+xc)/3;
                fxc=get_gamma(xc,Tc);g_calls++;
			}
        } else {
            sign -=1;
            if (sign<5){
                sign=0;
                xc=(2*xa+xc)/3;
                fxc=get_gamma(xc,Tc);g_calls++;
			}
		}
		if (g_calls%40==0 || g_calls%41==0) {
			cout <<"Restart gamma iteration" << endl;
			return zero_gamma(xc,Tc);
		}
		if (g_calls%g_info==0 && g_calls>4*g_info) cout<<"g_calls = " << g_calls << " theta surfactant =" << xc <<  " gamma = "<< fxc<< endl;
	}
    return xc;
}

Real Microemulsion::zero_J0(Real guessTs, Real GuessTc) {
	if (debug) cout << "In microemulsion: zero_j0 " << endl;
    Real ThetaS=zero_gamma(guessTs,GuessTc);
    Real dx=1;
    Real xa=GuessTc;
    int j_calls=0;
    Real fxa=Sys[0]->GetSpontaneousCurvature();
    if (abs(fxa) < j_tolerance){
        //cout << "kJ0=" << fxa << endl;
        return GuessTc;
	}

    Real xb;
    if (fxa<0) xb=xa-dx; else xb=xa+dx;
    if (xb>0) ThetaS=zero_gamma(ThetaS-0.1,xb); else {
		dx=dx/2;
		ThetaS=zero_gamma(ThetaS-0.1,dx);
	}
    Real fxb=Sys[0]->GetSpontaneousCurvature();
    while (fxa*fxb>0){
		j_calls++;
        xa=xb;
        fxa=fxb;
        if (fxa<0) xb=xa-dx; else xb=xa+dx;
        if (xb>0) ThetaS=zero_gamma(ThetaS-0.1,xb); else {
			dx=dx/2; xb=dx;
			ThetaS=zero_gamma(ThetaS-0.1,xb);
		}
        fxb=Sys[0]->GetSpontaneousCurvature();
        if (j_calls%j_info==0) cout << "j_calls = " << j_calls << " theta cosolvent = " << xb <<  " kJ0 = " << fxb << endl;
        if (xb<0.01) {
			cout <<"Amount of co-solvent tends to go negative" <<endl;
			return -1;
		}
	}
    //cout << "here we go...for J0..." << endl;
    Real xc=(xa*fxb-xb*fxa)/(fxb-fxa);
    if (xc>0) ThetaS=zero_gamma(ThetaS-0.1,xc); else ThetaS=zero_gamma(ThetaS-0.1,0.01);
    Real fxc=Sys[0]->GetSpontaneousCurvature();
    j_calls=0;
    while (abs(fxc)>j_tolerance) {
		j_calls++;
        if (fxa*fxc<0){
            xb=xc;
        } else {
            xa=xc;
		}
        xc=(xa*fxb-xb*fxa)/(fxb-fxa);
        if (xc>0) ThetaS=zero_gamma(ThetaS-0.1,xc); else {
			cout <<" cosolvent tend to go negative. Problem terminated " << endl;
			return xc;
		}
        fxc=Sys[0]->GetSpontaneousCurvature();
        if (j_calls%j_info==0 && j_calls>3*j_info) cout << "j_calls = " << j_calls << " theta cosolvent = " << xc << " kJ0 = " << fxc << endl;
        if (j_calls%20==0 || j_calls%21==0) {
			cout <<"Restart J0 iteration" << endl;
			return zero_gamma(ThetaS,xc);
		}
	}
    return xc;
}

bool Microemulsion:: WriteResults() {
	if (debug) cout << "In microemulsions writeResults " << endl;
	int n_out = In[0]->OutputList.size();

	for (int ii = 0; ii < n_out; ii++) {
		Out.push_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[ii], ii, n_out));
		if (!Out[ii]->CheckInput(start)){
			cout << "input_error in output " << endl;
			return 0;
		} else {
			if (Out[ii]->name=="kal") { //this is to make sure that append is set to true when 'kal'-file is not initiated for the first time.
				if (kal_append) Out[ii]->append=true;
				else kal_append=true;
			}
		}
	}
	New[0]->PushOutput();
	cout <<"theta surfactant = " << Mol[surfactant]->theta << " theta co-solvent = " << Mol[co_solvent]->theta << " kBar = " << Sys[0]->GetKBar() << endl;
	for (int ii = 0; ii < n_out; ii++){
		Out[ii]->WriteOutput(subloop);
	}
	if (Sys[0]->final_guess == "file"){ //if iv have changed: see namics variant.
		Lat[0]->StoreGuess(Sys[0]->guess_outputfile, New[0]->xx , METHOD, MONLIST, STATELIST, CHARGED, start);
	}
	subloop++;
	return true;
}

